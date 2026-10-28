#include "utils/incident_recorder.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <vector>
#include <fcntl.h>
#include <sys/file.h>
#include <sys/resource.h>
#include <sys/syscall.h>
#include <unistd.h>
#include <glog/logging.h>
#include <rosbag2_cpp/writer.hpp>

namespace lightning::diagnostics {
namespace {
namespace fs = std::filesystem;
constexpr std::uint64_t MiB = 1024 * 1024;
constexpr char Owner[] = "LIGHTNING_LM_INCIDENT_V1";
std::int64_t Steady() {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now().time_since_epoch()).count();
}
std::int64_t Wall() {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::system_clock::now().time_since_epoch()).count();
}
double Age(std::int64_t now, std::int64_t then) { return detail::IncidentAgeSeconds(now, then); }
std::string Json(const std::string& input) {
    std::ostringstream out;
    out << '"';
    for (unsigned char c : input) {
        if (c == '"' || c == '\\') out << '\\' << c;
        else if (c < 32 || c == '<' || c == '>' || c == '&')
            out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << int(c) << std::dec;
        else out << c;
    }
    return out.str() + '"';
}
std::string Number(double value) {
    if (!std::isfinite(value)) return "null";
    std::ostringstream out; out << std::setprecision(10) << value; return out.str();
}
void AtomicWrite(const fs::path& path, const std::string& contents) {
    const auto temporary = path.string() + ".tmp";
    { std::ofstream file(temporary, std::ios::trunc); file << contents; file.flush();
      if (!file) throw std::runtime_error("cannot write " + temporary); }
    fs::rename(temporary, path);
}
std::uint64_t DirectoryBytes(const fs::path& path) {
    std::uint64_t result = 0;
    for (const auto& entry : fs::recursive_directory_iterator(path))
        if (!entry.is_symlink() && entry.is_regular_file()) result += entry.file_size();
    return result;
}
void CopySmallFile(const fs::path& source, const fs::path& target) {
    if (!source.empty() && fs::is_regular_file(source) && fs::file_size(source) <= MiB)
        fs::copy_file(source, target, fs::copy_options::overwrite_existing);
}
enum Reason : unsigned {
    Raw2 = 1, Raw4 = 2, Raw9 = 4, Startup = 8, Output = 16, Input = 32, FrozenStamp = 64
};
unsigned FaultBit(int code) { return code == 2 ? unsigned(Raw2) : code == 4 ? unsigned(Raw4) : code == -9 ? unsigned(Raw9) : 0u; }
std::string Reasons(unsigned bits) {
    std::string result = "[";
    for (const auto& item : std::vector<std::pair<unsigned, std::string>>{
             {Raw2, "raw_2"}, {Raw4, "raw_4"}, {Raw9, "raw_-9"},
             {Startup, "startup_without_output"}, {Output, "publish_completion_timeout"},
             {Input, "input_callback_silence"}, {FrozenStamp, "input_stamp_not_advancing"}}) {
        if (bits & item.first) { if (result.size() > 1) result += ','; result += Json(item.second); }
    }
    return result + ']';
}
}  // namespace

IncidentRecorder::Options IncidentRecorder::Options::FromYaml(const YAML::Node& n) {
    Options o;
    if (!n) return o;
#define READ_OPTION(name) if (n[#name]) o.name = n[#name].as<decltype(o.name)>()
    READ_OPTION(enabled); READ_OPTION(pre_seconds); READ_OPTION(post_seconds);
    READ_OPTION(max_event_seconds); READ_OPTION(startup_timeout_seconds);
    READ_OPTION(output_timeout_seconds); READ_OPTION(input_timeout_seconds);
    READ_OPTION(cooldown_seconds); READ_OPTION(write_mb_per_second);
    READ_OPTION(memory_mb); READ_OPTION(max_message_mb); READ_OPTION(event_mb);
    READ_OPTION(disk_mb); READ_OPTION(min_free_mb); READ_OPTION(segment_mb); READ_OPTION(max_events);
#undef READ_OPTION
    o.Validate();
    return o;
}
void IncidentRecorder::Options::Validate() const {
    for (double value : {post_seconds, max_event_seconds, startup_timeout_seconds,
                         output_timeout_seconds, input_timeout_seconds, write_mb_per_second})
        if (!std::isfinite(value) || value <= 0 || value > 3600)
            throw std::invalid_argument("incident_capture: positive finite timing/rate values required");
    if (!std::isfinite(pre_seconds) || pre_seconds < 0 || pre_seconds > 120 ||
        !std::isfinite(cooldown_seconds) || cooldown_seconds < 0 || cooldown_seconds > 3600 ||
        max_event_seconds < post_seconds || max_event_seconds > 120 ||
        memory_mb < 1 || memory_mb > 65536 || max_message_mb < 1 || max_message_mb > memory_mb ||
        event_mb < 1 || event_mb > 65536 || disk_mb < event_mb + 1 || disk_mb > 1048576 ||
        segment_mb < 1 || segment_mb > event_mb || min_free_mb > 1048576 ||
        max_events < 1 || max_events > 1000)
        throw std::invalid_argument("incident_capture: inconsistent or unbounded resource limits");
}

class IncidentRecorder::Impl {
 public:
    struct Stream {
        std::string topic, type;
        Kind kind;
        Serializer serialize;
        std::atomic<std::int64_t> received_at{0}, advanced_at{0}, sensor_stamp{0};
        std::atomic<std::uint64_t> count{0};
    };
    struct Sample {
        std::int64_t steady = 0, wall = 0;
        Safety safety;
        Pipeline pipeline;
        double imu_age = -1, lidar_age = -1, stamp_age = -1;
        double prepared_age = -1, publish_begin_age = -1, publish_end_age = -1, pipeline_age = -1;
        bool publishing = false;
        std::uint64_t dropped = 0;
    };
    struct Budget { std::atomic<std::uint64_t> held{0}, peak{0}; };
    struct Record {
        std::shared_ptr<Budget> budget;
        std::shared_ptr<const void> message;
        std::uint64_t cost = 0;
        std::int64_t steady = 0, wall = 0;
        int stream = -1;  // -1 is a lightweight timeline sample.
        ~Record() { if (budget) budget->held.fetch_sub(cost, std::memory_order_relaxed); }
    };
    using Row = std::shared_ptr<Record>;
    struct Event {
        std::int64_t trigger = 0, wall = 0, until = 0;
        std::atomic<unsigned> reasons{0};
        double pre_coverage = 0;
        std::deque<Row> pending;
        bool closed = false;
        std::atomic<bool> shortened{false};
        fs::path path;
        std::vector<Sample> samples;
        std::vector<std::string> segments;
        std::vector<std::uint64_t> counts;
        std::vector<std::int64_t> first_written, last_written;
        std::vector<double> trigger_callback_age, trigger_stamp_age;
        std::uint64_t bytes = 0, omitted = 0;
        std::uint64_t closed_disk_bytes = 0;
        std::string error;
    };
    Options options;
    fs::path directory, config;
    std::vector<std::shared_ptr<Stream>> streams;
    std::shared_ptr<Budget> budget = std::make_shared<Budget>();
    std::atomic<std::uint64_t> received{0}, admitted{0}, busy{0}, full{0}, oversized{0};
    std::atomic<std::uint64_t> events{0}, merged{0}, suppressed{0}, errors{0}, written{0};
    std::atomic<std::int64_t> prepared{0}, publish_begin{0}, publish_end{0}, pipeline_at{0};
    std::atomic<bool> publishing{false}, running{false}, stopping{false}, stop_monitor{false};
    std::atomic<unsigned> fault_edges{0};
    std::atomic<int> last_raw{0};
    std::mutex snapshot_mutex, mutex;
    Safety latest_safety;
    Pipeline latest_pipeline;
    std::deque<Row> history;
    std::shared_ptr<Event> active;
    std::condition_variable wake;
    std::thread monitor, writer;
    std::int64_t started = 0, last_finished = 0;
    int lock_fd = -1;

    Impl(Options o, std::string d, std::string c)
        : options(o), directory(std::move(d)), config(std::move(c)) { options.Validate(); }
    ~Impl() { Stop(); if (lock_fd >= 0) close(lock_fd); }

    Stats StatsNow() const {
        return {budget->held.load(), budget->peak.load(), received.load(), admitted.load(),
                busy.load(), full.load(), oversized.load(), events.load(), merged.load(),
                suppressed.load(), errors.load(), written.load()};
    }
    bool Add(int stream, std::shared_ptr<const void> message, std::uint64_t bytes,
             std::int64_t steady, std::int64_t wall) {
        ++received;
        if (bytes > options.max_message_mb * MiB) { ++oversized; return false; }
        // Never evict/destroy large messages on an algorithm callback.
        std::unique_lock<std::mutex> lock(mutex, std::try_to_lock);
        if (!lock.owns_lock()) { ++busy; return false; }
        if (stop_monitor) return false;
        bytes += sizeof(Record) + 64;
        if (budget->held + bytes > options.memory_mb * MiB) { ++full; return false; }
        try {
            auto row = std::make_shared<Record>();
            row->message = std::move(message); row->stream = stream;
            row->steady = steady; row->wall = wall;
            row->cost = bytes; row->budget = budget;
            const auto held = budget->held.fetch_add(bytes) + bytes;
            budget->peak = std::max(budget->peak.load(), held);
            history.push_back(row);
            if (active && !active->closed) { active->pending.push_back(row); wake.notify_one(); }
            ++admitted;
            return true;
        } catch (...) { ++full; return false; }
    }
    bool Start() {
        if (!options.enabled || running || stopping) return false;
        try {
            if (directory.empty() || !directory.is_absolute())
                throw std::runtime_error("an absolute incident output directory is required");
            fs::create_directories(directory);
            directory = fs::canonical(directory);
            lock_fd = open((directory / ".capture.lock").c_str(), O_CREAT | O_RDWR | O_CLOEXEC, 0600);
            if (lock_fd < 0 || flock(lock_fd, LOCK_EX | LOCK_NB) != 0)
                throw std::runtime_error("incident directory is already in use or inaccessible");
            started = Steady(); running = true;
            writer = std::thread([this] {
                try { WriterLoop(); }
                catch (const std::exception& e) { ++errors; stop_monitor = true; LOG(ERROR) << "INCIDENT_CAPTURE writer stopped: " << e.what(); }
            });
            monitor = std::thread([this] {
                try { MonitorLoop(); }
                catch (const std::exception& e) {
                    ++errors; stop_monitor = true;
                    { std::lock_guard<std::mutex> lock(mutex);
                      if (active) { active->closed = true; active->shortened = true; } }
                    wake.notify_all();
                    LOG(ERROR) << "INCIDENT_CAPTURE monitor stopped: " << e.what();
                }
            });
            LOG(INFO) << "INCIDENT_CAPTURE started directory=" << directory
                      << " memory_mb=" << options.memory_mb << " pre_seconds=" << options.pre_seconds
                      << " post_seconds=" << options.post_seconds << " extra_dds_subscriptions=0";
            return true;
        } catch (const std::exception& e) {
            LOG(ERROR) << "INCIDENT_CAPTURE disabled: " << e.what();
            ++errors; Stop(); return false;
        }
    }
    void Timing(Sample& s, std::int64_t now) {
        s.steady = now; s.wall = Wall();
        s.prepared_age = Age(now, prepared); s.publish_begin_age = Age(now, publish_begin);
        s.publish_end_age = Age(now, publish_end); s.pipeline_age = Age(now, pipeline_at);
        s.publishing = publishing; s.dropped = busy + full + oversized;
        for (const auto& stream : streams) {
            if (stream->kind != Kind::kLidar && stream->kind != Kind::kImu) continue;
            const auto at = stream->received_at.load();
            const double age = Age(now, at > 0 ? at : started);
            if (stream->kind == Kind::kImu) s.imu_age = age;
            else s.lidar_age = std::max(s.lidar_age, age);
            const auto advanced = stream->advanced_at.load();
            if (at > 0) s.stamp_age = std::max(s.stamp_age, Age(now, advanced > 0 ? advanced : started));
        }
    }
    Sample Snapshot(std::int64_t now) {
        Sample s;
        { std::lock_guard<std::mutex> lock(snapshot_mutex); s.safety = latest_safety; s.pipeline = latest_pipeline; }
        Timing(s, now);
        return s;
    }
    void Trigger(unsigned reason, std::int64_t now) {
        std::lock_guard<std::mutex> lock(mutex);
        if (active) {
            if (!active->closed) {
                active->reasons.fetch_or(reason);
                active->until = std::min(active->trigger + std::int64_t(options.max_event_seconds * 1e9),
                                         now + std::int64_t(options.post_seconds * 1e9));
                ++merged;
            } else ++suppressed;
            return;
        }
        if (last_finished && Age(now, last_finished) < options.cooldown_seconds) { ++suppressed; return; }
        auto event = std::make_shared<Event>();
        event->trigger = now; event->wall = Wall(); event->reasons = reason;
        event->until = now + std::int64_t(options.post_seconds * 1e9);
        event->pending = history;  // references only, not a point-cloud copy
        if (!history.empty()) event->pre_coverage = Age(now, history.front()->steady);
        event->counts.resize(streams.size());
        event->first_written.resize(streams.size()); event->last_written.resize(streams.size());
        for (const auto& stream : streams) {
            event->trigger_callback_age.push_back(Age(now, stream->received_at));
            event->trigger_stamp_age.push_back(Age(now, stream->advanced_at));
        }
        active = std::move(event); ++events; wake.notify_one();
    }
    void MonitorLoop() {
        unsigned prior_issues = 0;
        while (!stop_monitor) {
            const auto now = Steady();
            {
                std::lock_guard<std::mutex> lock(mutex);
                while (!history.empty() &&
                       (Age(now, history.front()->steady) > options.pre_seconds ||
                        budget->held > options.memory_mb * MiB * 9 / 10)) history.pop_front();
                if (active && now >= active->until) { active->closed = true; wake.notify_one(); }
            }
            auto sample = std::make_shared<Sample>(Snapshot(now));
            Add(-1, sample, sizeof(Sample), sample->steady, sample->wall);
            unsigned issues = 0;
            if (sample->publish_end_age < 0) {
                if (Age(now, started) > options.startup_timeout_seconds) issues |= Startup;
            } else if (sample->publish_end_age > options.output_timeout_seconds) issues |= Output;
            if (Age(now, started) > options.startup_timeout_seconds || publish_end > 0) {
                if (std::max(sample->imu_age, sample->lidar_age) > options.input_timeout_seconds) issues |= Input;
                if (sample->stamp_age > options.input_timeout_seconds) issues |= FrozenStamp;
            }
            const auto reasons = (issues & ~prior_issues) | fault_edges.exchange(0);
            if (reasons) Trigger(reasons, now);
            prior_issues = issues;
            std::unique_lock<std::mutex> lock(mutex);
            wake.wait_for(lock, std::chrono::milliseconds(100), [this] { return stop_monitor.load(); });
        }
    }
    std::string Stage(const Sample& s) const {
        if (s.publishing && s.publish_begin_age > options.output_timeout_seconds) return "publish_call_pending";
        if (s.imu_age > options.input_timeout_seconds) return "primary_imu_callback_silence";
        if (s.lidar_age > options.input_timeout_seconds) return "lidar_callback_silence";
        if (s.stamp_age > options.input_timeout_seconds) return "input_stamp_not_advancing";
        if (s.pipeline.sensor_lag > options.input_timeout_seconds &&
            (s.pipeline.sensor_pending || s.pipeline.localization_pending)) return "processing_backlog";
        if (s.pipeline.match_stale && s.pipeline_age >= 0 && s.pipeline_age < 1) return "lidar_match_stale";
        if (s.safety.raw == 4 && s.safety.uncertainty == 2) return "uncertainty_check";
        if (s.safety.raw == 2) return "radar_diagnostic";
        if (s.publish_end_age > options.output_timeout_seconds) return "output_stale_unknown_upstream";
        return "unknown";
    }
    std::string StatsJson() const {
        const auto s = StatsNow();
        std::ostringstream out;
        out << "{\"held_bytes\":" << s.held_bytes << ",\"peak_bytes\":" << s.peak_bytes
            << ",\"received\":" << s.received << ",\"admitted\":" << s.admitted
            << ",\"dropped_busy\":" << s.dropped_busy << ",\"dropped_budget\":" << s.dropped_budget
            << ",\"dropped_oversize\":" << s.dropped_oversize << ",\"events\":" << s.events
            << ",\"merged\":" << s.merged << ",\"suppressed\":" << s.suppressed
            << ",\"writer_errors\":" << s.writer_errors << ",\"written_messages\":" << s.written << '}';
        return out.str();
    }
    void Report(const Event& event, const std::string& status) {
        std::string first = "unknown";
        double first_time = 0;
        for (const auto& s : event.samples) {
            const auto stage = Stage(s);
            if (stage != "unknown") { first = stage; first_time = (s.steady - event.trigger) * 1e-9; break; }
        }
        std::ostringstream json, rows, links;
        json << "{\"schema\":1,\"owner\":" << Json(Owner) << ",\"status\":" << Json(status)
             << ",\"trigger_wall_ns\":" << event.wall << ",\"reasons\":" << Reasons(event.reasons)
             << ",\"observed_first_stage\":" << Json(first) << ",\"first_stage_relative_sec\":" << Number(first_time)
             << ",\"requested_pre_seconds\":" << options.pre_seconds
             << ",\"retained_history_seconds\":" << Number(event.pre_coverage)
             << ",\"shutdown_shortened_window\":" << (event.shortened ? "true" : "false")
             << ",\"serialized_bytes\":" << event.bytes << ",\"omitted_messages\":" << event.omitted
             << ",\"closed_disk_bytes\":" << event.closed_disk_bytes
             << ",\"capture_loss_seen\":" << (busy + full + oversized > 0 ? "true" : "false")
             << ",\"error\":" << Json(event.error) << ",\"collector_totals\":" << StatsJson()
             << ",\"streams\":[";
        for (std::size_t i = 0; i < streams.size(); ++i) {
            if (i) json << ',';
            json << "{\"topic\":" << Json(streams[i]->topic) << ",\"type\":" << Json(streams[i]->type)
                 << ",\"written\":" << event.counts[i]
                 << ",\"callback_age_at_trigger\":" << Number(event.trigger_callback_age[i])
                 << ",\"stamp_age_at_trigger\":" << Number(event.trigger_stamp_age[i])
                 << ",\"first_relative_sec\":" << (event.first_written[i] ? Number((event.first_written[i] - event.trigger) * 1e-9) : "null")
                 << ",\"last_relative_sec\":" << (event.last_written[i] ? Number((event.last_written[i] - event.trigger) * 1e-9) : "null") << '}';
        }
        json << "],\"closed_segments\":[";
        for (std::size_t i = 0; i < event.segments.size(); ++i) {
            if (i) json << ',';
            json << Json(event.segments[i]);
            links << "<li><a href=\"" << event.segments[i] << "/metadata.yaml\">" << event.segments[i] << "</a></li>";
        }
        json << "]}\n";
        rows << '[';
        for (std::size_t i = 0; i < event.samples.size(); ++i) {
            const auto& s = event.samples[i]; if (i) rows << ',';
            rows << '[' << Number((s.steady - event.trigger) * 1e-9) << ',' << s.safety.raw << ',' << s.safety.published
                 << ',' << Number(s.imu_age) << ',' << Number(s.lidar_age) << ',' << Number(s.publish_end_age)
                 << ',' << s.pipeline.sensor_pending << ',' << s.pipeline.localization_pending
                 << ',' << Number(s.pipeline.sensor_lag) << ',' << Number(s.pipeline.match_age)
                 << ',' << s.dropped << ',' << Json(Stage(s)) << ']';
        }
        rows << ']';
        AtomicWrite(event.path / "summary.json", json.str());
        std::ostringstream html;
        html << R"HTML(<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>定位故障事件</title>
<style>body{font:16px system-ui;max-width:1200px;margin:32px auto;background:#111827;color:#e5e7eb;padding:20px}a{color:#93c5fd}canvas{width:100%;height:250px;background:#1f2937}table{border-collapse:collapse;width:100%;font-size:13px}td,th{padding:7px;border-bottom:1px solid #374151;text-align:left}tr:hover{background:#374151}pre{white-space:pre-wrap}input{padding:8px;width:300px}</style>
<h1>定位故障事件</h1><p id="headline"></p><p>最早观测异常用于缩小范围，不代表已证明硬件或网络根因。发布完成不代表下游已收到。</p>
<p><a href="summary.json">结构化摘要</a> · <a href="timeline.csv">完整时间线 CSV</a> · <a href="config.yaml">运行配置</a></p>
<details><summary>各话题实际保存范围与触发时输入状态</summary><pre id="coverage"></pre></details>
<canvas id="plot" width="1100" height="250"></canvas><p>时间为相对触发秒数。蓝：IMU 回调间隔；橙：最慢雷达回调间隔；绿：位姿发布完成间隔。负值表示尚未观测。</p>
<input id="filter" placeholder="筛选环节，例如 imu、backlog、publish"><table><thead><tr><th>秒</th><th>原始/上报</th><th>IMU/雷达间隔</th><th>输出间隔</th><th>输入/定位队列</th><th>输入落后</th><th>累计采集丢弃</th><th>观测环节</th></tr></thead><tbody id="rows"></tbody></table>
<h2>已关闭且可独立读取的 MCAP 分段</h2><ul>)HTML" << links.str()
             << "</ul><p>突然强杀可能留下未关闭分段及未导出的内存数据。未列出的分段需单独检查；summary 的状态可能停留在采集中。</p><script>const summary="
             << json.str() << ";const samples=" << rows.str() << R"JS(;
document.getElementById('headline').textContent=`状态 ${summary.status}；触发 ${summary.reasons.join(', ')}；最早观测 ${summary.observed_first_stage} (${summary.first_stage_relative_sec.toFixed(2)}s)；导出省略 ${summary.omitted_messages} 条；本次进程采集出现过丢弃：${summary.capture_loss_seen?'是':'否'}`;
const fmt=x=>typeof x==='number'?x.toFixed(3):'未知';
document.getElementById('coverage').textContent=summary.streams.map(s=>`${s.topic}: ${s.written} 条，范围 ${fmt(s.first_relative_sec)} 到 ${fmt(s.last_relative_sec)} s；触发时回调年龄 ${fmt(s.callback_age_at_trigger)} s，时间戳前进年龄 ${fmt(s.stamp_age_at_trigger)} s`).join('\n');
function render(){const q=document.getElementById('filter').value.toLowerCase();const root=document.getElementById('rows');root.replaceChildren();for(const r of samples){if(!r[11].includes(q))continue;const tr=document.createElement('tr');for(const v of [fmt(r[0]),`${r[1]}/${r[2]}`,`${fmt(r[3])}/${fmt(r[4])}`,fmt(r[5]),`${r[6]}/${r[7]}`,fmt(r[8]),r[10],r[11]]){const td=document.createElement('td');td.textContent=v;tr.append(td)}root.append(tr)}}document.getElementById('filter').oninput=render;render();
const c=document.getElementById('plot'),ctx=c.getContext('2d');if(samples.length>1){const a=samples[0][0],b=samples.at(-1)[0];const max=Math.max(1,...samples.flatMap(r=>r.slice(3,6).filter(x=>Number.isFinite(x)&&x>=0)));for(const [k,color] of [[3,'#60a5fa'],[4,'#fb923c'],[5,'#34d399']]){ctx.strokeStyle=color;ctx.beginPath();let started=false;for(const r of samples){if(!Number.isFinite(r[k])||r[k]<0){started=false;continue}const x=30+(r[0]-a)/Math.max(.01,b-a)*1040,y=220-r[k]/max*190;if(started)ctx.lineTo(x,y);else ctx.moveTo(x,y);started=true}ctx.stroke()}ctx.fillStyle='#e5e7eb';ctx.fillText(`${a.toFixed(1)}s → ${b.toFixed(1)}s; y最大 ${max.toFixed(2)}s`,30,240)}
</script></html>)JS";
        AtomicWrite(event.path / "summary.html", html.str());
    }
    void Prepare(Event& event) {
        std::vector<std::pair<fs::path, std::uint64_t>> old;
        std::uint64_t total = 0;
        for (const auto& entry : fs::directory_iterator(directory)) {
            if (entry.is_symlink() || !entry.is_directory() || entry.path().filename().string().rfind("event_", 0) != 0) continue;
            std::ifstream marker(entry.path() / "OWNER"); std::string owner; std::getline(marker, owner);
            if (owner != Owner) continue;
            auto bytes = DirectoryBytes(entry.path()); old.emplace_back(entry.path(), bytes); total += bytes;
        }
        std::sort(old.begin(), old.end());
        std::size_t count = old.size();
        for (const auto& item : old) {
            if (count < options.max_events && total + (options.event_mb + 1) * MiB <= options.disk_mb * MiB) break;
            // Only immediate, owned, non-symlink children of the locked archive.
            if (fs::weakly_canonical(item.first).parent_path() != directory) throw std::runtime_error("unsafe retention path");
            fs::remove_all(item.first); total -= item.second; --count;
        }
        if (fs::space(directory).available < (options.min_free_mb + std::max(options.segment_mb, options.max_message_mb) + 1) * MiB)
            throw std::runtime_error("insufficient free disk space for incident capture");
        event.path = directory / ("event_" + std::to_string(event.wall) + "_" + std::to_string(getpid()));
        if (!fs::create_directory(event.path)) throw std::runtime_error("incident directory already exists");
        AtomicWrite(event.path / "OWNER", std::string(Owner) + '\n');
        CopySmallFile(config, event.path / "config.yaml");
        if (!config.empty()) CopySmallFile(config.parent_path() / "run_metadata.txt", event.path / "run_metadata.txt");
        Report(event, "capturing");
    }
    void WriteEvent(const std::shared_ptr<Event>& event) {
        std::unique_ptr<rosbag2_cpp::Writer> bag;
        std::uint64_t segment_bytes = 0;
        std::int64_t opened_at = 0, pace_until = Steady();
        std::string segment;
        std::ofstream csv;
        auto close_segment = [&] {
            if (!bag) return;
            bag->close(); bag.reset();
            event->segments.push_back(segment); segment_bytes = 0;
            event->closed_disk_bytes = DirectoryBytes(event->path);
            csv.flush(); Report(*event, "capturing");
        };
        try {
            Prepare(*event);
            csv.open(event->path / "timeline.csv");
            csv << "relative_sec,wall_ns,raw,selected,published,strict,imu_callback_age,lidar_callback_age,stamp_advance_age,prepared_age,publish_begin_age,publish_end_age,publishing,pipeline_age,sensor_pending,localization_pending,high_frequency_pending,sensor_dropped,localization_dropped,sensor_lag,match_age,match_stale,output_gate,missing_lidars,lio_ms,major_sigma,lateral_sigma,yaw_sigma,covariance_valid,unobservable,collector_dropped,observed_stage,consistency,uncertainty,translation_condition,yaw_condition,confirmed_failure,unknown_tolerated,optional_can_unavailable,output_stale,localization_lost,radar_all_failed,non_sync_fault,non_sync_all_pass,bounded_policy_enabled,time_sync_fault_mask,bounded_sync_warning_mask,error_x,error_y,error_yaw_rad,output_age_sec,imu_arrival_age_sec,imu_stamp_lag_sec,can_arrival_age_sec,can_stamp_lag_sec,yaw_pass_age_sec\n";
            for (;;) {
                Row row; bool done;
                {
                    std::unique_lock<std::mutex> lock(mutex);
                    wake.wait_for(lock, std::chrono::milliseconds(100), [&] { return !event->pending.empty() || event->closed; });
                    if (!event->pending.empty()) { row = std::move(event->pending.front()); event->pending.pop_front(); }
                    done = event->closed && event->pending.empty() && !row;
                }
                if (done) break;
                if (bag && Age(Steady(), opened_at) >= 2) close_segment();
                if (!row) continue;
                if (row->stream < 0) {
                    const auto& s = *static_cast<const Sample*>(row->message.get());
                    event->samples.push_back(s);
                    csv << std::setprecision(12) << (s.steady - event->trigger) * 1e-9 << ',' << s.wall
                        << ',' << s.safety.raw << ',' << s.safety.selected << ',' << s.safety.published << ',' << s.safety.strict
                        << ',' << s.imu_age << ',' << s.lidar_age << ',' << s.stamp_age << ',' << s.prepared_age
                        << ',' << s.publish_begin_age << ',' << s.publish_end_age << ',' << s.publishing << ',' << s.pipeline_age
                        << ',' << s.pipeline.sensor_pending << ',' << s.pipeline.localization_pending << ',' << s.pipeline.high_frequency_pending
                        << ',' << s.pipeline.sensor_dropped << ',' << s.pipeline.localization_dropped << ',' << s.pipeline.sensor_lag
                        << ',' << s.pipeline.match_age << ',' << s.pipeline.match_stale << ',' << s.pipeline.output_gate << ',' << s.pipeline.missing_lidars
                        << ',' << s.pipeline.lio_ms
                        << ',' << s.safety.major_sigma << ',' << s.safety.lateral_sigma << ',' << s.safety.yaw_sigma
                        << ',' << s.safety.covariance_valid << ',' << s.safety.unobservable << ',' << s.dropped << ',' << Stage(s)
                        << ',' << s.safety.consistency << ',' << s.safety.uncertainty
                        << ',' << s.safety.translation_condition << ',' << s.safety.yaw_condition
                        << ',' << s.safety.confirmed_failure << ',' << s.safety.unknown_tolerated << ',' << s.safety.optional_can_unavailable
                        << ',' << s.safety.output_stale << ',' << s.safety.localization_lost << ',' << s.safety.radar_all_failed
                        << ',' << s.safety.non_sync_fault << ',' << s.safety.non_sync_all_pass << ',' << s.safety.bounded_policy_enabled
                        << ',' << s.safety.time_sync_fault_mask << ',' << s.safety.bounded_sync_warning_mask
                        << ',' << s.safety.error_x << ',' << s.safety.error_y << ',' << s.safety.error_yaw_rad
                        << ',' << s.safety.output_age_sec << ',' << s.safety.imu_arrival_age_sec << ',' << s.safety.imu_stamp_lag_sec
                        << ',' << s.safety.can_arrival_age_sec << ',' << s.safety.can_stamp_lag_sec << ',' << s.safety.yaw_pass_age_sec << '\n';
                    continue;
                }
                if (event->bytes >= options.event_mb * MiB || event->closed_disk_bytes >= options.event_mb * MiB) {
                    ++event->omitted; continue;
                }
                const auto& stream = *streams[row->stream];
                auto serialized = stream.serialize(row->message.get());
                const auto bytes = serialized->size();
                if (bytes > options.max_message_mb * MiB || event->bytes + bytes > options.event_mb * MiB) {
                    ++event->omitted; continue;
                }
                if (bag && segment_bytes + bytes > options.segment_mb * MiB) close_segment();
                if (event->closed_disk_bytes >= options.event_mb * MiB) { ++event->omitted; continue; }
                if (!bag) {
                    // Reserve a complete segment (or one oversized-to-segment message),
                    // not only the first small message in it. Filesystem use is checked
                    // at every close; disk_mb is a retention target with one-segment slack.
                    const auto reserve = std::max(options.segment_mb * MiB, bytes) + MiB;
                    if (fs::space(directory).available < options.min_free_mb * MiB + reserve)
                        throw std::runtime_error("free disk floor reached while exporting");
                    segment = "segment_" + std::to_string(event->segments.size());
                    rosbag2_storage::StorageOptions storage;
                    storage.uri = (event->path / segment).string(); storage.storage_id = "mcap";
                    storage.storage_preset_profile = "fastwrite"; storage.max_cache_size = 0;
                    bag = std::make_unique<rosbag2_cpp::Writer>(); bag->open(storage);
                    opened_at = Steady();
                }
                bag->write(serialized, stream.topic, stream.type, rclcpp::Time(row->wall));
                event->bytes += bytes; segment_bytes += bytes; ++event->counts[row->stream]; ++written;
                if (!event->first_written[row->stream]) event->first_written[row->stream] = row->steady;
                event->last_written[row->stream] = row->steady;
                // Average export-rate limit; only the writer sleeps. No compression.
                pace_until = std::max(pace_until, Steady()) + std::int64_t(bytes / (options.write_mb_per_second * MiB) * 1e9);
                std::this_thread::sleep_until(std::chrono::steady_clock::time_point(std::chrono::nanoseconds(pace_until)));
            }
            close_segment(); csv.flush();
            if (!csv) throw std::runtime_error("timeline write failed");
            Report(*event, event->omitted || event->shortened || busy + full + oversized > 0 ? "partial" : "complete");
            LOG(INFO) << "INCIDENT_CAPTURE saved=" << event->path << " bytes=" << event->bytes
                      << " omitted=" << event->omitted;
        } catch (const std::exception& e) {
            ++errors; event->error = e.what();
            LOG(ERROR) << "INCIDENT_CAPTURE export failed: " << e.what();
            try { close_segment(); if (!event->path.empty()) Report(*event, "error"); } catch (...) {}
        }
    }
    void WriterLoop() {
        setpriority(PRIO_PROCESS, static_cast<id_t>(syscall(SYS_gettid)), 10);
        std::int64_t last_status = 0;
        while (true) {
            std::shared_ptr<Event> event;
            {
                std::unique_lock<std::mutex> lock(mutex);
                wake.wait_for(lock, std::chrono::seconds(1), [this] { return active || stopping; });
                event = active;
                if (!event && stopping) break;
            }
            if (event) {
                WriteEvent(event);
                std::lock_guard<std::mutex> lock(mutex);
                active.reset(); last_finished = Steady();
            }
            if (Age(Steady(), last_status) > 5 || !last_status) {
                try { AtomicWrite(directory / "collector_status.json", StatsJson()); }
                catch (const std::exception& e) { ++errors; LOG_EVERY_N(ERROR, 10) << "INCIDENT_CAPTURE status: " << e.what(); }
                last_status = Steady();
            }
        }
    }
    void Stop() {
        stop_monitor = true; wake.notify_all();
        if (monitor.joinable()) monitor.join();
        { std::lock_guard<std::mutex> lock(mutex);
          if (active && !active->closed) { active->closed = true; active->shortened = true; }
          stopping = true; }
        wake.notify_all();
        if (writer.joinable()) writer.join();
        history.clear(); active.reset(); running = false;
        if (lock_fd >= 0) { flock(lock_fd, LOCK_UN); close(lock_fd); lock_fd = -1; }
    }
};

IncidentRecorder::IncidentRecorder(Options options, std::string directory, std::string config_path)
    : impl_(std::make_unique<Impl>(options, std::move(directory), std::move(config_path))) {}
IncidentRecorder::~IncidentRecorder() = default;
int IncidentRecorder::Register(const std::string& topic, const std::string& type, Kind kind, Serializer serializer) {
    if (impl_->running) throw std::logic_error("register incident streams before Start");
    auto stream = std::make_shared<Impl::Stream>();
    stream->topic = topic; stream->type = type; stream->kind = kind; stream->serialize = std::move(serializer);
    impl_->streams.push_back(stream); return static_cast<int>(impl_->streams.size() - 1);
}
bool IncidentRecorder::Start() { return impl_->Start(); }
bool IncidentRecorder::Capture(int stream, std::shared_ptr<const void> message, std::uint64_t bytes,
                               std::int64_t stamp) {
    if (!impl_->running || impl_->stop_monitor || !message || stream < 0 || std::size_t(stream) >= impl_->streams.size()) return false;
    const auto now = Steady(); const auto& s = impl_->streams[stream];
    s->received_at = now; ++s->count;
    const auto previous = s->sensor_stamp.exchange(stamp);
    if (stamp > 0 && stamp > previous) s->advanced_at = now;
    return impl_->Add(stream, std::move(message), bytes, now, Wall());
}
void IncidentRecorder::ObserveSafety(const Safety& snapshot) {
    if (!impl_->running) return;
    const bool changed = snapshot.raw != impl_->last_raw.exchange(snapshot.raw);
    if (changed) impl_->fault_edges.fetch_or(FaultBit(snapshot.raw));
    Impl::Sample sample;
    {
        std::unique_lock<std::mutex> lock(impl_->snapshot_mutex, std::try_to_lock);
        if (!lock.owns_lock()) { ++impl_->busy; return; }
        impl_->latest_safety = snapshot;
        if (changed) { sample.safety = snapshot; sample.pipeline = impl_->latest_pipeline; }
    }
    // Preserve short transitions between monitor ticks, including masked faults.
    if (changed) {
        impl_->Timing(sample, Steady());
        impl_->Add(-1, std::make_shared<Impl::Sample>(sample), sizeof(sample), sample.steady, sample.wall);
    }
}
void IncidentRecorder::ObservePipeline(const Pipeline& snapshot) {
    if (!impl_->running) return;
    std::unique_lock<std::mutex> lock(impl_->snapshot_mutex, std::try_to_lock);
    if (lock.owns_lock()) { impl_->latest_pipeline = snapshot; impl_->pipeline_at = Steady(); } else ++impl_->busy;
}
void IncidentRecorder::PosePrepared() { impl_->prepared = Steady(); }
void IncidentRecorder::PublishBegin() { impl_->publish_begin = Steady(); impl_->publishing = true; }
void IncidentRecorder::PublishEnd() { impl_->publish_end = Steady(); impl_->publishing = false; }
void IncidentRecorder::BeginShutdown() { impl_->stop_monitor = true; impl_->wake.notify_all(); }
void IncidentRecorder::Stop() { impl_->Stop(); }
IncidentRecorder::Stats IncidentRecorder::GetStats() const { return impl_->StatsNow(); }
}  // namespace lightning::diagnostics
