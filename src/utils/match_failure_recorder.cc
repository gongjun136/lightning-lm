#include "utils/match_failure_recorder.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstring>
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
#include <pthread.h>
#include <sched.h>
#include <sys/file.h>
#include <sys/resource.h>
#include <unistd.h>
#include <openssl/evp.h>
#include <glog/logging.h>

namespace lightning::diagnostics {
namespace {
namespace fs = std::filesystem;
constexpr std::uint64_t MiB = 1024 * 1024;
constexpr unsigned MaxValidation = 2000, MaxSource = 1000000;
constexpr char Owner[] = "LIGHTNING_LM_MATCH_FAILURE_V1";
using Recorder = MatchFailureRecorder;
std::int64_t Now(bool wall = false) {
    return wall ? std::chrono::duration_cast<std::chrono::nanoseconds>(
                      std::chrono::system_clock::now().time_since_epoch()).count()
                : std::chrono::duration_cast<std::chrono::nanoseconds>(
                      std::chrono::steady_clock::now().time_since_epoch()).count();
}
bool Finite(const Recorder::Point& p) { return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z); }
std::array<double, 3> Transform(const Recorder::Point& p, const Recorder::Pose& m) {
    return {m[0]*p.x+m[1]*p.y+m[2]*p.z+m[3], m[4]*p.x+m[5]*p.y+m[6]*p.z+m[7],
            m[8]*p.x+m[9]*p.y+m[10]*p.z+m[11]};
}
std::string Json(const std::string& text) {
    std::ostringstream out; out << '"';
    for (unsigned char c : text) {
        if (c == '"' || c == '\\') out << '\\' << c;
        else if (c < 32) out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << int(c) << std::dec;
        else out << c;
    }
    return out.str() + '"';
}
std::string Number(double n) {
    if (!std::isfinite(n)) return "null";
    std::ostringstream out; out << std::setprecision(17) << n; return out.str();
}
template <std::size_t N> std::string Array(const std::array<double, N>& values) {
    std::string out = "[";
    for (std::size_t i = 0; i < N; ++i) { if (i) out += ','; out += Number(values[i]); }
    return out + ']';
}
const char* ReasonName(Recorder::Reason reason) {
    switch (reason) {
        case Recorder::Reason::kOverlap: return "map_overlap";
        case Recorder::Reason::kBounds: return "map_bounds";
        case Recorder::Reason::kGravity: return "map_gravity";
        case Recorder::Reason::kOdometry: return "raw_odometry_disagreement";
        default: return "registration_rejected";
    }
}
void Sync(const fs::path& path, bool directory = false) {
    int fd = open(path.c_str(), O_RDONLY | O_CLOEXEC | (directory ? O_DIRECTORY : 0));
    if (fd < 0 || fsync(fd) != 0) { if (fd >= 0) close(fd); throw std::runtime_error("fsync failed: " + path.string()); }
    close(fd);
}
void WriteText(const fs::path& path, const std::string& data) {
    { std::ofstream out(path); out << data; out.flush(); if (!out) throw std::runtime_error("write failed: " + path.string()); }
    Sync(path);
}
std::uint64_t DirectoryBytes(const fs::path& path) {
    std::uint64_t bytes = 0;
    for (const auto& item : fs::recursive_directory_iterator(path))
        if (!item.is_symlink() && item.is_regular_file()) bytes += item.file_size();
    return bytes;
}
std::string ReadFirstLine(const fs::path& path) { std::ifstream f(path); std::string s; std::getline(f, s); return s; }
class Digest {
 public:
    Digest() : ctx(EVP_MD_CTX_new(), EVP_MD_CTX_free) {
        if (!ctx || EVP_DigestInit_ex(ctx.get(), EVP_sha256(), nullptr) != 1) throw std::runtime_error("SHA256 init failed");
    }
    void Add(const void* data, std::size_t size) {
        if (EVP_DigestUpdate(ctx.get(), data, size) != 1) throw std::runtime_error("SHA256 update failed");
    }
    std::string Finish() {
        unsigned char bytes[EVP_MAX_MD_SIZE]; unsigned length = 0;
        if (EVP_DigestFinal_ex(ctx.get(), bytes, &length) != 1) throw std::runtime_error("SHA256 final failed");
        std::ostringstream out; for (unsigned i = 0; i < length; ++i) out << std::hex << std::setw(2) << std::setfill('0') << unsigned(bytes[i]);
        return out.str();
    }
 private:
    std::unique_ptr<EVP_MD_CTX, decltype(&EVP_MD_CTX_free)> ctx;
};
void Put32(char* out, std::uint32_t value) { for (int k = 0; k < 4; ++k) out[k] = static_cast<char>((value >> (8*k)) & 255); }
void PutFloat(char* out, float value) { std::uint32_t bits; std::memcpy(&bits, &value, 4); Put32(out, bits); }
void EncodePoint(char* out, const Recorder::Point& p) { PutFloat(out, p.x); PutFloat(out+4, p.y); PutFloat(out+8, p.z); }
std::string ConfigureWriter() {
    sched_param param{};
    if (pthread_setschedparam(pthread_self(), SCHED_OTHER, &param) != 0) throw std::runtime_error("cannot set writer SCHED_OTHER");
    const int nice_value = std::max(10, getpriority(PRIO_PROCESS, 0));
    if (setpriority(PRIO_PROCESS, 0, nice_value) != 0) throw std::runtime_error("cannot lower writer priority");
    cpu_set_t allowed; CPU_ZERO(&allowed);
    if (sched_getaffinity(0, sizeof(allowed), &allowed) != 0) throw std::runtime_error("cannot read writer affinity");
    CPU_CLR(10, &allowed); CPU_CLR(11, &allowed);
    std::ifstream isolated("/sys/devices/system/cpu/isolated"); std::string list; std::getline(isolated, list);
    std::istringstream ranges(list); std::string range;
    while (std::getline(ranges, range, ',')) {
        if (range.empty()) continue;
        const auto dash = range.find('-'); int a = std::stoi(range), b = dash == std::string::npos ? a : std::stoi(range.substr(dash+1));
        for (int cpu = std::max(0, a); cpu <= b && cpu < CPU_SETSIZE; ++cpu) CPU_CLR(cpu, &allowed);
    }
    if (!CPU_COUNT(&allowed) || sched_setaffinity(0, sizeof(allowed), &allowed) != 0)
        throw std::runtime_error("no permitted non-isolated CPU for writer; affinity was not widened");
    std::ostringstream out; out << "{\"policy\":\"SCHED_OTHER\",\"nice\":" << nice_value << ",\"cpus\":[";
    bool first = true; for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) if (CPU_ISSET(cpu, &allowed)) { if (!first) out << ','; first = false; out << cpu; }
    return out.str() + "]}";
}
}  // namespace

MatchFailureRecorder::Options MatchFailureRecorder::Options::FromYaml(const YAML::Node& node) {
    Options options; if (!node) return options;
#define READ(name) if (node[#name]) options.name = node[#name].as<decltype(options.name)>()
    READ(enabled); READ(cooldown_seconds); READ(recovery_seconds); READ(protect_days); READ(write_mb_per_second);
    READ(disk_mb); READ(min_free_mb); READ(event_mb); READ(max_scan_points);
#undef READ
    options.Validate(); return options;
}
void MatchFailureRecorder::Options::Validate() const {
    if (!std::isfinite(cooldown_seconds) || cooldown_seconds < 1 || cooldown_seconds > 3600 ||
        !std::isfinite(recovery_seconds) || recovery_seconds < 1 || recovery_seconds > 60 ||
        !std::isfinite(write_mb_per_second) || write_mb_per_second < 1 || write_mb_per_second > 1000 ||
        protect_days < 7 || protect_days > 3650 || disk_mb < 1 || disk_mb > 1048576 ||
        min_free_mb > 1048576 || event_mb < 1 || event_mb > 64 || disk_mb < event_mb ||
        max_scan_points < 1 || max_scan_points > 50000)
        throw std::invalid_argument("match_failure_capture: invalid resource limits or protection shorter than 7 days");
}

class MatchFailureRecorder::Impl {
 public:
    struct SavedPoint { Point point; std::uint32_t index = 0, validation = 0; };
    struct Slot {
        enum State { kEmpty, kReady, kWriting } state = kEmpty;
        Frame frame; std::vector<SavedPoint> points;
        std::uint64_t source_points = 0, sequence = 0, validation_points = 0;
        std::int64_t wall = 0; bool downsampled = false;
    };
    Options options; Context context; CloudView map;
    std::array<Slot, 2> slots;
    std::mutex mutex; std::condition_variable wake; std::thread writer;
    std::atomic<bool> running{false}, stopping{false};
    std::atomic<std::int64_t> last_capture{0}, good_since{0};
    std::atomic<std::uint64_t> received{0}, admitted{0}, saved{0}, suppressed{0}, busy{0}, queue_full{0}, oversize{0}, disk_rejected{0}, map_rejected{0}, errors{0};
    const std::int64_t session = Now(true);
    int lock_fd = -1; std::string scheduling = "null", map_identity, last_error;
    explicit Impl(Options o, Context c, CloudView m) : options(o), context(std::move(c)), map(m) {
        options.Validate();
        if (context.config_snapshot.size() > MiB) throw std::invalid_argument("match failure configuration snapshot exceeds 1 MiB");
        for (auto& slot : slots) slot.points.reserve(options.max_scan_points + MaxValidation);
    }
    ~Impl() { Stop(); if (lock_fd >= 0) close(lock_fd); }
    Stats GetStats() const { return {received.load(), admitted.load(), saved.load(), suppressed.load(), busy.load(), queue_full.load(), oversize.load(), disk_rejected.load(), map_rejected.load(), errors.load()}; }
    bool Start() {
        if (!options.enabled || running || stopping || writer.joinable()) return false;
        if (!map.data || !map.point || !map.size || context.directory.empty()) return false;
        try {
            fs::path directory(context.directory);
            if (!directory.is_absolute()) throw std::runtime_error("match failure directory must be absolute");
            fs::create_directories(directory); context.directory = fs::canonical(directory).string();
            lock_fd = open((fs::path(context.directory)/".capture.lock").c_str(), O_CREAT|O_RDWR|O_CLOEXEC, 0600);
            if (lock_fd < 0 || flock(lock_fd, LOCK_EX|LOCK_NB) != 0) throw std::runtime_error("match failure directory is already in use");
            running = true; writer = std::thread([this] { Loop(); }); return true;
        } catch (const std::exception& e) {
            running = false;
            if (lock_fd >= 0) { close(lock_fd); lock_fd = -1; }
            ++errors; LOG(ERROR) << "MATCH_FAILURE_CAPTURE start failed: " << e.what(); return false;
        }
    }
    void Observe(bool accepted) {
        if (!running) return;
        if (!accepted) { good_since = 0; return; }
        const auto now = Now(); std::int64_t zero = 0; good_since.compare_exchange_strong(zero, now);
        if (now - good_since.load() >= options.recovery_seconds * 1e9) last_capture = 0;
    }
    bool Capture(const CloudView& input, const Frame& frame) {
        if (!running || stopping) return false;
        ++received; const auto now = Now();
        if (last_capture && now - last_capture.load() < options.cooldown_seconds * 1e9) { ++suppressed; return false; }
        if (!input.data || !input.point || !input.size || input.size > MaxSource ||
            frame.validation_points > MaxValidation || !frame.validation_stride) { ++oversize; return false; }
        std::unique_lock<std::mutex> lock(mutex, std::try_to_lock);
        if (!lock.owns_lock()) { ++busy; return false; }
        if (last_capture && now - last_capture.load() < options.cooldown_seconds * 1e9) { ++suppressed; return false; }
        Slot* slot = nullptr; for (auto& item : slots) if (item.state == Slot::kEmpty) { slot = &item; break; }
        if (!slot) { ++queue_full; return false; }
        slot->points.clear(); slot->validation_points = 0;
        const auto display_stride = std::max<std::size_t>(1, (input.size + options.max_scan_points - 1) / options.max_scan_points);
        std::size_t finite_index = 0;
        for (std::size_t i = 0; i < input.size; ++i) {
            const Point point = input.point(input.data, i); if (!Finite(point)) continue;
            const bool validation = frame.evaluated && frame.validation_points > 0 && finite_index % frame.validation_stride == 0;
            const bool display = finite_index % display_stride == 0; ++finite_index;
            if (!display && !validation) continue;
            if (slot->points.size() >= options.max_scan_points + MaxValidation) { ++oversize; slot->points.clear(); return false; }
            if (validation) ++slot->validation_points;
            slot->points.push_back({point, static_cast<std::uint32_t>(i), validation ? 1u : 0u});
        }
        if (slot->validation_points != (frame.evaluated ? frame.validation_points : 0) || slot->points.empty()) {
            ++oversize; slot->points.clear(); return false;
        }
        slot->frame = frame; slot->source_points = input.size; slot->downsampled = display_stride > 1;
        slot->wall = Now(true); slot->sequence = ++admitted; slot->state = Slot::kReady;
        last_capture = now; wake.notify_one(); return true;
    }
    void Status() {
        const Stats s = GetStats();
        std::ostringstream out; out << "{\"schema\":1,\"kind\":\"lightning_match_failure_collector\",\"pid\":" << getpid()
            << ",\"session_start_wall_ns\":" << session << ",\"run_dir\":" << Json(context.run_dir)
            << ",\"running\":" << (running && !stopping ? "true" : "false")
            << ",\"received\":" << s.received << ",\"admitted\":" << s.admitted << ",\"saved\":" << s.saved
            << ",\"suppressed\":" << s.suppressed << ",\"dropped_busy\":" << s.dropped_busy
            << ",\"dropped_queue\":" << s.dropped_queue << ",\"dropped_oversize\":" << s.dropped_oversize
            << ",\"rejected_disk\":" << s.rejected_disk << ",\"rejected_map\":" << s.rejected_map
            << ",\"writer_errors\":" << s.writer_errors << ",\"last_error\":" << Json(last_error)
            << ",\"protection_days\":" << options.protect_days << ",\"disk_mb\":" << options.disk_mb
            << ",\"writer_scheduling\":" << scheduling << "}\n";
        const fs::path path = fs::path(context.directory)/"collector_status.json";
        WriteText(path.string()+".tmp", out.str()); fs::rename(path.string()+".tmp", path);
    }
    bool Reserve(std::uint64_t needed) {
        struct Event { fs::path path; std::uint64_t bytes; std::int64_t wall; };
        std::vector<Event> events;
        const fs::path root(context.directory);
        // Crashed .partial packages, foreign files and rendering caches consume
        // the same quota. Only explicitly owned, expired events may be removed.
        std::uint64_t used = DirectoryBytes(root);
        for (const auto& item : fs::directory_iterator(root)) {
            if (item.is_symlink() || !item.is_directory() || item.path().filename().string().rfind("match_failure_", 0) != 0 ||
                fs::is_symlink(item.path()/"OWNER") ||
                ReadFirstLine(item.path()/"OWNER") != Owner) continue;
            const auto bytes = DirectoryBytes(item.path());
            std::int64_t wall = 0;
            try { wall = std::stoll(ReadFirstLine(item.path()/"captured_wall_ns")); } catch (...) {}
            events.push_back({item.path(), bytes, wall});
        }
        std::sort(events.begin(), events.end(), [](const Event& a, const Event& b) { return a.wall < b.wall; });
        const auto cutoff = Now(true) - static_cast<std::int64_t>(options.protect_days)*86400*1000000000LL;
        for (const auto& event : events) {
            if (used + needed <= options.disk_mb*MiB && fs::space(root).available >= options.min_free_mb*MiB+needed) break;
            // Missing/invalid timestamps and all protected events are retained.
            if (event.wall <= 0 || event.wall > cutoff) continue;
            if (fs::weakly_canonical(event.path).parent_path() != root) throw std::runtime_error("unsafe match evidence retention path");
            fs::remove_all(event.path); used -= event.bytes;
        }
        return used + needed <= options.disk_mb*MiB && fs::space(root).available >= options.min_free_mb*MiB+needed;
    }
    void Save(const Slot& slot) {
        const Frame& f = slot.frame;
        std::array<double,3> low, high; low.fill(std::numeric_limits<double>::infinity()); high.fill(-std::numeric_limits<double>::infinity());
        for (const auto& p : slot.points) for (const auto* pose : {&f.guess, &f.ndt, &f.check}) {
            const auto transformed = Transform(p.point, *pose);
            for (unsigned axis = 0; axis < 3; ++axis) {
                if (!std::isfinite(transformed[axis])) throw std::runtime_error("nonfinite captured pose");
                low[axis] = std::min(low[axis], transformed[axis]); high[axis] = std::max(high[axis], transformed[axis]);
            }
        }
        // Include float conversion error in PCL's transformed XYZ as well as the
        // accepted radius. This crop preserves every possible nearest hit.
        for (unsigned k = 0; k < 3; ++k) {
            const double margin = f.nearest_neighbor_distance_m + 1e-4 + 2*std::numeric_limits<float>::epsilon()*std::max(std::abs(low[k]), std::abs(high[k]));
            low[k] -= margin; high[k] += margin;
        }
        auto inside = [&](const Point& p) { return Finite(p) && p.x >= low[0] && p.x <= high[0] && p.y >= low[1] && p.y <= high[1] && p.z >= low[2] && p.z <= high[2]; };
        std::uint64_t map_points = 0;
        for (std::size_t i = 0; i < map.size; ++i) if (inside(map.point(map.data, i))) ++map_points;
        const std::uint64_t needed = slot.points.size()*24 + map_points*12 + context.config_snapshot.size() + 65536;
        if (needed > options.event_mb*MiB) { ++map_rejected; last_error = "local map exceeds complete event budget; evidence was not truncated"; LOG(WARNING) << "MATCH_FAILURE_CAPTURE " << last_error; return; }
        if (!Reserve(needed)) { ++disk_rejected; last_error = "protected evidence/disk budget or minimum free space prevents new capture"; LOG(WARNING) << "MATCH_FAILURE_CAPTURE " << last_error; return; }
        const fs::path root(context.directory);
        const std::string name = "match_failure_" + std::to_string(slot.wall) + "_" + std::to_string(getpid()) + "_" + std::to_string(slot.sequence);
        const fs::path temporary = root/("."+name+".partial"), final = root/name;
        if (!fs::create_directory(temporary)) throw std::runtime_error("duplicate evidence staging path");
        try {
            WriteText(temporary/"OWNER", std::string(Owner)+"\n");
            WriteText(temporary/"captured_wall_ns", std::to_string(slot.wall)+"\n");
            if (!context.config_snapshot.empty()) WriteText(temporary/"config.yaml",context.config_snapshot);
            Digest scan_hash, map_hash;
            const auto write_started = std::chrono::steady_clock::now();
            std::uint64_t written_bytes = 0, next_rate_check = 65536;
            auto throttle = [&](std::uint64_t bytes) {
                written_bytes += bytes;
                if (written_bytes < next_rate_check) return;
                next_rate_check = written_bytes + 65536;
                const auto target = write_started + std::chrono::duration_cast<std::chrono::steady_clock::duration>(
                    std::chrono::duration<double>(written_bytes/(options.write_mb_per_second*MiB)));
                std::this_thread::sleep_until(target);
            };
            { std::ofstream scan(temporary/"scan.bin", std::ios::binary);
              for (const auto& p : slot.points) { char bytes[24]; EncodePoint(bytes, p.point); Put32(bytes+12, p.point.lidar_id); Put32(bytes+16, p.index); Put32(bytes+20, p.validation); scan.write(bytes, sizeof(bytes)); scan_hash.Add(bytes,sizeof(bytes)); throttle(sizeof(bytes)); }
              scan.flush(); if (!scan) throw std::runtime_error("scan write failed"); }
            { std::ofstream target(temporary/"map_local.bin", std::ios::binary);
              for (std::size_t i = 0; i < map.size; ++i) { const auto p = map.point(map.data, i); if (!inside(p)) continue; char bytes[12]; EncodePoint(bytes,p); target.write(bytes,sizeof(bytes)); map_hash.Add(bytes,sizeof(bytes)); throttle(sizeof(bytes)); }
              target.flush(); if (!target) throw std::runtime_error("map write failed"); }
            Sync(temporary/"scan.bin"); Sync(temporary/"map_local.bin");
            std::ostringstream s;
            s << "{\"schema\":1,\"kind\":\"lightning_match_failure\",\"status\":\"complete\",\"wall_ns\":" << slot.wall
              << ",\"sensor_stamp_ns\":" << f.sensor_stamp_ns << ",\"pid\":" << getpid() << ",\"sequence\":" << slot.sequence
              << ",\"session_start_wall_ns\":" << session << ",\"config_path\":" << Json(context.config_path) << ",\"run_dir\":" << Json(context.run_dir)
              << ",\"synthetic_fixture\":" << (context.synthetic_fixture ? "true" : "false")
              << ",\"reason\":" << Json(ReasonName(f.reason)) << ",\"build\":{\"commit\":" << Json(context.build_commit) << ",\"dirty\":" << Json(context.build_dirty) << ",\"built_utc\":" << Json(context.build_utc) << "}"
              << ",\"poses\":{\"guess\":" << Array(f.guess) << ",\"ndt\":" << Array(f.ndt) << ",\"check\":" << Array(f.check) << "}"
              << ",\"pose_frame\":\"internal_map_from_primary_lidar_at_scan_end\",\"map\":{\"path\":" << Json(context.map_path)
              << ",\"source_points\":" << map.size << ",\"bounds_min\":" << Array(context.map_min) << ",\"bounds_max\":" << Array(context.map_max)
              << ",\"local_bounds_min\":" << Array(low) << ",\"local_bounds_max\":" << Array(high) << ",\"local_points\":" << map_points
              << ",\"coverage\":\"all_possible_neighbors_for_saved_scan_and_all_three_poses\",\"identity_sha256\":" << Json(map_identity)
              << ",\"identity_encoding\":\"loaded_map_ordered_xyz_f32_le\",\"record_bytes\":12}"
              << ",\"scan\":{\"source_points\":" << slot.source_points << ",\"saved_points\":" << slot.points.size() << ",\"validation_points\":" << slot.validation_points
              << ",\"display_downsampled\":" << (slot.downsampled ? "true" : "false") << ",\"record_bytes\":24,\"encoding\":\"little_endian_xyz_f32_lidar_u32_source_index_u32_validation_u32\"}"
              << ",\"criterion\":{\"evaluated\":" << (f.evaluated ? "true" : "false")
              << ",\"sample_replay_available\":" << (f.evaluated && f.validation_points > 0 ? "true" : "false")
              << ",\"sampling_stride\":" << f.validation_stride << ",\"points\":" << f.validation_points
              << ",\"overlap\":" << Number(f.overlap) << ",\"inside_xy\":" << Number(f.inside_xy) << ",\"inside_xyz\":" << Number(f.inside_xyz) << ",\"gravity_cos\":" << Number(f.gravity_cos)
              << ",\"min_overlap\":" << Number(f.min_overlap) << ",\"min_inside_xy\":" << Number(f.min_inside_xy) << ",\"min_gravity_cos\":" << Number(f.min_gravity_cos)
              << ",\"nearest_neighbor_distance_m\":" << Number(f.nearest_neighbor_distance_m) << ",\"bounds_margin_m\":" << Number(f.bounds_margin_m) << "}"
              << ",\"registration\":{\"confidence\":" << Number(f.confidence) << ",\"iterations\":" << f.iterations << ",\"raw_odom_valid\":" << (f.raw_odom_valid ? "true" : "false") << ",\"raw_odom_delta\":" << Number(f.raw_odom_delta) << "}"
              << ",\"files\":{\"scan\":\"scan.bin\",\"map\":\"map_local.bin\",\"config\":" << (context.config_snapshot.empty() ? "null" : "\"config.yaml\"") << "},\"files_sha256\":{\"scan\":" << Json(scan_hash.Finish()) << ",\"map\":" << Json(map_hash.Finish()) << "}"
              << ",\"protection_days\":" << options.protect_days << ",\"writer_scheduling\":" << scheduling << "}\n";
            WriteText(temporary/"summary.json", s.str()); Sync(temporary,true);
            fs::rename(temporary, final); Sync(root,true); ++saved; last_error.clear();
            LOG(INFO) << "MATCH_FAILURE_CAPTURE saved=" << final << " sensor_stamp_ns=" << f.sensor_stamp_ns << " overlap=" << f.overlap << " validation_points=" << slot.validation_points;
        } catch (...) {
            // An incomplete package is never presented as complete evidence.
            std::error_code ec; fs::remove_all(temporary,ec); throw;
        }
    }
    void Loop() {
        try {
            scheduling = ConfigureWriter();
            Digest identity;
            for (std::size_t i = 0; i < map.size; ++i) { char bytes[12]; EncodePoint(bytes,map.point(map.data,i)); identity.Add(bytes,sizeof(bytes)); }
            map_identity = identity.Finish(); Status();
            while (true) {
                Slot* slot = nullptr;
                {
                    std::unique_lock<std::mutex> lock(mutex);
                    wake.wait_for(lock,std::chrono::seconds(1),[&] { return stopping || slots[0].state == Slot::kReady || slots[1].state == Slot::kReady; });
                    for (auto& item : slots) if (item.state == Slot::kReady && (!slot || item.sequence < slot->sequence)) slot = &item;
                    if (slot) slot->state = Slot::kWriting;
                    else if (stopping) break;
                }
                if (slot) {
                    try { Save(*slot); } catch (const std::exception& e) { ++errors; last_error = e.what(); LOG(ERROR) << "MATCH_FAILURE_CAPTURE writer: " << e.what(); }
                    std::lock_guard<std::mutex> lock(mutex); slot->points.clear(); slot->state = Slot::kEmpty;
                }
                Status();
            }
        } catch (const std::exception& e) { ++errors; last_error = e.what(); LOG(ERROR) << "MATCH_FAILURE_CAPTURE writer stopped: " << e.what(); }
        running = false;
        try { Status(); } catch (...) { ++errors; }
    }
    void Stop() {
        stopping = true; wake.notify_all(); if (writer.joinable()) writer.join(); running = false;
        if (lock_fd >= 0) { close(lock_fd); lock_fd = -1; }
    }
};

MatchFailureRecorder::MatchFailureRecorder(Options options, Context context, CloudView map)
    : impl_(std::make_unique<Impl>(options,std::move(context),map)) {}
MatchFailureRecorder::~MatchFailureRecorder() = default;
bool MatchFailureRecorder::Start() { return impl_->Start(); }
void MatchFailureRecorder::ObserveTracking(bool accepted) { impl_->Observe(accepted); }
bool MatchFailureRecorder::Capture(const CloudView& input,const Frame& frame) { return impl_->Capture(input,frame); }
void MatchFailureRecorder::Stop() { impl_->Stop(); }
MatchFailureRecorder::Stats MatchFailureRecorder::GetStats() const { return impl_->GetStats(); }
}  // namespace lightning::diagnostics
