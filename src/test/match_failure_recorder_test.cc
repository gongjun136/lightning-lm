#include "utils/match_failure_recorder.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <thread>
#include <vector>
#include <unistd.h>
#include <glog/logging.h>

namespace {
namespace fs = std::filesystem;
using R = lightning::diagnostics::MatchFailureRecorder;
using Cloud = std::vector<R::Point>;
constexpr std::uint64_t MiB = 1024*1024;
void Check(bool condition, const std::string& message) { if (!condition) throw std::runtime_error(message); }
R::Point PointAt(const void* data, std::size_t i) { return static_cast<const Cloud*>(data)->at(i); }
R::CloudView View(const Cloud& cloud) { return {&cloud,cloud.size(),PointAt}; }
R::Pose Pose(double x = 0) { return {1,0,0,x,0,1,0,0,0,0,1,0,0,0,0,1}; }
std::int64_t Wall() { return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count(); }
struct Fixture {
    fs::path root;
    Fixture() { char path[] = "/tmp/lightning_match_failure.XXXXXX"; const auto p = mkdtemp(path); Check(p,"mkdtemp"); root = p; }
    ~Fixture() { std::error_code ec; fs::remove_all(root,ec); }
};
R::Options Options() {
    R::Options o; o.enabled = true; o.cooldown_seconds = 1; o.recovery_seconds = 1;
    o.min_free_mb = 0; o.disk_mb = 16; o.event_mb = 8; o.write_mb_per_second = 1000;
    return o;
}
R::Context Context(const fs::path& root) {
    R::Context c; c.directory = (root/"evidence").string(); c.config_path = (root/"config.yaml").string();
    c.run_dir = root.string(); c.map_path = (root/"map/global.pcd").string(); c.map_min = {-10,-10,-10}; c.map_max = {100,100,100};
    c.build_commit = "native-test"; c.build_dirty = "1"; c.build_utc = "test"; c.synthetic_fixture = true; return c;
}
R::Frame Frame(const Cloud& cloud) {
    R::Frame f; f.sensor_stamp_ns = 1793348973600000000LL; f.guess = Pose(1); f.ndt = Pose(2); f.check = Pose();
    f.reason = R::Reason::kOverlap; f.evaluated = true; f.validation_stride = std::max<std::size_t>(1,(cloud.size()+1999)/2000);
    std::size_t finite = 0;
    for (auto p : cloud) if (std::isfinite(p.x)&&std::isfinite(p.y)&&std::isfinite(p.z)) {
        if (finite % f.validation_stride == 0) ++f.validation_points;
        ++finite;
    }
    f.overlap = 1; f.inside_xy = f.inside_xyz = f.gravity_cos = 1; f.raw_odom_valid = true;
    return f;
}
std::vector<fs::path> Events(const fs::path& root) {
    std::vector<fs::path> found;
    if (!fs::is_directory(root)) return found;
    for (const auto& item : fs::directory_iterator(root))
        if (item.is_directory() && item.path().filename().string().rfind("match_failure_",0) == 0 && fs::exists(item.path()/"summary.json")) found.push_back(item.path());
    std::sort(found.begin(),found.end()); return found;
}
std::vector<char> Read(const fs::path& path) { std::ifstream file(path,std::ios::binary); return {std::istreambuf_iterator<char>(file),{}}; }
std::uint32_t U32(const char* bytes) { std::uint32_t value = 0; for (int i=0;i<4;++i) value |= std::uint32_t(static_cast<unsigned char>(bytes[i])) << (i*8); return value; }
float F32(const char* bytes) { const auto bits = U32(bytes); float value; std::memcpy(&value,&bits,4); return value; }
void FakeOwned(const fs::path& root, const std::string& name, std::int64_t wall, std::uint64_t size) {
    fs::create_directories(root/name);
    std::ofstream(root/name/"OWNER") << "LIGHTNING_LM_MATCH_FAILURE_V1\n";
    std::ofstream(root/name/"captured_wall_ns") << wall << '\n';
    { std::ofstream file(root/name/"data"); file.seekp(size-1); file.put('x'); }
}

fs::path ExactSelection(const fs::path& root) {
    Cloud cloud;
    for (unsigned i=0;i<6503;++i) cloud.push_back({float(i%100)*.1f,float(i/100)*.1f,float(i%7)*.1f,i%3});
    cloud[3].x = cloud[7].x = std::numeric_limits<float>::quiet_NaN();
    Cloud map;
    for (auto p : cloud) if (std::isfinite(p.x)) { map.push_back(p); p.x += 2; map.push_back(p); }
    map.push_back({1000,1000,1000,0});
    auto o = Options(); o.max_scan_points = 4;
    const auto context = Context(root); R recorder(o,context,View(map));
    Check(recorder.Start(),"exact capture start"); const auto frame = Frame(cloud);
    Check(recorder.Capture(View(cloud),frame),"first failure was not admitted");
    Check(!recorder.Capture(View(cloud),frame),"continuous failure bypassed cooldown"); recorder.Stop();
    Check(recorder.GetStats().saved == 1,"exact capture did not save");
    const auto events = Events(context.directory); Check(events.size()==1,"event count");
    const auto path = events.front(); const auto metadata = YAML::LoadFile((path/"summary.json").string());
    Check(metadata["status"].as<std::string>()=="complete","partial exposed");
    Check(metadata["run_dir"].as<std::string>()==root.string(),"run identity missing");
    Check(metadata["map"]["local_points"].as<std::size_t>()==map.size()-1,"three-pose local map coverage wrong");
    Check(metadata["map"]["identity_sha256"].as<std::string>().size()==64,"map identity missing");
    const auto bytes=Read(path/"scan.bin"); Check(bytes.size()%24==0,"scan layout");
    std::vector<unsigned> selected;
    for (std::size_t i=0;i<bytes.size();i+=24) {
        const auto index=U32(&bytes[i+16]); Check(index<cloud.size(),"invalid source index");
        Check(F32(&bytes[i])==cloud[index].x && F32(&bytes[i+4])==cloud[index].y && F32(&bytes[i+8])==cloud[index].z,"scan altered");
        if (U32(&bytes[i+20])) selected.push_back(index);
    }
    std::vector<unsigned> expected; unsigned finite=0;
    for (unsigned i=0;i<cloud.size();++i) if (std::isfinite(cloud[i].x)) { if (finite%frame.validation_stride==0) expected.push_back(i); ++finite; }
    Check(selected==expected,"exact validation sampling lost after display downsampling/NaNs");
    Check(metadata["scan"]["validation_points"].as<std::size_t>()==expected.size(),"validation metadata mismatch");
    Check(metadata["writer_scheduling"]["policy"].as<std::string>()=="SCHED_OTHER","writer inherited realtime policy");
    for (const auto& cpu : metadata["writer_scheduling"]["cpus"]) Check(cpu.as<int>()!=10 && cpu.as<int>()!=11,"writer occupies isolated critical CPU");
    return path;
}

void ProtectionAndRetention() {
    Fixture f; auto context=Context(f.root); auto options=Options(); options.disk_mb=options.event_mb=1;
    const fs::path archive(context.directory); FakeOwned(archive,"match_failure_protected",Wall(),MiB);
    Cloud cloud{{0,0,0,0},{.1f,0,0,0}};
    {
        R recorder(options,context,View(cloud)); Check(recorder.Start(),"protected start"); Check(recorder.Capture(View(cloud),Frame(cloud)),"protected admission"); recorder.Stop();
        Check(recorder.GetStats().rejected_disk==1,"protected evidence not budget rejected");
        Check(fs::exists(archive/"match_failure_protected/data"),"protected evidence deleted");
        Check(Events(archive).empty(),"full budget exposed incomplete event");
    }
    const auto old=Wall()-8LL*86400*1000000000LL;
    std::ofstream(archive/"match_failure_protected/captured_wall_ns")<<old;
    fs::create_directories(archive/"foreign"); std::ofstream(archive/"foreign/keep")<<"user evidence";
    const auto external=f.root/"external"; fs::create_directories(external); std::ofstream(external/"keep")<<"outside";
    fs::create_directory_symlink(external,archive/"match_failure_symlink");
    {
        R recorder(options,context,View(cloud)); Check(recorder.Start(),"retention restart"); Check(recorder.Capture(View(cloud),Frame(cloud)),"retention admission"); recorder.Stop();
        Check(recorder.GetStats().saved==1,"expired owned evidence not reclaimed");
        Check(!fs::exists(archive/"match_failure_protected"),"expired evidence retained under pressure");
        Check(fs::exists(archive/"foreign/keep") && fs::exists(external/"keep") && fs::is_symlink(archive/"match_failure_symlink"),"retention touched unowned paths");
    }
}

void OversizeMapAndQueue() {
    Fixture f; auto context=Context(f.root); auto options=Options(); options.event_mb=1;
    Cloud map(100000,{0,0,0,0}), cloud{{0,0,0,0}};
    R recorder(options,context,View(map)); Check(recorder.Start(),"map cap start");
    Check(recorder.Capture(View(cloud),Frame(cloud)),"map cap admission"); recorder.Stop();
    Check(recorder.GetStats().rejected_map==1 && Events(context.directory).empty(),"oversize map truncated or published");

    Fixture f2; options=Options(); options.event_mb=8; options.write_mb_per_second=1;
    context=Context(f2.root); Cloud large_map(400000,{0,0,0,0});
    R slow(options,context,View(large_map)); Check(slow.Start(),"queue test start");
    auto frame=Frame(cloud); Check(slow.Capture(View(cloud),frame),"queue first capture");
    std::this_thread::sleep_for(std::chrono::milliseconds(1100));
    frame.reason=R::Reason::kBounds;
    Check(slow.Capture(View(cloud),frame),"queue second capture");
    std::this_thread::sleep_for(std::chrono::milliseconds(1100));
    frame.reason=R::Reason::kGravity;
    const auto begin=std::chrono::steady_clock::now();
    Check(!slow.Capture(View(cloud),frame),"queue failed to remain bounded");
    const auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
    Check(elapsed<.1,"full queue blocked algorithm callback");
    Check(slow.GetStats().dropped_queue==1,"queue loss not visible"); slow.Stop();
    Check(slow.GetStats().saved==2,"normal stop did not drain evidence");
}

void WaitSaved(R& recorder,std::uint64_t count) {
    const auto deadline=std::chrono::steady_clock::now()+std::chrono::seconds(5);
    while (recorder.GetStats().saved<count && std::chrono::steady_clock::now()<deadline)
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    Check(recorder.GetStats().saved==count,"asynchronous evidence save timed out");
}

void EpisodeAndUnavailable() {
    Fixture f; Cloud cloud{{0,0,0,0},{.1f,0,0,0}};
    auto context=Context(f.root); auto options=Options(); options.recovery_seconds=2;
    R recorder(options,context,View(cloud)); Check(recorder.Start(),"episode start");
    auto frame=Frame(cloud); frame.phase=R::Phase::kStartup; frame.stage=R::Stage::kSearch;
    frame.reason=R::Reason::kSearch; frame.candidate_available=false; frame.pose_available={{false,false,false}};
    frame.evaluated=false; frame.validation_points=0; frame.refinement_backend=R::Backend::kNone;
    std::strcpy(frame.rejection_reason.data(),"geometric_registration_failed");
    Check(recorder.Capture(View(cloud),frame),"first startup failure missing"); WaitSaved(recorder,1);
    for (int i=0;i<120;++i) { ++frame.sensor_stamp_ns; Check(!recorder.Capture(View(cloud),frame),"same startup reason produced periodic capture"); }
    auto paths=Events(context.directory); auto meta=YAML::LoadFile((paths[0]/"summary.json").string());
    Check(meta["schema"].as<int>()==2 && meta["poses"]["guess"].IsNull() && meta["poses"]["ndt"].IsNull() && meta["poses"]["check"].IsNull(),"unavailable poses fabricated");
    Check(fs::file_size(paths[0]/"map_local.bin")==0 && meta["map"]["coverage"].as<std::string>()=="unavailable_no_candidate","no-candidate map was fabricated");
    frame.phase=R::Phase::kGlobalRelocalization; frame.stage=R::Stage::kMapPrecheck;
    frame.reason=R::Reason::kOverlap; frame.candidate_available=true; frame.pose_available={{true,false,true}};
    frame.guess=frame.check=Pose(); frame.evaluated=true; frame.validation_stride=1; frame.validation_points=2;
    frame.overlap=.50; frame.retrieval_score=.6;
    std::strcpy(frame.rejection_reason.data(),"candidate_map_precheck_rejected");
    Check(recorder.Capture(View(cloud),frame),"changed reason not captured"); WaitSaved(recorder,2);
    frame.overlap=.54; Check(!recorder.Capture(View(cloud),frame),"minor score drift used evidence slot");
    frame.overlap=.56; Check(recorder.Capture(View(cloud),frame),"meaningful candidate improvement not captured"); WaitSaved(recorder,3);
    frame.reason=R::Reason::kBounds; frame.stage=R::Stage::kMapValidation;
    for (int i=0;i<120;++i) { ++frame.sensor_stamp_ns; Check(!recorder.Capture(View(cloud),frame),"episode exceeded three snapshots"); }
    // A single apparently successful candidate does not reset the quota.
    recorder.ObserveTracking(true); recorder.ObserveTracking(false);
    Check(!recorder.Capture(View(cloud),frame),"single success reset episode");
    // Two isolated successes separated by missing/invalid inputs must not be
    // mistaken for a continuously healthy two-second recovery interval.
    recorder.ObserveTracking(true);
    std::this_thread::sleep_for(std::chrono::milliseconds(2100));
    recorder.ObserveTracking(true); recorder.ObserveTracking(false);
    Check(!recorder.Capture(View(cloud),frame),"missing observation interval reset episode");
    // Repeated success notifications must measure from the first success, not
    // perpetually reset the start time on each frame.
    recorder.ObserveTracking(true);
    for (int i=0;i<11;++i) { std::this_thread::sleep_for(std::chrono::milliseconds(200)); recorder.ObserveTracking(true); }
    recorder.ObserveTracking(false);
    Check(recorder.Capture(View(cloud),frame),"continuous successful recovery did not reset episode"); recorder.Stop();
    Check(recorder.GetStats().saved==4,"new episode count wrong");
    paths=Events(context.directory); meta=YAML::LoadFile((paths.back()/"summary.json").string());
    Check(meta["episode"]["id"].as<int>()==2 && meta["episode"]["index"].as<int>()==1,"new episode metadata wrong");
    Check(meta["episode"]["recovery_seconds"].as<double>()==2 && meta["episode"]["recovery_max_gap_sec"].as<double>()==.5,"recovery continuity policy missing");
    const auto status=YAML::LoadFile((fs::path(context.directory)/"collector_status.json").string());
    Check(status["latest_observation"]["sensor_stamp_ns"].as<std::int64_t>()==frame.sensor_stamp_ns,"latest lightweight observation lost after cap");
    Check(status["episode"]["suppressed_at_limit"].as<int>()>=120,"episode limit not visible");
}

void FullGlobalSamples() {
    Fixture f; Cloud cloud;
    for (unsigned i=0;i<60001;++i) cloud.push_back({float(i%200)*.1f,float(i/200)*.1f,0,0});
    cloud[4].x=std::numeric_limits<float>::quiet_NaN();
    auto options=Options(); options.max_scan_points=4;
    auto frame=Frame(cloud); frame.phase=R::Phase::kGlobalRelocalization; frame.stage=R::Stage::kMapPrecheck;
    frame.validation_stride=2; frame.validation_points=30000;
    auto context=Context(f.root); R recorder(options,context,View(cloud)); Check(recorder.Start(),"global samples start");
    Check(recorder.Capture(View(cloud),frame),"global >2000 validation samples rejected"); recorder.Stop();
    auto paths=Events(context.directory); Check(paths.size()==1,"global validation package missing");
    const auto bytes=Read(paths[0]/"scan.bin"); std::vector<unsigned> actual,expected;
    for (std::size_t i=0;i<bytes.size();i+=24) if(U32(&bytes[i+20])) actual.push_back(U32(&bytes[i+16]));
    unsigned finite=0; for (unsigned i=0;i<cloud.size();++i) if(std::isfinite(cloud[i].x)) { if(finite%2==0)expected.push_back(i); ++finite; }
    Check(actual==expected && actual.size()==30000,"global exact finite-index sampling changed");
}

void ExportRelocalizationFixtures(const fs::path& root) {
    Cloud cloud{{0,0,0,0},{.1f,.1f,0,0},{0,.1f,0,0}};
    auto context=Context(root); R recorder(Options(),context,View(cloud)); Check(recorder.Start(),"relocalization fixture start");
    auto frame=Frame(cloud); frame.phase=R::Phase::kStartup; frame.stage=R::Stage::kSearch;
    frame.reason=R::Reason::kSearch; frame.evaluated=false; frame.validation_points=0;
    frame.candidate_available=false; frame.pose_available={{false,false,false}}; frame.refinement_backend=R::Backend::kNone;
    std::strcpy(frame.rejection_reason.data(),"geometric_registration_failed");
    Check(recorder.Capture(View(cloud),frame),"no-candidate fixture capture"); WaitSaved(recorder,1);
    frame=Frame(cloud); frame.phase=R::Phase::kGlobalRelocalization; frame.stage=R::Stage::kMapPrecheck;
    frame.reason=R::Reason::kBounds; frame.candidate_available=true; frame.candidate_id=42; frame.retrieval_score=.66;
    frame.pose_available={{true,false,true}}; frame.guess=frame.check=Pose(500); frame.refinement_backend=R::Backend::kNone;
    frame.overlap=frame.inside_xy=frame.inside_xyz=0;
    std::strcpy(frame.rejection_reason.data(),"candidate_map_precheck_rejected");
    Check(recorder.Capture(View(cloud),frame),"outside-map fixture capture"); recorder.Stop();
    Check(recorder.GetStats().saved==2,"relocalization fixture save");
    for (const auto& event:Events(context.directory)) std::cout<<event<<'\n';
}
void ValidationAndLock() {
    auto options=Options(); options.protect_days=6;
    bool rejected=false; try { options.Validate(); } catch (const std::invalid_argument&) { rejected=true; }
    Check(rejected,"protection shorter than 7 days accepted");
    Fixture f; Cloud cloud{{0,0,0,0}}; auto context=Context(f.root); options=Options();
    R first(options,context,View(cloud)), second(options,context,View(cloud));
    Check(first.Start() && !second.Start(),"shared archive writer lock failed");
    auto frame=Frame(cloud); ++frame.validation_points;
    Check(!first.Capture(View(cloud),frame),"incorrect validation point count accepted");
    Check(first.GetStats().dropped_oversize==1,"sampling mismatch not reported"); first.Stop();
}

void ZeroSamplesAndPartialQuota() {
    Fixture f; Cloud cloud{{0,0,0,0}}; auto context=Context(f.root); auto options=Options();
    R recorder(options,context,View(cloud)); Check(recorder.Start(),"zero sample start");
    auto frame=Frame(cloud); frame.validation_points=0; frame.overlap=0;
    Check(recorder.Capture(View(cloud),frame),"zero-sample failed criterion lost all evidence"); recorder.Stop();
    auto paths=Events(context.directory); Check(paths.size()==1,"zero sample event absent");
    auto meta=YAML::LoadFile((paths[0]/"summary.json").string());
    Check(!meta["criterion"]["sample_replay_available"].as<bool>() && meta["scan"]["validation_points"].as<int>()==0,"zero sample falsely exact");

    Fixture f2; context=Context(f2.root); options.event_mb=options.disk_mb=1;
    const auto partial=fs::path(context.directory)/".crashed.partial"; fs::create_directories(partial);
    { std::ofstream out(partial/"unknown.bin"); out.seekp(MiB-1); out.put('x'); }
    R limited(options,context,View(cloud)); Check(limited.Start(),"partial quota start");
    Check(limited.Capture(View(cloud),Frame(cloud)),"partial quota admission"); limited.Stop();
    Check(limited.GetStats().rejected_disk==1 && fs::exists(partial/"unknown.bin"),"partial bypassed quota or was deleted");
}

fs::path ExportFixture(const fs::path& root) {
    // Synthetic two static walls plus a central object absent from the map.
    Cloud map, scan;
    for (int side : {-1,1}) for (int x=0;x<40;++x) for (int z=0;z<4;++z)
        map.push_back({float(x)*.2f,float(side)*4.f,float(z)*.3f,0});
    scan=map;
    for (int x=0;x<12;++x) for (int y=0;y<12;++y)
        scan.push_back({2.f+float(x)*.1f,-.6f+float(y)*.1f,1.f,1});
    auto frame=Frame(scan); frame.guess=frame.ndt=frame.check=Pose();
    frame.overlap=double(map.size())/scan.size();
    auto context=Context(root); context.map_min={0,-4,0}; context.map_max={7.8,4,1.5};
    R recorder(Options(),context,View(map)); Check(recorder.Start(),"fixture start");
    Check(recorder.Capture(View(scan),frame),"fixture capture"); recorder.Stop();
    const auto events=Events(context.directory); Check(events.size()==1,"fixture save"); return events.front();
}
}  // namespace
int main(int argc,char** argv) {
    google::InitGoogleLogging(argv[0]);
    try {
        if (argc==3 && (std::string(argv[1])=="--output" || std::string(argv[1])=="--export-fixture" || std::string(argv[1])=="--export-relocalization-fixtures")) {
            const auto root=fs::absolute(argv[2]); fs::create_directories(root);
            Check(!fs::exists(root/"evidence"),"example output must be new");
            if (std::string(argv[1])=="--export-relocalization-fixtures") ExportRelocalizationFixtures(root);
            else std::cout<<ExportFixture(root)<<'\n';
            return 0;
        }
        Fixture f; ExactSelection(f.root); ProtectionAndRetention(); OversizeMapAndQueue(); ValidationAndLock(); ZeroSamplesAndPartialQuota(); EpisodeAndUnavailable(); FullGlobalSamples();
        std::cout<<"match_failure_recorder_test passed\n"; return 0;
    } catch (const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
