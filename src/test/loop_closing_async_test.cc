#include <algorithm>
#include <chrono>
#include <condition_variable>
#include <future>
#include <iostream>
#include <memory>
#include <mutex>
#include <set>
#include <vector>

#include "core/loop_closing/loop_closing.h"

namespace {

// Expose the existing protected queue to control callback timing without NDT work.
class TestLoopClosing : public lightning::LoopClosing {
   public:
    using LoopClosing::LoopClosing;
    using LoopClosing::kf_thread_;
    using LoopClosing::optimizer_;
};

bool CheckQueueAndDestruction(int queued_frames) {
    using namespace std::chrono_literals;
    std::mutex mutex;
    std::condition_variable cv;
    bool first_started = false;
    bool release_first = false;
    std::vector<unsigned long> processed;

    lightning::LoopClosing::Options options;
    options.online_mode_ = true;
    auto loop = std::make_unique<TestLoopClosing>(options);
    loop->Init("");
    // Keep Init's capacity, replacing only the callback after joining its worker.
    loop->kf_thread_.Quit();
    loop->kf_thread_.SetProcFunc([&](const lightning::Keyframe::Ptr& kf) {
        if (kf->GetID() == 0) {
            std::unique_lock<std::mutex> lock(mutex);
            first_started = true;
            cv.notify_all();
            cv.wait(lock, [&] { return release_first; });
        }
        processed.push_back(kf->GetID());
    });
    loop->kf_thread_.Start();
    auto add_frame = [&](int id) {
        loop->AddKF(std::make_shared<lightning::Keyframe>(id, nullptr, lightning::NavState{}));
    };
    add_frame(0);
    {
        std::unique_lock<std::mutex> lock(mutex);
        if (!cv.wait_for(lock, 2s, [&] { return first_started; })) {
            release_first = true;
            cv.notify_all();
            std::cerr << "loop closing worker did not start" << std::endl;
            return false;
        }
    }
    for (int id = 1; id <= queued_frames; ++id) add_frame(id);
    const size_t retained = std::min(queued_frames, 40);
    const size_t dropped = queued_frames - retained;
    const bool counts_ok = loop->kf_thread_.PendingCount() == retained + 1 &&
                           loop->kf_thread_.DroppedCount() == dropped;

    // The owner must wait for the in-flight callback and drain the retained tail.
    auto closed = std::async(std::launch::async, [owner = std::move(loop)]() mutable { owner.reset(); });
    const bool waited = closed.wait_for(20ms) == std::future_status::timeout;
    {
        std::lock_guard<std::mutex> lock(mutex);
        release_first = true;
    }
    cv.notify_all();
    closed.get();

    std::vector<unsigned long> expected{0};
    for (int id = static_cast<int>(dropped) + 1; id <= queued_frames; ++id) expected.push_back(id);
    if (!counts_ok || !waited || processed != expected) {
        std::cerr << "loop closing queue/destructor contract failed for " << queued_frames
                  << " queued frames: counts_ok=" << counts_ok << ", waited=" << waited
                  << ", processed=" << processed.size() << ", expected=" << expected.size() << std::endl;
        return false;
    }
    return true;
}

bool CheckKeyframeIdGaps() {
    lightning::LoopClosing::Options options;
    options.online_mode_ = true;
    options.verbose_ = false;
    options.with_height_ = false;
    TestLoopClosing loop(options);
    loop.Init("");
    // These are the actual IDs retained after earlier frames were evicted.
    // ID 80 also enters the candidate/submap lookup path with an empty cloud.
    for (int id : {0, 11, 12, 80}) {
        loop.AddKF(std::make_shared<lightning::Keyframe>(
            id, std::make_shared<lightning::PointCloudType>(), lightning::NavState{}));
    }
    loop.kf_thread_.Quit();
    std::set<std::pair<int, int>> edges;
    for (const auto& edge : loop.optimizer_->GetEdges()) {
        edges.emplace(edge->GetVertex(0)->GetId(), edge->GetVertex(1)->GetId());
    }
    const std::set<std::pair<int, int>> expected{{0, 11}, {0, 12}, {11, 12}, {11, 80}, {12, 80}};
    if (loop.optimizer_->GetVertices().size() != 4 || edges != expected) {
        std::cerr << "loop closing lost odometry constraints across evicted keyframe IDs" << std::endl;
        return false;
    }
    return true;
}

}  // namespace

int main() {
    // Both the non-overloaded tail and the newest 40 of 50 waiting frames drain.
    if (!CheckQueueAndDestruction(10) || !CheckQueueAndDestruction(50) || !CheckKeyframeIdGaps()) return 1;
    return 0;
}
