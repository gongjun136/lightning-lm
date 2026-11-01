//
// Created by gx on 23-11-23.
//

#include "utils/timer.h"

#include <algorithm>
#include <fstream>
#include <numeric>

namespace lightning {

std::map<std::string, Timer::TimerRecord> Timer::records_;

void Timer::PrintAll() {
    LOG(INFO) << (">>> ===== 运行耗时统计 =====");
    for (auto& r : records_) {
        auto& rec = r.second.time_usage_in_ms_;
        std::sort(rec.begin(), rec.end());

        LOG(INFO) << "耗时统计 > [" << r.first
                  << "] average time usage: " << std::accumulate(rec.begin(), rec.end(), 0.0) / double(rec.size())
                  << " ms, med: " << rec[rec.size() * 0.5] << " 95%: " << rec[rec.size() * 0.95]
                  << ", called times : " << r.second.time_usage_in_ms_.size();
    }
    LOG(INFO) << ">>> ===== 运行耗时统计结束 =====";
}

void Timer::DumpIntoFile(const std::string& file_name) {
    std::ofstream ofs(file_name, std::ios::out);
    if (!ofs.is_open()) {
        LOG(ERROR) << "event=TIMER_FILE_OPEN_FAILED 无法打开计时文件 path=" << file_name;
        return;
    } else {
        LOG(INFO) << "event=TIMER_RECORDS_SAVED 计时记录已写入 path=" << file_name;
    }

    size_t max_length = 0;
    for (const auto& iter : records_) {
        ofs << iter.first << ", ";
        if (iter.second.time_usage_in_ms_.size() > max_length) {
            max_length = iter.second.time_usage_in_ms_.size();
        }
    }
    ofs << std::endl;

    for (size_t i = 0; i < max_length; ++i) {
        for (const auto& iter : records_) {
            if (i < iter.second.time_usage_in_ms_.size()) {
                ofs << iter.second.time_usage_in_ms_[i] << ",";
            } else {
                ofs << ",";
            }
        }
        ofs << std::endl;
    }
    ofs.close();
}

double Timer::GetMeanTime(const std::string& func_name) {
    if (records_.find(func_name) == records_.end()) {
        return 0.0;
    }

    auto r = records_[func_name];
    return std::accumulate(r.time_usage_in_ms_.begin(), r.time_usage_in_ms_.end(), 0.0) /
           double(r.time_usage_in_ms_.size());
}
}  // namespace lightning
