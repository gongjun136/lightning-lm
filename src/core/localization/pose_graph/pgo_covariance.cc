#include "core/localization/pose_graph/pgo_impl.h"
#include "core/localization/pose_covariance.h"
#include "core/robust_kernel/robust_kernel.h"

namespace lightning::loc {
void PGOImpl::ComputePosteriorCovariance() {
    auto& latest = frames_.back();
    latest->pose_covariance_valid_ = false;
    latest->pose_unobservable_ = false;
    std::map<miao::Vertex*, int> offsets;
    for (size_t i = 0; i < frames_.size(); ++i) {
        auto vertex = optimizer_->GetVertex(frames_[i]->frame_id_);
        if (!vertex) return;
        offsets[vertex.get()] = 6 * i;
    }
    Eigen::MatrixXd h = Eigen::MatrixXd::Zero(6 * frames_.size(), 6 * frames_.size());
    // Linearize the physical factors at the converged poses, outside the LM solver.
    // Robust rho' preserves PSD and suppresses rejected measurements.
    for (auto* edge : optimizer_->ActiveEdges()) {
        if (edge->Level() != 0) continue;
        auto* typed = dynamic_cast<miao::BaseEdge<6, SE3>*>(edge);
        if (!typed) return;
        edge->ComputeError();
        const auto residual = typed->GetError();
        Mat6d information = typed->Information();
        if (!residual.allFinite() || !information.allFinite()) return;
        if (edge->GetRobustKernel()) {
            Vec3d rho;
            edge->GetRobustKernel()->Robustify(edge->Chi2(), rho);
            if (!rho.allFinite() || rho[1] < 0.0) return;
            information *= rho[1];
        }
        std::vector<Mat6d> jacobians;
        for (auto* raw : edge->GetVertices()) {
            auto* vertex = dynamic_cast<miao::VertexSE3*>(raw);
            if (!vertex || offsets.count(raw) == 0) return;
            const SE3 saved = vertex->Estimate();
            Mat6d j;
            for (int k = 0; k < 6; ++k) {
                Vec6d delta = Vec6d::Zero();
                delta[k] = 1e-6;
                vertex->SetEstimate(saved);
                vertex->OplusImpl(delta.data());
                edge->ComputeError();
                const Vec6d plus = typed->GetError();
                vertex->SetEstimate(saved);
                delta[k] = -1e-6;
                vertex->OplusImpl(delta.data());
                edge->ComputeError();
                j.col(k) = (plus - typed->GetError()) / 2e-6;
            }
            vertex->SetEstimate(saved);
            if (!j.allFinite()) { edge->ComputeError(); return; }
            jacobians.push_back(j);
        }
        edge->ComputeError();
        for (size_t i = 0; i < jacobians.size(); ++i)
            for (size_t j = 0; j < jacobians.size(); ++j)
                h.block<6, 6>(offsets.at(edge->GetVertex(i)), offsets.at(edge->GetVertex(j))) +=
                    jacobians[i].transpose() * information * jacobians[j];
    }
    latest->pose_covariance_valid_ = MarginalPoseCovariance(h, latest->pose_covariance_);
    latest->pose_unobservable_ = h.allFinite() && !latest->pose_covariance_valid_;
}
}
