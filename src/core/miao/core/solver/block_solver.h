//
// Created by xiang on 24-4-26.
//

#ifndef MIAO_BLOCK_SOLVER_H
#define MIAO_BLOCK_SOLVER_H

/// 模板化的solver
/// Block solver 是 Solver 的子类，实际完成solver中的各种计算
/// 同时持有linear solver的指针，用linear solver计算Ax=GetB
/// NOTE: 在没有pose, landmark概念的问题中，尽量不要定义pose, landmark的概念

#include "core/common/config.h"
#include "core/sparse/sparse_block_matrix_diagonal.h"
#include "linear_solver.h"
#include "solver.h"

#include <memory>

namespace lightning::miao {

/**
 * \brief traits to summarize the properties of the fixed size optimization
 * problem
 */
template <int _PoseDim, int _LandmarkDim>
struct BlockSolverTraits {
    static const int PoseDim = _PoseDim;  ///< 位姿块的维度，动态配置时为 Eigen::Dynamic。
    static const int LandmarkDim = _LandmarkDim;  ///< 路标块的维度，动态配置时为 Eigen::Dynamic。
    typedef Eigen::Matrix<double, PoseDim, PoseDim> PoseMatrixType;
    typedef Eigen::Matrix<double, LandmarkDim, LandmarkDim> LandmarkMatrixType;
    typedef Eigen::Matrix<double, PoseDim, LandmarkDim> PoseLandmarkMatrixType;
    typedef Eigen::Matrix<double, PoseDim, 1> PoseVectorType;
    typedef Eigen::Matrix<double, LandmarkDim, 1> LandmarkVectorType;

    typedef SparseBlockMatrix<PoseMatrixType> PoseHessianType;
    typedef SparseBlockMatrix<LandmarkMatrixType> LandmarkHessianType;
    typedef SparseBlockMatrix<PoseLandmarkMatrixType> PoseLandmarkHessianType;
    typedef LinearSolver<PoseMatrixType> LinearSolverType;
};

/// Block Solver X 的特化
/**
 * \brief traits to summarize the properties of the dynamic size optimization
 * problem
 */
template <>
struct BlockSolverTraits<Eigen::Dynamic, Eigen::Dynamic> {
    static const int PoseDim = Eigen::Dynamic;  ///< 位姿块的维度，动态配置时为 Eigen::Dynamic。
    static const int LandmarkDim = Eigen::Dynamic;  ///< 路标块的维度，动态配置时为 Eigen::Dynamic。
    typedef lightning::MatrixX PoseMatrixType;
    typedef lightning::MatrixX LandmarkMatrixType;
    typedef lightning::MatrixX PoseLandmarkMatrixType;
    typedef lightning::VectorX PoseVectorType;
    typedef lightning::VectorX LandmarkVectorType;

    typedef SparseBlockMatrix<PoseMatrixType> PoseHessianType;
    typedef SparseBlockMatrix<LandmarkMatrixType> LandmarkHessianType;
    typedef SparseBlockMatrix<PoseLandmarkMatrixType> PoseLandmarkHessianType;
    typedef LinearSolver<PoseMatrixType> LinearSolverType;
};

/**
 * \brief Implementation of a Solver operating on the blocks of the Hessian
 *
 * 带有pose和landmark的概念，默认Landmark会被marg，而pose不会
 * 在没有pose, landmark概念的问题中，可以简单地以是否marg来区分pose或landmark
 */
template <typename Traits>
class BlockSolver : public Solver {
   public:
    static const int PoseDim = Traits::PoseDim;  ///< 位姿块的维度，动态配置时为 Eigen::Dynamic。
    static const int LandmarkDim = Traits::LandmarkDim;  ///< 路标块的维度，动态配置时为 Eigen::Dynamic。
    typedef typename Traits::PoseMatrixType PoseMatrixType;
    typedef typename Traits::LandmarkMatrixType LandmarkMatrixType;
    typedef typename Traits::PoseLandmarkMatrixType PoseLandmarkMatrixType;
    typedef typename Traits::PoseVectorType PoseVectorType;
    typedef typename Traits::LandmarkVectorType LandmarkVectorType;

    typedef typename Traits::PoseHessianType PoseHessianType;
    typedef typename Traits::LandmarkHessianType LandmarkHessianType;
    typedef typename Traits::PoseLandmarkHessianType PoseLandmarkHessianType;
    typedef typename Traits::LinearSolverType LinearSolverType;

   public:
    /**
     * allocate a Block solver ontop of the underlying linear Solver.
     * NOTE: The BlockSolver assumes exclusive access to the linear Solver and
     * will therefore free the pointer in its destructor.
     */
    explicit BlockSolver(std::unique_ptr<LinearSolverType> linearSolver);
    ~BlockSolver() override;

    bool Init(Optimizer* optimizer) override;

    /// 初始化时调用，建立各矩阵块的数据
    bool BuildStructure(bool zeroBlocks) override;

    /// 原生模式下建立问题结构
    bool BuildStructureFromRaw(bool zero_blocks);

    /// 增量模式下建立问题结构
    bool BuildStructureInc(bool zero_blocks);

    /// 每次迭代时调用的函数,计算当前这次迭代中的jacobian
    bool BuildSystem() override;

    /// 建立H, b矩阵并交给linear solver
    bool Solve() override;

    bool SetLambda(double lambda, bool backup) override;
    void RestoreDiagonal() override;
    bool SupportsSchur() override { return true; }
    bool Schur() override { return do_schur_; }
    void SetSchur(bool s) override { do_schur_ = s; }

    void MultiplyHessian(double* dest, const double* src) const override {
        Hpp_->MultiplySymmetricUpperTriangle(dest, src);
    }

   protected:
    void Resize(const std::vector<int>& blockPoseIndices, const std::vector<int>& blockLandmarkIndices, int totalDim);

    void Deallocate();

    OptimizerConfig config_;  ///< 块求解器的运行配置。

    /// pose to pose 的矩阵块，解稠密问题，只求解这一部分的计算
    std::unique_ptr<SparseBlockMatrix<PoseMatrixType>> Hpp_;  // pose to pose

    std::unique_ptr<SparseBlockMatrix<LandmarkMatrixType>> Hll_;      ///< landmark to landmark
    std::unique_ptr<SparseBlockMatrix<PoseLandmarkMatrixType>> Hpl_;  ///< pose to landmark

    std::unique_ptr<SparseBlockMatrix<PoseMatrixType>> Hschur_;  ///< 消去路标后构建的位姿 Schur 补矩阵。
    std::unique_ptr<SparseBlockMatrixDiagonal<LandmarkMatrixType>> D_inv_schur_;  ///< 路标对角 Hessian 块的逆，用于 Schur 消元。

    std::unique_ptr<SparseBlockMatrixCCS<PoseLandmarkMatrixType>> Hpl_CCS_;  ///< 位姿与路标交叉块的压缩列索引。
    std::unique_ptr<SparseBlockMatrixCCS<PoseMatrixType>> H_schur_transpose_CCS_;  ///< 构建 Schur 补使用的转置压缩列索引。

    std::unique_ptr<LinearSolverType> linear_solver_;  ///< 实际求解块正规方程的线性求解器。

    std::vector<PoseVectorType> diagonal_backup_pose_;  ///< 施加阻尼前备份的位姿对角块。
    std::vector<LandmarkVectorType> diagonal_backup_landmark_;  ///< 施加阻尼前备份的路标对角块。

    bool do_schur_ = true;  ///< 是否使用 Schur 消元。

    lightning::VecXd coeffs_;  ///< 求解过程使用的临时系数向量。
    lightning::VecXd b_schur_;  ///< Schur 消元后的位姿右端项。

    int num_poses_ = 0 /**< 参与求解的位姿顶点数量。 */, num_landmarks_ = 0 /**< 参与求解的路标顶点数量。 */;
    int size_poses_ = 0 /**< 位姿变量的总向量维度。 */, size_landmarks_ = 0 /**< 路标变量的总向量维度。 */;
};

template <int p, int l>
using BlockSolverPL = BlockSolver<BlockSolverTraits<p, l>>;

}  // namespace lightning::miao

#include "block_solver.hpp"

#endif  // MIAO_BLOCK_SOLVER_H
