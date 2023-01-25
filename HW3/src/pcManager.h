//
// Created by bettaz on 18/01/23.
//

#ifndef THIRD_ASSIGNMENT_PCMANAGER_H
#define THIRD_ASSIGNMENT_PCMANAGER_H
#include <random>
#include <array>
#include "Eigen/Dense"
#include "nanoflann.hpp"

// Adapter from OpenGP because the nanoflann works well with horizontal vectors making the code ugly
namespace nanoflann {
    /// KD-tree adaptor for working with data directly stored in an Eigen Matrix, without duplicating the data storage.
    /// This code is adapted from the KDTreeEigenMatrixAdaptor class of nanoflann.hpp
    template <class MatrixType, int DIM = -1, class Distance = nanoflann::metric_L2, typename IndexType = int>
    struct KDTreeAdaptor {
        typedef KDTreeAdaptor<MatrixType,DIM,Distance> self_t;
        typedef typename MatrixType::Scalar              num_t;
        typedef typename Distance::template traits<num_t,self_t>::distance_t metric_t;
        typedef KDTreeSingleIndexAdaptor< metric_t,self_t,DIM,IndexType>  index_t;
        index_t* index;
        KDTreeAdaptor(const MatrixType &mat, const int leaf_max_size = 10) : m_data_matrix(mat) {
            const size_t dims = mat.rows();
            index = new index_t( dims, *this, nanoflann::KDTreeSingleIndexAdaptorParams(leaf_max_size) );
            index->buildIndex();
        }
        ~KDTreeAdaptor() {delete index;}
        const MatrixType &m_data_matrix;
        /// Query for the num_closest closest points to a given point (entered as query_point[0:dim-1]).
        inline void query(const num_t *query_point, const size_t num_closest, IndexType *out_indices, num_t *out_distances_sq) const {
            nanoflann::KNNResultSet<typename MatrixType::Scalar,IndexType> resultSet(num_closest);
            resultSet.init(out_indices, out_distances_sq);
            index->findNeighbors(resultSet, query_point, nanoflann::SearchParameters());
        }
        /// Query for the closest points to a given point (entered as query_point[0:dim-1]).
        inline IndexType closest(const num_t *query_point) const {
            IndexType out_indices;
            num_t out_distances_sq;
            query(query_point, 1, &out_indices, &out_distances_sq);
            return out_indices;
        }
        const self_t & derived() const {return *this;}
        self_t & derived() {return *this;}
        inline size_t kdtree_get_point_count() const {return m_data_matrix.cols();}
        /// Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
        inline num_t kdtree_distance(const num_t *p1, const size_t idx_p2,size_t size) const {
            num_t s=0;
            for (size_t i=0; i<size; i++) {
                const num_t d= p1[i]-m_data_matrix.coeff(i,idx_p2);
                s+=d*d;
            }
            return s;
        }
        /// Returns the dim'th component of the idx'th point in the class:
        inline num_t kdtree_get_pt(const size_t idx, int dim) const {
            return m_data_matrix.coeff(dim,idx);
        }
        /// Optional bounding-box computation: return false to default to a standard bbox computation loop.
        template <class BBOX> bool kdtree_get_bbox(BBOX&) const {return false;}
    };
}

namespace pc{
    using PointCloud = Eigen::Matrix3Xd;
    using happlyPC = std::vector<std::array<double,3>>;
    using Point = Eigen::Vector3d;
    using KdTree = nanoflann::KDTreeAdaptor<Eigen::Matrix3Xd, 3, nanoflann::metric_L2_Simple>;
    using Translation = Eigen::Vector3d;
    using Affine = Eigen::Affine3d;


    inline PointCloud
    vecToEigen(const happlyPC &src){
        PointCloud pcOut(3,src.size());
        for (size_t i = 0; i < src.size(); i++)
            pcOut.col(i) = Point(src[i].data());
        return pcOut;
    };


    inline happlyPC
    EigenToVec(const PointCloud &src){
        happlyPC pcOut;
        for (auto line : src.colwise()) {
            pcOut.push_back({{line(0),line(1),line(2)}});
        }
        return pcOut;
    };


    inline PointCloud
    addNoise(const PointCloud &src, const double &mean, const double & stdDev) {
        PointCloud noisy(src);
        if(stdDev != 0.0){
            std::default_random_engine generator;
            std::normal_distribution<double> dist(mean, stdDev);
            for(auto it : noisy.colwise()){
                it = it + Point(dist(generator),dist(generator),dist(generator));
            }
        }
        return noisy;
    };


    inline PointCloud
    rotoTranslatePointCloud(const PointCloud &src, const Affine & rotoTranslation) {
        return rotoTranslation*src;
    };


    inline double pointDistance(const Point & p1, const Point &p2){
        return sqrt((p1-p2).norm());
    }


    inline Affine
    interpolate(PointCloud &X, PointCloud &Q ) {
        /// Calculates the centroids and centers both PCs
        Point X_mean = X.rowwise().sum()/X.cols();
        Point Q_mean = Q.rowwise().sum()/Q.cols();
        X.colwise() -= X_mean;
        Q.colwise() -= Q_mean;
        /// Compute transformation through the SVD method
        Affine transformation;
        Eigen::Matrix3d covariance;
        if(X.cols()==Q.cols()){
            covariance= X * Q.transpose();
        }
        else{
            X.block(0,0,3,Q.cols())*Q.transpose();
        }

        Eigen::JacobiSVD<Eigen::Matrix3d> svd(covariance, Eigen::ComputeFullU | Eigen::ComputeFullV);
        if(svd.matrixU().determinant()*svd.matrixV().determinant() < 0.0) {
            Eigen::Vector3d S = Eigen::Vector3d::Ones(); S(2) = -1.0;
            transformation.linear().noalias() = svd.matrixV()*S.asDiagonal()*svd.matrixU().transpose();
        } else {
            transformation.linear().noalias() = svd.matrixV()*svd.matrixU().transpose();
        }
        transformation.translation().noalias() = Q_mean - transformation.linear()*X_mean;
        /// Apply transformation
        X = transformation*X;
        /// Restores the previous position
        X.colwise() += X_mean;
        Q.colwise() += Q_mean;
        return transformation;
    }


    inline double computeError(const PointCloud &base, const PointCloud &reference) {
        /// Computes distance between every point assuming they are ordered by it
        PointCloud diff = base - reference;
        return diff.cwiseAbs2().colwise().sum().cwiseSqrt().mean();
    }


    inline PointCloud
    ICP(const PointCloud &reference, const PointCloud &source, Affine &estRotTransformation, int &iterationCount, double & finalError, const double &epsilon=0.5, const size_t &maxIteration=1000) {
        std::cout<<"Started ICP"<<std::endl;
        /// Temporary storage point-clouds and variables
        PointCloud ref(reference),src(source);
        Eigen::Matrix3Xd toInterpolate = Eigen::Matrix3Xd::Zero(3, ref.cols());
        Eigen::Matrix3Xd outPC = ref;
        Affine iterationRotoTranslation;
        estRotTransformation.linear() = Eigen::Matrix3d::Identity();
        estRotTransformation.translation() = Point(0.0,0.0,0.0);
        double prevIterationError= epsilon+1.0;
        finalError = prevIterationError;
        /// Generates a KdTree from the source point-cloud
        KdTree kdTree(src);
        int stopping = 0;
        /// ICP
        for(int iteration=0; iteration < maxIteration; ++iteration) {
            std::cout << "Iteration " << iteration << " / " << maxIteration << std::endl;
            /// Find the closest point
            for(int i=0; i < ref.cols(); ++i) {
                toInterpolate.col(i) = src.col(kdTree.closest(ref.col(i).data()));
            }
            /// Stopping criteria
            /// Compute rotation and translation
            iterationRotoTranslation = interpolate(ref, toInterpolate);
            finalError = computeError(ref, src);
            std::cout << "Error " << iteration << ": " << finalError << std::endl;
            estRotTransformation = iterationRotoTranslation * estRotTransformation;
            outPC = ref;
            /// Checks convergence
            if(prevIterationError == finalError || finalError < epsilon){
                break;
            }
            prevIterationError = finalError;
            iterationCount = iteration;
            toInterpolate.setZero();
        }
        std::cout<<"End ICP"<<std::endl;
        return outPC;
    };

    /// Custom PointCloud based on the Eigen::Matrix3Xd that gathers each node information to trim the pointcloud
    struct DistanceResult{
        double distance;
        size_t index_s;
        size_t index_d;
    }Distance;
    struct
    {
        bool operator()(DistanceResult a, DistanceResult b) const { return a.distance < b.distance; }
    }DistanceCriteria;


    inline PointCloud trim(PointCloud & ref, const PointCloud & src, std::vector<DistanceResult> &distances, const float &xsi){
        std::sort(distances.begin(), distances.end(),DistanceCriteria);
        const int total = distances.size();
        const int index = (xsi * (float)total)-1;
        distances.erase(distances.begin()+index,distances.end());
        PointCloud newRef(3, total);
        PointCloud outSrc(3, total);
        for(size_t i = 0; i<total; i++){
            newRef.col(i) = ref.col(distances[i].index_d);
            outSrc.col(i) = src.col(distances[i].index_s);
        }
        ref = newRef;
        return outSrc;
    }


    inline PointCloud
    trICP(const PointCloud &reference, const PointCloud &source, Affine &estRotTransformation, int &iterationCount, double & finalError, const double &epsilon=0.5, const float &xsi = 0.8, const size_t &maxIteration=1000){
        std::cout<<"Started Tr-ICP"<<std::endl;
        /// Temporary storage point-clouds and variables
        PointCloud ref(reference),src(source);
        PointCloud toInterpolate = Eigen::Matrix3Xd::Zero(3, ref.cols());
        PointCloud outPC = ref;
        Affine iterationRotoTranslation;
        estRotTransformation.linear() = Eigen::Matrix3d::Identity();
        estRotTransformation.translation() = Point(0.0,0.0,0.0);
        double prevIterationError= epsilon+1.0;
        finalError = prevIterationError;
        /// Generates a KdTree from the source point-cloud
        KdTree kdTree(src);
        /// ICP
        for(int iteration=0; iteration < maxIteration; ++iteration) {
            std::cout << "Iteration " << iteration << " / " << maxIteration << std::endl;
            toInterpolate = PointCloud::Zero(3, ref.cols());
            /// Find the closest point and stores data for trimming
            std::vector<DistanceResult> dist;
            for(size_t i=0; i < ref.cols(); ++i) {
                size_t closest = kdTree.closest(ref.col(i).data());
                dist.push_back(DistanceResult({.distance=pointDistance(ref.col(i), src.col(closest)), .index_s=i, .index_d=static_cast<size_t>(closest)}));
            }
            /// Trimming and interpolation process
            PointCloud trimRef(ref);
            toInterpolate = trim(trimRef, src , dist, xsi);
            /// Compute rotation and translation
            iterationRotoTranslation = interpolate(trimRef, toInterpolate);
            estRotTransformation = iterationRotoTranslation * estRotTransformation;
            /// Apply transform to the original reference PC
            ref = iterationRotoTranslation * ref;
            outPC = ref;
            /// Stopping criteria
            finalError = computeError(ref, src);
            std::cout << "Error " << iteration << ": " << finalError << std::endl;
            if(prevIterationError == finalError || finalError < epsilon || prevIterationError+epsilon< finalError){
                break;
            }
            prevIterationError = finalError;
            iterationCount = iteration;
            dist.empty();
        }
        std::cout<<"End Tr-ICP!"<<std::endl;
        return outPC;
    }
}
#endif //THIRD_ASSIGNMENT_PCMANAGER_H
