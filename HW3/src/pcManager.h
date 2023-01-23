//
// Created by bettaz on 18/01/23.
//

#ifndef THIRD_ASSIGNMENT_PCMANAGER_H
#define THIRD_ASSIGNMENT_PCMANAGER_H
#include <random>
#include <array>
#include "Eigen/Dense"
#include "nanoflann.hpp"

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
    using Point = Eigen::Vector3d;
    using KdTree = nanoflann::KDTreeAdaptor<Eigen::Matrix<double, 3, -1, 0, 3, -1>>;
    using Translation = Eigen::Vector3d;
    using Affine = Eigen::Affine3d;


    inline PointCloud
    vecToEigen(const std::vector<std::array<double, 3>> &src){
        PointCloud pcOut(3,src.size());
        for (size_t i = 0; i < src.size(); i++)
            pcOut.col(i) = Point(src[i].data());
        return pcOut;
    };


    inline std::vector<std::array<double, 3>>
    EigenToVec(const PointCloud &src){
        std::vector<std::array<double, 3>> pcOut;
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
                it = it + Eigen::Vector3d(dist(generator),dist(generator),dist(generator));
            }
        }
        return noisy;
    };


    inline PointCloud
    rotoTranslatePointCloud(const PointCloud &src, const Affine & rotoTranslation) {
        return rotoTranslation*src;
    };

    inline double pointDistance(const Point & p1, const Point &p2){
        return sqrt((p1-p2).cwiseAbs2().sum());
    }


    inline Affine
    interpolate(PointCloud reference, PointCloud source) {
        Eigen::Matrix3d covariance;
        Point center1, center2;

        size_t pc1Count = reference.cols(), pc2Count = source.cols();
        // Calculate the center of pcl1
        //std::cout<<"Ref:"<<std::endl<<reference<<std::endl;
        center1 = reference.rowwise().mean();
        //std::cout<<"Cen:"<<std::endl<<center1<<std::endl;

        // Calculate the center of pcl2
        //std::cout<<"Src:"<<std::endl<<source<<std::endl;
        center2 = source.rowwise().mean();
        //std::cout<<"Cen:"<<std::endl<<center2<<std::endl;

        // Move pcl1 to the center
        //reference = rotoTranslatePointCloud(reference, Rotation(Eigen::Matrix3d::Identity()), -center1);
        reference.colwise()-=center1;
        //std::cout<<"RotRef:"<<std::endl<<reference<<std::endl;
        // Move pcl2 to the center

        //source = rotoTranslatePointCloud(source, Rotation(Eigen::Matrix3d::Identity()), -center2);
        source.colwise() -= center2;
        //std::cout<<"Src:"<<std::endl<<source<<std::endl;

        // Calculate covariance matrix
        /*for (int i = 0; i < pc2Count; i++) {
            covariance += source.col(i) * reference.col(i).transpose();
        }*/
        covariance = reference * source.transpose();
        //covariance /= pc2Count;
        Eigen::MatrixXd u, v;
        Eigen::JacobiSVD<Eigen::MatrixXd> jacobi(covariance, Eigen::ComputeFullU | Eigen::ComputeFullV);
        jacobi.compute(covariance);
        u=jacobi.matrixU();
        v=jacobi.matrixV();
        Affine rotoTransform;
        if(u.determinant()*v.determinant() < 0.0) {
            Eigen::Vector3d S = Eigen::Vector3d::Ones(); S(2) = -1.0;
            rotoTransform.linear().noalias() = v * S.asDiagonal() * u.transpose();
        } else {
            rotoTransform.linear().noalias() = v * u.transpose();
        }
        //rotoTransform.translation().noalias() = Y_mean - rotoTransform.linear()*X_mean;

        // Calculate the rotation matrix
        /*Eigen::Matrix3d rot = v * u.transpose();
        if (rot.determinant() < 0.0) {
            v.row(2) *= -1.0;
            rot = v * u.transpose();
        }*/

        // Calculate the translation matrix
        rotoTransform.translation()= center2 - (rotoTransform.linear() * center1);
        return rotoTransform;
    }

    double computeError(const PointCloud &base, const PointCloud &reference) {
        PointCloud diff = base - reference;
        return diff.cwiseAbs2().colwise().sum().cwiseSqrt().sum();;
    }

    inline PointCloud
    ICP(PointCloud &reference, PointCloud &source, Affine &estRotTras, int &iterationCount, const double &epsilon=5, const size_t &maxIteration=100) {
        std::cout<<"Started ICP"<<std::endl;
        double error = epsilon+1.0;
        PointCloud out(source);
        PointCloud ref(reference);
        estRotTras = Affine();
        for (int iter = 0; iter < maxIteration && error>epsilon; ++iter) {
            // Nanoflann KD-Tree
            KdTree tree(out);
            Eigen::Index closest;
            PointCloud toInterpolate;
            double distance;
            for (int i=0; i < ref.cols(); i++){
                closest = tree.closest(ref.col(i).data());
                distance = pointDistance(ref.col(i), out.col(closest));
                toInterpolate.col(i) = out.col(closest);
            }
            Affine iterRotoTranslation = interpolate(reference, toInterpolate);
            source = rotoTranslatePointCloud(source, iterRotoTranslation);

            toInterpolate.colwise() += toInterpolate.rowwise().mean();
            reference.colwise() += reference.rowwise().mean();

            happly::PLYData iterPLY;
            std::vector<std::array<double,3>> data = pc::EigenToVec(out);
            std::string fn = "./data/iteration" + std::to_string(iter) + ".ply";
            iterPLY.addVertexPositions(data);
            iterPLY.write(fn, happly::DataFormat::ASCII);

            error = computeError(reference, out);
            std::cout<<"Error it "<<iter<<": "<<error<<std::endl;
            estRotTras = iterRotoTranslation * estRotTras;
            iterationCount = iter;
        }
        std::cout<<"Finished ICP!"<<std::endl;
        return out;
    };

    struct DistanceResult{
        double distance;
        size_t index;
    }Distance;
    struct
    {
        bool operator()(DistanceResult a, DistanceResult b) const { return a.distance < b.distance; }
    }DistanceCriteria;
    class TrPointCloud : public PointCloud {
    public:
        std::vector<DistanceResult> distances;
        TrPointCloud(): PointCloud(){
            distances = std::vector<DistanceResult>();
        };
        TrPointCloud(PointCloud p): PointCloud(p){
            distances = std::vector<DistanceResult>(p.rows());
        };
        TrPointCloud(const int &x, const int &y) : PointCloud(x,y){
            distances = std::vector<DistanceResult>(x);
        }
        TrPointCloud getTrimmedPointCloud(const float &xsi){
            std::sort(distances.begin(), distances.end(),DistanceCriteria);
            const int total = distances.size();
            const int index = xsi * (float)total;
            int removeCount = 0;
            distances.erase(distances.begin()+index,distances.end());
            TrPointCloud out(3,index);
            for (auto entry : distances) {
                out.col(removeCount) = this->col(entry.index);
                removeCount++;
            }
            out.distances = this->distances;
            return out;
        }
    };

    inline PointCloud
    trICP( PointCloud &reference, PointCloud &source, Affine & estRotTran, int &iterationCount, const float & xsi = 0.8, const double &epsilon=0.001, const size_t &maxIteration=100){
        std::cout<<"Started Tr-ICP"<<std::endl;
        double error = epsilon+1.0;
        TrPointCloud out(source);
        estRotTran = Affine ();
        for (int iter = 0; iter < maxIteration && error>epsilon; ++iter) {
            Affine iterRotTran;
            // construct a kd-tree index:
            KdTree tree(out);
            TrPointCloud toInterpolate (3 ,reference.cols());
            Eigen::Index closer;
            double distance;
            for (int i=0; i<reference.cols(); i++){
                closer = tree.closest(reference.col(i).data());
                toInterpolate.distances[i] = DistanceResult({.distance=pointDistance(reference.col(i),out.col(closer)), .index=static_cast<size_t>(closer)});
            }
            toInterpolate= toInterpolate.getTrimmedPointCloud(xsi);
            iterRotTran = interpolate(reference, toInterpolate );
            reference = rotoTranslatePointCloud(reference, iterRotTran);
            error = computeError(reference,out);
            std::cout<<"Pointcloud error:"<<error<<std::endl;
            estRotTran = iterRotTran * estRotTran;
            iterationCount = iter;
        }
        std::cout<<"Finished Tr-ICP!"<<std::endl;
        return out;
    };
}
#endif //THIRD_ASSIGNMENT_PCMANAGER_H
