//
// Created by bettaz on 18/01/23.
//

#ifndef THIRD_ASSIGNMENT_PCMANAGER_H
#define THIRD_ASSIGNMENT_PCMANAGER_H
#include <random>
#include <array>
#include "Eigen/Dense"
#include <nanoflann.hpp>

namespace pc{
    using PointCloud = Eigen::MatrixX3d;
    using Quaternion = Eigen::Quaterniond;
    using KdTree = nanoflann::KDTreeEigenMatrixAdaptor<PointCloud>;
    using Rotation = Eigen::Transform<double,3,0>;
    using Translation = Eigen::Translation3d;


    inline PointCloud
    vecToEigen(const std::vector<std::array<double, 3>> &src){
        PointCloud pcOut(src.size(),3);
        for (size_t i = 0; i < src.size(); i++) {
            pcOut(i, 0) = src[i][0];
            pcOut(i, 1) = src[i][1];
            pcOut(i, 2) = src[i][2];
        }
        return pcOut;
    };


    inline std::vector<std::array<double, 3>>
    EigenToVec(const PointCloud &src){
        std::vector<std::array<double, 3>> pcOut;
        for (auto line : src.rowwise()) {
            pcOut.push_back({{line(0),line(1),line(2)}});
        }
        return pcOut;
    };


    inline PointCloud
    addNoise(const PointCloud &src, const double &mean, const double & stdDev) {
        PointCloud noisy;
        noisy = src;
        if(mean != 0.0 && stdDev != 0.0){
            std::default_random_engine generator;
            std::normal_distribution<double> dist(mean, stdDev);
            for(auto it : noisy.rowwise()){
                Eigen::RowVector3d noise;
                noise << dist(generator), dist(generator), dist(generator);
                it = it + noise;
            }
        }
        return noisy;
    };


    inline PointCloud
    rotatePointCloud(const PointCloud &src, const Rotation &rotation, const Translation & translation) {
        PointCloud out(src);
        for(auto line : out.rowwise())
            line = (rotation * (translation * line.transpose())).transpose();
        return out;
    };


    inline void
    interpolate(PointCloud reference, PointCloud source, Rotation & rotation, Translation & translation) {
        Eigen::Matrix3d covariance;
        Eigen::Vector3d center1, center2;

        size_t pc1Rows = reference.rows(), pc2Rows = source.rows();
        // Calculate the center of pcl1
        center1 = reference.colwise().mean();

        // Calculate the center of pcl2
        center2 = source.colwise().mean();

        // Move pcl1 to the center
        reference = rotatePointCloud(reference,Rotation(Eigen::Matrix3d::Identity()), Translation(-center1));

        // Move pcl2 to the center
        source = rotatePointCloud(source,Rotation(Eigen::Matrix3d::Identity()), Translation(-center2));


        // Calculate covariance matrix
        for (int i = 0; i < pc2Rows; i++) {
            covariance(0, 0) += reference.row(i).x() * source.row(i).x();
            covariance(0, 1) += reference.row(i).x() * source.row(i).y();
            covariance(0, 2) += reference.row(i).x() * source.row(i).z();
            covariance(1, 0) += reference.row(i).y() * source.row(i).x();
            covariance(1, 1) += reference.row(i).y() * source.row(i).y();
            covariance(1, 2) += reference.row(i).y() * source.row(i).z();
            covariance(2, 0) += reference.row(i).z() * source.row(i).x();
            covariance(2, 1) += reference.row(i).z() * source.row(i).y();
            covariance(2, 2) += reference.row(i).z() * source.row(i).z();
        }
        covariance /= pc1Rows;
        Eigen::MatrixXd w, u, vt;
        Eigen::JacobiSVD<Eigen::MatrixXd> jacobi(covariance, Eigen::ComputeThinV | Eigen::ComputeThinU);
        jacobi.compute(covariance);
        u=jacobi.matrixU();
        vt=jacobi.matrixV();

        // Calculate the rotation matrix
        Eigen::Matrix3d rot = vt * u;
        if (rot.determinant() < 0.0) {
            vt.row(2) *= -1.0;
            rot = vt.transpose() * u.transpose();
        }

        // Calculate the translation matrix
        rotation = Rotation ();//(rot);
        translation =Translation (center2 - (rotation * center1));
    }

    double computeError(const PointCloud &matrix, const PointCloud &matrix1) {
        PointCloud diff = matrix-matrix1;
        return diff.cwiseAbs2().colwise().sum().cwiseSqrt().sum();;
    }

    inline PointCloud
    ICP(const PointCloud &reference, const PointCloud &source, Rotation & estRot, Translation &estTras, int &iterationCount, const double &epsilon=0.001, const size_t &maxIteration=100) {
        std::cout<<"Started ICP"<<std::endl;
        double error = epsilon+1.0;
        PointCloud out(source);
        estRot = Rotation();
        estTras = Translation();
        for (int iter = 0; iter < maxIteration && error>epsilon; ++iter) {
            Rotation iterRotation;
            Translation iterTranslation;
            // construct a kd-tree index:
            KdTree tree(3, std::cref(out), 2000);
            PointCloud toInterpolate (reference.rows(),3);
            Eigen::Index closer;
            double distance;
            for (int i=0; i<reference.rows(); i++){
                tree.index_->knnSearch(reference.row(i).data(), 1, &closer, &distance);
                toInterpolate.row(i) = source.row(closer);
            }
            interpolate(reference, toInterpolate, estRot, estTras);
            out = rotatePointCloud(out,iterRotation, iterTranslation);
            error = computeError(reference,out);
            std::cout<<"Error:"<<error<<std::endl;
            estRot = estRot.rotation() * iterRotation;
            estTras.translation() = estTras.translation() + iterTranslation.translation();
            iterationCount = iter;
        }
        std::cout<<"Finished ICP!"<<std::endl;
        return out;
    };

    struct DistanceResult{
        double distance;
        int index;
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
            TrPointCloud out(index,3);
            for (auto entry : distances) {
                out.row(removeCount) = this->row(entry.index);
                removeCount++;
            }
            out.distances = this->distances;
            return out;
        }
    };

    inline PointCloud
    trICP(const PointCloud &reference, const PointCloud &source, Rotation & estRot, Translation &estTras, int &iterationCount, const float & xsi = 0.8, const double &epsilon=0.001, const size_t &maxIteration=100){
        std::cout<<"Started Tr-ICP"<<std::endl;
        double error = epsilon+1.0;
        TrPointCloud out(source);
        estRot = Rotation();
        estTras = Translation();
        for (int iter = 0; iter < maxIteration && error>epsilon; ++iter) {
            Rotation iterRotation;
            Translation iterTranslation;
            // construct a kd-tree index:
            KdTree tree(3, std::cref(out), 2000);
            TrPointCloud toInterpolate (reference.rows(),3);
            Eigen::Index closer;
            double distance;
            for (int i=0; i<reference.rows(); i++){
                tree.index_->knnSearch(reference.row(i).data(), 1, &closer, &distance);
                toInterpolate.row(i) = source.row(closer);
                toInterpolate.distances[i] = DistanceResult({.distance=distance, .index=i});
            }
            toInterpolate= toInterpolate.getTrimmedPointCloud(xsi);
            interpolate(reference, toInterpolate, estRot, estTras);
            out = rotatePointCloud(out,iterRotation, iterTranslation);
            error = computeError(reference,out);
            std::cout<<"Error:"<<error<<std::endl;
            estRot = estRot.rotation() * iterRotation;
            estTras.translation() = estTras.translation() + iterTranslation.translation();
            iterationCount = iter;
        }
        std::cout<<"Finished Tr-ICP!"<<std::endl;
        return PointCloud();
    };
}
#endif //THIRD_ASSIGNMENT_PCMANAGER_H
