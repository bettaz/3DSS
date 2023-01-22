#include "happly/happly.h"
#include "Eigen/Geometry"
#include "pcManager.h"
using Eigen::MatrixX3d;
using Eigen::AngleAxisd;
using Eigen::Quaterniond;
using Eigen::Matrix3d;

int main(int argc, char **argv){
    std::vector<std::string> fileNames = {"airplane","ant","beethoven"};
    std::vector<double> sigmas = {0.0,0.5,1.0};
    std::string ext = ".ply";
    std::string dataFolder = "data";
    const double mean = 0.0;
    double angleX = 10.0*M_PI/180.0;
    AngleAxisd xRot = AngleAxisd(angleX,Eigen::Vector3d ::UnitX());
    double angleY = 2.0*M_PI/180.0;
    AngleAxisd yRot = AngleAxisd(angleY,Eigen::Vector3d ::UnitY());
    double angleZ = 3.0*M_PI/180.0;
    AngleAxisd zRot = AngleAxisd(angleZ,Eigen::Vector3d ::UnitZ());
    pc::Rotation totalRotation;
    totalRotation = zRot * yRot * xRot;
    pc::Translation translation(1.0,2.0,5.0);
    for(std::string fileName :fileNames){
        std::string abs = "./"+dataFolder+"/"+fileName+ext;
        happly::PLYData gtPLYPointCloud(abs);
        pc::PointCloud gt = pc::vecToEigen(gtPLYPointCloud.getVertexPositions());
        for(double stdDev : sigmas){
            // Add Gaussian noise to GT
            pc::PointCloud noisyGt = pc::addNoise(gt,mean,stdDev);
            // Rotate
            pc::PointCloud rotGt = pc::rotatePointCloud(gt, totalRotation, translation);
            // Add Gaussian noise to the rotated one
            pc::PointCloud noisyRot = pc::addNoise(rotGt, mean, stdDev);
            pc::Rotation estRot;
            int iterations;
            pc::PointCloud icpOut = pc::ICP(noisyGt, noisyRot, estRot,translation,iterations);
            pc::PointCloud trIcpOut = pc::trICP(noisyGt, noisyRot,estRot,translation,iterations,0.8);
        }
    }
    return 0;
}
