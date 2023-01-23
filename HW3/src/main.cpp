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
    pc::Affine rotoTranslation;
    rotoTranslation.linear() = (zRot * yRot * xRot).toRotationMatrix();
    rotoTranslation.translation() = pc::Translation (0.0,0.0,0.0);
    for(std::string fileName :fileNames){
        std::string abs = "./"+dataFolder+"/"+fileName+ext;
        happly::PLYData gtPLYPointCloud(abs);
        pc::PointCloud gt = pc::vecToEigen(gtPLYPointCloud.getVertexPositions());
        for(double stdDev : sigmas){
            // Add Gaussian noise to GT
            pc::PointCloud noisyGt = pc::addNoise(gt,mean,stdDev);
            happly::PLYData noisyGtPLY;
            std::vector<std::array<double,3>> out = pc::EigenToVec(noisyGt);
            std::string noisyGtFileName = "./" + dataFolder + "/" + fileName + "_noisy" + ext;
            noisyGtPLY.addVertexPositions(out);
            noisyGtPLY.write(noisyGtFileName, happly::DataFormat::ASCII);
            // Rotate
            pc::PointCloud rotGt = pc::rotoTranslatePointCloud(gt, rotoTranslation);
            // Add Gaussian noise to the rotated one
            pc::PointCloud noisyRot = pc::addNoise(rotGt, mean, stdDev);
            happly::PLYData noisyRotPLY;
            out = pc::EigenToVec(noisyRot);
            std::string noisyRotFileName = "./" + dataFolder + "/" + fileName + "_noisy_rotated" + ext;
            noisyRotPLY.addVertexPositions(out);
            noisyRotPLY.write(noisyRotFileName, happly::DataFormat::ASCII);

            pc::Affine estRotTran;
            int iterations;
            pc::PointCloud icpOut = pc::ICP(noisyGt, noisyRot, estRotTran,iterations, 0.5);
            happly::PLYData ICPPLY;
            out = pc::EigenToVec(icpOut);
            std::string ICPFileName = "./" + dataFolder + "/" + fileName + "ICP" + ext;
            ICPPLY.addVertexPositions(out);
            ICPPLY.write(ICPFileName, happly::DataFormat::ASCII);

            /*pc::PointCloud trIcpOut = pc::trICP(noisyGt, noisyRot,estRot,translation,iterations,0.8);
            happly::PLYData TrICPPLY;
            out = pc::EigenToVec(trIcpOut);
            std::string TrICPFileName = "./" + dataFolder + "/" + fileName + "_TrICP" + ext;
            TrICPPLY.addVertexPositions(out);
            TrICPPLY.write(TrICPFileName, happly::DataFormat::ASCII);*/

        }
    }
    return 0;
}
