#include <chrono>
#include "happly/happly.h"
#include "Eigen/Geometry"
#include "pcManager.h"

struct ICPResult{
    std::string filename;
    bool isTrimmed;
    int iterations;
    size_t PCsize;
    double rotationX;
    double rotationY;
    double rotationZ;
    double rotEstErrX;
    double rotEstErrY;
    double rotEstErrZ;
    double transX;
    double transY;
    double transZ;
    double tranEstErrX;
    double tranEstErrY;
    double tranEstErrZ;
    double rotPerc;
    double transPerc;
    double noiseVariance;
    double totalError;
    size_t executionTime;
}Result;
void writeCSVHeader(std::ofstream &FILE) {
    const std::string out = "\"fileName\",\"isTrimmed\";\"iterations\";\"PCSize\";\"rotationX\";\"rotationY\";\"rotationZ\";\"rotEstErrX\";\"rotEstErrY\";\"rotEstErrZ\";\"transX\";\" transY\";\"transZ\";\"tranEstErrX\";\"tranEstErrY\";\"tranEstErrZ\";\"rotationPercentage\";\"translationPercentage\";\"noiseVariance\";\"totalError\";\"executionTime\"";
    FILE<<out<<std::endl;
}
void writeRecord(std::ofstream &FILE, const ICPResult & result){
    FILE<<"\""<<result.filename<<"\";"
            <<"\""<<(result.isTrimmed?"Tr-ICP":"ICP")<<"\";"
            <<result.iterations<<";"
            <<result.PCsize<<";"
            <<std::to_string(result.rotationX)<<";"
            <<std::to_string(result.rotationY)<<";"
            <<std::to_string(result.rotationZ)<<";"
            <<std::to_string(result.rotEstErrX)<<";"
            <<std::to_string(result.rotEstErrY)<<";"
            <<std::to_string(result.rotEstErrZ)<<";"
            <<std::to_string(result.transX)<<";"
            <<std::to_string(result.transY)<<";"
            <<std::to_string(result.transZ)<<";"
            <<std::to_string(result.tranEstErrX)<<";"
            <<std::to_string(result.tranEstErrY)<<";"
            <<std::to_string(result.tranEstErrZ)<<";"
            <<std::to_string(result.rotPerc)<<";"
            <<std::to_string(result.transPerc)<<";"
            <<std::to_string(result.noiseVariance)<<";"
            <<std::to_string(result.totalError)<<";"
            <<std::to_string(result.executionTime)<<std::endl;
}
std::array<double,6> calculateRotTransErr(const pc::Affine & rotation, const pc::Affine estimation){
    std::array<double,6> out;
    Eigen::Quaterniond rot(rotation.rotation());
    Eigen::AngleAxisd est(estimation.rotation());
    out.at(0) = (rot.toRotationMatrix().eulerAngles(0,1,2).coeff(0)-est.toRotationMatrix().eulerAngles(0,1,2).coeff(0));
    out.at(1) = (rot.toRotationMatrix().eulerAngles(0,1,2).coeff(1)-est.toRotationMatrix().eulerAngles(0,1,2).coeff(1));
    out.at(2) = (rot.toRotationMatrix().eulerAngles(0,1,2).coeff(2)-est.toRotationMatrix().eulerAngles(0,1,2).coeff(2));
    out.at(3) = abs(rotation.translation().coeff(0)-estimation.translation().coeff(0));
    out.at(4) = abs(rotation.translation().coeff(1)-estimation.translation().coeff(2));
    out.at(5) = abs(rotation.translation().coeff(1)-estimation.translation().coeff(2));
    return out;
}
int main(int argc, char **argv){
    std::vector<std::string> fileNames = {"airplane","ant","beethoven"};
    std::vector<double> noiseInfluences = {0.0,0.02,0.07};
    std::vector<double> rotations = {5.0,20.0,50.0};
    std::vector<double> translations = {0.0,0.01,0.02};
    std::string ext = ".ply";
    std::string dataFolder = "data";
    std::ofstream file("./data/results.csv");
    writeCSVHeader(file);
    ICPResult record;
    srand( (unsigned)time( NULL ) );
    //For each file
    for(std::string fileName :fileNames){
        std::string abs = "./"+dataFolder+"/"+fileName+ext;
        happly::PLYData gtPLYPointCloud(abs);
        pc::PointCloud gt = pc::vecToEigen(gtPLYPointCloud.getVertexPositions());
        record.filename = fileName;
        record.PCsize = gt.cols();
        pc::Point spaceRange = gt.rowwise().maxCoeff()-gt.rowwise().minCoeff();
        //For each rotational max
        for (auto rotLimit: rotations) {
            record.rotPerc = rotLimit;
            pc::Affine rotoTranslation;
            double angleX = static_cast<double>(rand())/RAND_MAX*rotLimit*M_PI/180.0;
            record.rotationX = angleX;
            Eigen::AngleAxisd xRot(angleX,Eigen::Vector3d ::UnitX());
            double angleY = static_cast<double>(rand())/RAND_MAX
                            *rotLimit*M_PI/180.0;
            record.rotationY = angleY;
            Eigen::AngleAxisd yRot(angleY,Eigen::Vector3d ::UnitY());
            double angleZ = static_cast<double>(rand())/RAND_MAX
                            *rotLimit*M_PI/180.0;
            record.rotationZ = angleZ;
            Eigen::AngleAxisd zRot(angleZ,Eigen::Vector3d ::UnitZ());
            rotoTranslation.linear() = (zRot * yRot * xRot).toRotationMatrix();
            // For each translational maximum
            for (auto transLimit: translations) {
                record.transPerc = transLimit;
                rotoTranslation.translation() = pc::Translation (
                        (static_cast<double>(rand())/RAND_MAX
                        )*transLimit*spaceRange(0),
                        (static_cast<double>(rand())/RAND_MAX
                        )*transLimit*spaceRange(1),
                        (static_cast<double>(rand())/RAND_MAX
                        )*transLimit)*spaceRange(2);
                record.transX = rotoTranslation.translation().coeff(0);
                record.transY = rotoTranslation.translation().coeff(1);
                record.transZ = rotoTranslation.translation().coeff(2);
                //For each noise variance
                for(double influence : noiseInfluences){
                    // Add Gaussian noise to GT proportionally to the spacial size of the pointcloud
                    double stdDev = spaceRange.mean()*influence;
                    record.noiseVariance = stdDev;
                    pc::PointCloud gt2(gt);
                    happly::PLYData GtPLY;
                    pc::happlyPC out = pc::EigenToVec(gt);
                    std::string noisyGtFileName = "./" + dataFolder + "/" + fileName +"R"+ std::to_string((int)rotLimit)+"T"+ std::to_string((int)(transLimit*100.0))+"N"+std::to_string((int)(influence*100.0))+ "_noisy" + ext;
                    GtPLY.addVertexPositions(out);
                    GtPLY.write(noisyGtFileName, happly::DataFormat::ASCII);
                    // Rotate
                    pc::PointCloud rotGt = pc::rotoTranslatePointCloud(gt, rotoTranslation);
                    // Add Gaussian noise to the rotated one
                    pc::PointCloud noisyRot = pc::addNoise(rotGt, 0.0, stdDev);
                    happly::PLYData noisyRotPLY;
                    out = pc::EigenToVec(noisyRot);
                    std::string noisyRotFileName = "./" + dataFolder + "/" + fileName +"R"+ std::to_string((int)rotLimit)+"T"+ std::to_string((int)(transLimit*100.0))+"N"+std::to_string((int)(influence*100.0))+ "_noisy_rotated" + ext;
                    noisyRotPLY.addVertexPositions(out);
                    noisyRotPLY.write(noisyRotFileName, happly::DataFormat::ASCII);

                    pc::Affine estRotTran;
                    int iterations;
                    double error;
                    auto start_ICP = std::chrono::high_resolution_clock::now();
                    pc::PointCloud icpOut = pc::ICP(gt, noisyRot, estRotTran,iterations, error, 1e-5);
                    auto end_ICP = std::chrono::high_resolution_clock::now();
                    record.isTrimmed = false;
                    record.iterations = iterations;
                    record.totalError = error;
                    std::array<double,6> res = calculateRotTransErr(rotoTranslation, estRotTran);
                    record.rotEstErrX = res[0];
                    record.rotEstErrY = res[1];
                    record.rotEstErrZ = res[2];
                    record.tranEstErrX = res[3];
                    record.tranEstErrY = res[4];
                    record.tranEstErrZ = res[5];
                    record.executionTime = std::chrono::duration_cast<std::chrono::milliseconds>(end_ICP-start_ICP).count();
                    writeRecord(file,record);


                    happly::PLYData ICPPLY;
                    out = pc::EigenToVec(icpOut);
                    std::string ICPFileName = "./" + dataFolder + "/" + fileName +"R"+ std::to_string((int)rotLimit)+"T"+ std::to_string((int)(transLimit*100.0))+"N"+std::to_string((int)(influence*100.0))+ "ICP" + ext;
                    ICPPLY.addVertexPositions(out);
                    ICPPLY.write(ICPFileName, happly::DataFormat::ASCII);

                    auto start_Tr = std::chrono::high_resolution_clock::now();
                    pc::PointCloud trIcpOut = pc::trICP(gt2, noisyRot, estRotTran, iterations, error, 1e-5);
                    auto end_Tr = std::chrono::high_resolution_clock::now();
                    record.isTrimmed= true;
                    record.iterations = iterations;
                    record.totalError = error;
                    res = calculateRotTransErr(rotoTranslation, estRotTran);
                    record.rotEstErrX = res[0];
                    record.rotEstErrY = res[1];
                    record.rotEstErrZ = res[2];
                    record.tranEstErrX = res[3];
                    record.tranEstErrY = res[4];
                    record.tranEstErrZ = res[5];
                    record.executionTime = std::chrono::duration_cast<std::chrono::milliseconds>(end_Tr-start_Tr).count();
                    writeRecord(file,record);
                    happly::PLYData TrICPPLY;
                    out = pc::EigenToVec(trIcpOut);
                    std::string TrICPFileName = "./" + dataFolder + "/" + fileName +"R"+ std::to_string((int)rotLimit)+"T"+ std::to_string((int)(transLimit*100.0))+"N"+std::to_string((int)(influence*100.0))+ "_TrICP" + ext;
                    TrICPPLY.addVertexPositions(out);
                    TrICPPLY.write(TrICPFileName, happly::DataFormat::ASCII);
                    file.flush();
                }
            }

        }
    }
    file.close();
    return 0;
}
