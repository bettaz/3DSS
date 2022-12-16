//
// Created by bettaz on 26/10/22.
//
#include "PointCloud.h"

std::string PointCloud::point::formatPNY() {
    std::ostringstream s;
    s << this->x_val<< " " << this->y_val << " " << this->z_val<<" "<<this->R<<" "<<this->G<<" "<<this->B<<" "<<this->x_norm<< " " << this->y_norm << " " << this->z_norm;
    return s.str();
}

PointCloud::PointCloud(const int &width) : width(width), outDim(0) {
    this->matrix = std::vector<point>();
}

void PointCloud::pushPoint(const float &x, const float &y, const float &z, const int &R,const int &G, const int &B) {
    this->matrix.push_back(point(x,y,z,R,G,B));
    outDim++;
}

void PointCloud::calculateNormals() {
    std::vector<point> filtered(this->matrix);
    const int count = this->matrix.size();
    const int width = this->width;
    int halfWindowSize= 1;
    for (int i = halfWindowSize; i < (count/width)-halfWindowSize; ++i) {
        for (int j = halfWindowSize; j < (count % width) - halfWindowSize; ++j) {
            cv::Mat x, y, z, A, b, unk;
            A = cv::Mat::zeros(3, 3, CV_64F);
            b = cv::Mat::zeros(3, 1, CV_64F);
            x = cv::Mat(3, 3, CV_64F);
            y = cv::Mat(3, 3, CV_64F);
            z = cv::Mat(3, 3, CV_64F);
            unk = cv::Mat(3,1,CV_64F);
            for (int k = -halfWindowSize; k < halfWindowSize + 1; ++k) {
                for (int l = -halfWindowSize; l < halfWindowSize + 1; ++l) {
                    x.at<double>(k + halfWindowSize, l + halfWindowSize) = this->matrix[((i + k) * width) + j +
                                                                                        l].x_val;
                    y.at<double>(k + halfWindowSize, l + halfWindowSize) = this->matrix[((i + k) * width) + j +
                                                                                        l].y_val;
                    z.at<double>(k + halfWindowSize, l + halfWindowSize) = this->matrix[((i + k) * width) + j +
                                                                                        l].z_val;
                }
            }
            A.at<double>(0, 0) = *(cv::sum(x.dot(x)).val);
            A.at<double>(0, 1) = *(cv::sum(x.dot(y)).val);
            A.at<double>(0, 2) = *(cv::sum(x).val);
            A.at<double>(1, 0) = *(cv::sum(x.dot(y)).val);
            A.at<double>(1, 1) = *(cv::sum(y.dot(y)).val);
            A.at<double>(1, 2) = *(cv::sum(y).val);
            A.at<double>(2, 0) = *(cv::sum(x).val);
            A.at<double>(2, 1) = *(cv::sum(y).val);
            A.at<double>(2, 2) = 9.0;

            b.at<double>(0, 0) = *(cv::sum(x.dot(z)).val);
            b.at<double>(1, 0) = *(cv::sum(y.dot(z)).val);
            b.at<double>(2, 0) = *(cv::sum(z).val);

            cv::solve(A, b, unk, cv::DECOMP_SVD);
            double normCoefficient = cv::norm(unk);
            filtered[(i * width) + j].x_norm = unk.at<double>(0, 0) / normCoefficient;
            filtered[(i * width) + j].y_norm = unk.at<double>(1, 0) / normCoefficient;
            filtered[(i * width) + j].z_norm = 0.0-std::abs(unk.at<double>(2, 0) / normCoefficient);
        }
    }
    this->matrix= filtered;
}

void PointCloud::saveXYZ(const std::string &fileName) {
    std::stringstream out3d;
    out3d << fileName << ".xyz";
    for (int i = 0; i < this->matrix.size(); ++i) {
        out3d << this->matrix[i].formatPNY() << std::endl;
    }
    std::ofstream outfile(fileName+".xyz");
    outfile<<out3d.rdbuf();
    outfile.close();
}

void PointCloud::savePLY(const std::string &fileName) {
    std::stringstream outPLY;
    outPLY<<"ply"<<std::endl;
    outPLY<<"format ascii 1.0"<<std::endl;
    outPLY<<"element vertex "<<this->outDim<<std::endl;
    outPLY<<"property float x"<<std::endl;
    outPLY<<"property float y"<<std::endl;
    outPLY<<"property float z"<<std::endl;
    outPLY<<"property uchar red"<<std::endl;
    outPLY<<"property uchar blue"<<std::endl;
    outPLY<<"property uchar green"<<std::endl;
    outPLY<<"property float nx"<<std::endl;
    outPLY<<"property float ny"<<std::endl;
    outPLY<<"property float nz"<<std::endl;
    outPLY<<"end_header"<<std::endl;
    for (int i = 0; i < this->matrix.size(); ++i) {
        outPLY << this->matrix[i].formatPNY() << std::endl;
    }
    std::ofstream outfile(fileName+".ply");
    outfile << outPLY.rdbuf();
    outfile.close();
}

void PointCloud::reconstructSurfaces(const std::string &outputFileName) {
    this->savePLY(outputFileName);
}
