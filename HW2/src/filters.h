//
// Created by Alessio Bettarello on 08/12/22.
//
#include "opencv2/opencv.hpp"

#ifndef SECOND_ASSIGNMENT_FILTERS_H
#define SECOND_ASSIGNMENT_FILTERS_H
namespace filters{
    void BilateralFilter(const cv::Mat &source, cv::Mat &output,const int &kernelSize,const float &sigmaSpace, const float &sigmaSpectral);
    void BilateralMedianFilter(const cv::Mat &source, cv::Mat &output, const int &kernelSize, const float & sigmaSpace, const float & sigmaSpectral);
    void JointBilateralFilter(const cv::Mat &source,const cv::Mat &joint, cv::Mat &output, int neighbourhoodDimension,
                              const float & sigmaSpace, const float & sigmaSpectral);
    void JointBilateralUpSampling(const cv::Mat &source, const cv::Mat &joint, cv::Mat &output, int neighbourhoodDimension, const float & sigmaSpace, const float &sigmaSpectral,const float & coefficient);
    void IterativeUpSampling(const cv::Mat &source, const cv::Mat &joint, cv::Mat &output, const int & windowSize, const float &sigmaSpace, const float & sigmaSpectral);
}
#endif //SECOND_ASSIGNMENT_FILTERS_H
