//
// Created by Alessio Bettarello on 21/10/22.
//

#ifndef FIRST_ASSIGNMENT_STEREOMATCHER_H
#define FIRST_ASSIGNMENT_STEREOMATCHER_H
#include <opencv2/opencv.hpp>
#include <opencv2/quality/qualityssim.hpp>
#include <string>
#include "PointCloud.h"
#include <fstream>
#include <chrono>

typedef enum matchMode{
    Naive = 0,
    DP=1
} matchMode;
typedef enum{
    lOcc = 2,
    rOcc = 1,
    match = 0
}directions;
struct duration {
        std::chrono::time_point<std::chrono::system_clock,std::chrono::system_clock::duration> start;
        std::chrono::time_point<std::chrono::system_clock,std::chrono::system_clock::duration> end;
};
struct historyEntry {
    int windowSize;
    float lambda;
    float scalingFactor;
    struct duration matchingDuration;
    struct duration gtDuration;
    struct duration reconstructionDuration;
    struct duration normalsDuration;
    float metrics[3];
    matchMode mode;
};

class StereoMatcher {
private:
    int windowSize;
    double dMin;
    double focalLength;
    double baseLine;
    float lambda;
    float scalingFactor;
    cv::Mat left;
    cv::Mat originalLeft;
    cv::Mat right;
    cv::Mat gt;
    cv::Mat disparities;
    int width;
    int height;
    std::string fileName;
    struct duration matchingDuration;
    struct duration gtDuration;
    struct duration reconstructionDuration;
    struct duration normalsDuration;
    matchMode matchingMode;

    float metrics[3];

    std::list<struct historyEntry> history;
    std::string durationFormatter(const struct duration& param);
    float costCalculation(const int i, const int j, const cv::Mat &dissimilarity, cv::Mat &C, cv::Mat &M, const float &lambda);
    cv::Mat getDissimilarityMatrix(const float &halfWindowSize, const int &row);
    void getNaiveDisparities();
    void computeMetrics();
    PointCloud reconstruct3DPointCloud(const std::string &filename);
    void saveMetrics();
    void scaleImages();
    void calculateGt(const std::string &fileName);
public:
    StereoMatcher(const std::string &leftFileName,const std::string &rightFileName, const std::string &gtFileName, const int windowSize, const double dmin, const double focalLenght, const double baseLine,const float &lambda, const float &resizeFactor);
    ~StereoMatcher();
    PointCloud stereoMatch(const std::string fileName, const matchMode &mode);
    void showMetrics();
    void saveCSVPerformances(std::ofstream &outfile);
    static std::string CSVHeader();
    void updateParams(const int windowSize, const float &lambda);
};

#endif //FIRST_ASSIGNMENT_STEREOMATCHER_H
