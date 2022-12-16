//
// Created by Alessio Bettarello on 12/12/22.
//

#ifndef SECOND_ASSIGNMENT_EVALUATOR_H
#define SECOND_ASSIGNMENT_EVALUATOR_H


#include <vector>
#include <string>
#include <fstream>
#include <opencv2/core.hpp>
#include <opencv2/quality.hpp>

enum Algorithm{
    B = 0,
    BM = 1,
    JB = 2,
    JBU = 3,
    IterativeUpSampling = 4
};
struct Record{
    std::string fileName;
    Algorithm algorithm;
    int windowSize;
    float spatialSigma;
    float spectralSigma;
    float SSD;
    float SSIM;
    float NCC;
    time_t executionTime;
}typedef Record;

class Evaluator {
private:
    const std::vector<std::pair<std::string,std::string>> &filenames;
    const float startingSpatialSigma;
    const float endingSpatialSigma;
    const float spatialSigmaPeace;
    const float startingSpectralSigma;
    const float endingSpectralSigma;
    const float spectralSigmaPeace;
    const int windowSize;
    std::vector<Record> records;
    Algorithm algo;
    bool withNoise = false;
    static void printCSVHeader(std::ofstream &FILE);
    static float * computeMetrics(const cv::Mat &first,const cv::Mat &second);
    void cropWindowEdges(const cv::Mat & src, const int & ws, cv::Mat &output);
public:
    Evaluator(const Algorithm &algorithm, const std::vector<std::pair<std::string,std::string>> & fileList, const int &ws, const float & spatialSigma0, const float & spatialSigma1, const float & spatialSigmaPeace, const float & spectralSigma0, const float &spectralSigma1, const float & spectralSigmaPeace): algo(algorithm), filenames(fileList), startingSpatialSigma(spatialSigma0), endingSpatialSigma(spatialSigma1), spatialSigmaPeace(spatialSigmaPeace), startingSpectralSigma(spectralSigma0), endingSpectralSigma(spectralSigma1), spectralSigmaPeace(spectralSigmaPeace), windowSize(ws){
    };
    void Evaluate();
    void SetNoise(const bool &val);
    void SetAlgo(const Algorithm &algorithm);
    void Save(const std::string &outputFile);
    void startGUI();
};


#endif //SECOND_ASSIGNMENT_EVALUATOR_H
