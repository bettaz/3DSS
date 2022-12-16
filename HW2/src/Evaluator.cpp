//
// Created by Alessio Bettarello on 12/12/22.
//
#include <opencv2/ximgproc/edge_filter.hpp>
#include <execinfo.h>
#include <csignal>
#include "Evaluator.h"
#include "filters.h"
std::ostream& operator<<(std::ostream& stream, Algorithm a) {
    switch (a) {
        case B:
            stream << "B";
            break;
        case BM:
            stream << "BM";
            break;
        case JB:
            stream << "JB";
            break;
        case JBU:
            stream << "JBU";
            break;
        case IterativeUpSampling:
            stream << "IterativeUpSampling";
            break;
    }
    return stream;
}

void Evaluator::Evaluate() {
    int completed = 0;
    int count=0;
    cv::Mat source, depth, out, reference, srcCmp, depthCmp, croppedOut, croppedDepth;
    auto start = std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();
    float * metrics;
    double coefficient = 0.5;
    int count1 = filenames.size();
    int count2 = ((endingSpectralSigma-startingSpectralSigma)/spectralSigmaPeace)+1;
    int count3 = ((endingSpatialSigma-startingSpatialSigma)/spatialSigmaPeace)+1;
    int total = count1*count2*count3;
    for (const std::pair<std::string,std::string> fileName : filenames) {
        for (float spectralSigma = this->startingSpectralSigma; spectralSigma <= this->endingSpectralSigma; spectralSigma= spectralSigma + this->spectralSigmaPeace) {
            for(float spatialSigma = this->startingSpatialSigma; spatialSigma <= this->endingSpatialSigma; spatialSigma= spatialSigma + this->spatialSigmaPeace){
                count++;
                completed = std::ceil(count * 100.0 / total);
                std::cout <<"specS " << spectralSigma << " spatS " << spatialSigma << " " << fileName.first << " " << fileName.second<< " Percentage of algo "<< this->algo<<" : "<<completed<< "%\r"<<std::flush;
                source = cv::imread(fileName.first,cv::IMREAD_GRAYSCALE);
                if(this->withNoise){
                    cv::Mat noise = cv::Mat::zeros(source.size(),source.type());
                    cv::randn(noise, 0, 25);
                    source += noise;
                }
                depth = cv::imread(fileName.second,cv::IMREAD_GRAYSCALE);
                switch (this->algo) {
                    case B:
                        start = std::chrono::high_resolution_clock::now();
                        filters::BilateralFilter(source, out, windowSize, spatialSigma, spectralSigma);
                        end = std::chrono::high_resolution_clock::now();
                        cv::bilateralFilter(source, reference, windowSize,spectralSigma, spatialSigma);
                        metrics = this->computeMetrics(out,reference);
                        cv::imwrite(fileName.first + "_specS_" + std::to_string(spectralSigma) + "_spatS_" + std::to_string(spatialSigma) + "_B.png", out);
                        cv::imwrite(fileName.first + "_specS_" + std::to_string(spectralSigma) + "_spatS_" + std::to_string(spatialSigma) + "_B_diff.png", (out - reference));
                        records.push_back(Record{.fileName=fileName.first,.algorithm=B,.windowSize=windowSize,.spatialSigma=spatialSigma,.spectralSigma=spectralSigma,.SSD=metrics[0],.SSIM=metrics[1],.NCC=metrics[2],.executionTime=(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count())});
                        break;
                    case BM:
                        start = std::chrono::high_resolution_clock::now();
                        filters::BilateralFilter(source, out,windowSize, spatialSigma,spectralSigma);
                        end = std::chrono::high_resolution_clock::now();
                        cv::bilateralFilter(source, reference, windowSize,spectralSigma,  spatialSigma);
                        metrics = this->computeMetrics(out,reference);
                        cv::imwrite(fileName.first + "_specS_" + std::to_string(spectralSigma) + "_spatS_" + std::to_string(spatialSigma) + "_BM.png", out);
                        cv::imwrite(fileName.first + "_specS_" + std::to_string(spectralSigma) + "_spatS_" + std::to_string(spatialSigma) + "_BM_diff.png", (out - reference));
                        records.push_back(Record{.fileName=fileName.first,.algorithm=BM,.windowSize=windowSize,.spatialSigma=spatialSigma,.spectralSigma=spectralSigma,.SSD=metrics[0],.SSIM=metrics[1],.NCC=metrics[2],.executionTime=(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count())});
                        break;
                    case JB:
                        start = std::chrono::high_resolution_clock::now();
                        filters::JointBilateralFilter(source, depth, out, windowSize, spatialSigma,spectralSigma);
                        end = std::chrono::high_resolution_clock::now();
                        cv::ximgproc::jointBilateralFilter(source, depth, reference, spectralSigma, 0.0, spatialSigma);
                        metrics = this->computeMetrics(out,reference);
                        cv::imwrite(fileName.first + "_specS_" + std::to_string(spectralSigma) + "_spatS_" + std::to_string(spatialSigma) + "_JB.png", out);
                        cv::imwrite(fileName.first + "_specS_" + std::to_string(spectralSigma) + "_spatS_" + std::to_string(spatialSigma) + "_JB_diff.png", (out - reference));
                        records.push_back(Record{.fileName=fileName.first,.algorithm=JB,.windowSize=windowSize,.spatialSigma=spatialSigma, .spectralSigma=spectralSigma,.SSD=metrics[0],.SSIM=metrics[1],.NCC=metrics[2],.executionTime=(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count())});
                        break;
                    case JBU:
                        cv::resize(depth,depthCmp,cv::Size(depth.size[1]*coefficient, depth.size[0]*coefficient));
                        start = std::chrono::high_resolution_clock::now();
                        filters::JointBilateralUpSampling(source, depthCmp, out, windowSize,  spatialSigma,spectralSigma, coefficient);
                        end = std::chrono::high_resolution_clock::now();
                        cropWindowEdges(out, spectralSigma, croppedOut);
                        cropWindowEdges(depth, spectralSigma, croppedDepth);
                        metrics = this->computeMetrics(croppedOut,croppedDepth);
                        cv::imwrite(fileName.first + "_specS_" + std::to_string(spectralSigma) + "_spatS_" + std::to_string(spatialSigma) + "JBU.png", croppedOut);
                        cv::imwrite(fileName.first + "_specS_" + std::to_string(spectralSigma) + "_spatS_" + std::to_string(spatialSigma) + "_JBU_diff.png", (croppedOut - croppedDepth));
                        records.push_back(Record{.fileName=fileName.first,.algorithm=JBU,.windowSize=windowSize,.spatialSigma=spatialSigma,.spectralSigma=spectralSigma,.SSD=metrics[0],.SSIM=metrics[1],.NCC=metrics[2],.executionTime=(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count())});
                        break;
                    case IterativeUpSampling:
                        cv::resize(depth,depthCmp,cv::Size(depth.size[1]*0.5, depth.size[0]*0.5));
                        start = std::chrono::high_resolution_clock::now();
                        filters::IterativeUpSampling(depthCmp,source, out, windowSize,  spatialSigma, spectralSigma);
                        cv::imshow("test",out);
                        cv::waitKey(0);
                        end = std::chrono::high_resolution_clock::now();
                        metrics = this->computeMetrics(out,depth);
                        cv::imwrite(fileName.first + "_specS_" + std::to_string(spectralSigma) + "_spatS_" + std::to_string(spatialSigma) + "_IteU.png", out);
                        cv::imwrite(fileName.first + "_specS_" + std::to_string(spectralSigma) + "_spatS_" + std::to_string(spatialSigma) + "_IteU_diff.png", (out - reference));
                        records.push_back(Record{.fileName=fileName.first,.algorithm=IterativeUpSampling,.windowSize=windowSize,.spatialSigma=spatialSigma,.spectralSigma=spectralSigma,.SSD=metrics[0],.SSIM=metrics[1],.NCC=metrics[2],.executionTime=(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count())});
                        break;
                }
            }
        }
    }
}

void Evaluator::Save(const std::string &outputFile) {
    std::ofstream file(outputFile);
    this->printCSVHeader(file);
    for (const Record &record : records) {
        file << record.fileName <<
             ";" << record.algorithm <<
             ";" << record.windowSize <<
             ";" << record.spatialSigma <<
             ";" << record.spectralSigma <<
             ";" << record.SSD<<
        ";"<<record.SSIM<<
        ";"<<record.NCC<<
        ";"<<record.executionTime<<
        std::endl;
    }
    file.close();
}

void Evaluator::printCSVHeader(std::ofstream &FILE) {
    const std::string out = "\"fileName\",\"algorithm\";\"window size\";\"spatialSigma\";\"spectralSigma\";\"SSD\";\"SSIM\";\"NCC\";\"duration\"";
    FILE<<out<<std::endl;
}

float *Evaluator::computeMetrics(const cv::Mat &first,const cv::Mat &second){
    static float out[3];
    cv::Mat diff = cv::Mat::zeros(first.size[0],first.size[1],CV_8UC1);
    // SSD calculation
    cv::absdiff(first, second, diff);
    cv::Scalar sqrSum = diff.dot(diff);
    out[0]=*sqrSum.val;
    // SSIM
    cv::Scalar ssim = cv::quality::QualitySSIM::compute(first, second, diff);
    out[1]=*ssim.val;
    // Normalised cross correlation
    cv::matchTemplate(first, second, diff, cv::TM_CCORR_NORMED);
    out[2]=*(cv::sum(diff).val);
    return out;
}

void Evaluator::SetNoise(const bool &val) {
    this->withNoise = val;
}

void Evaluator::cropWindowEdges(const cv::Mat &src, const int &ws, cv::Mat &output) {
    const int halfWs = ws/2;
    const int endHeight = src.size[0]-halfWs;
    const int endWidth = src.size[1]-halfWs;
    src(cv::Range(halfWs,endHeight),cv::Range(halfWs,endWidth)).copyTo(output);
}

void Evaluator::SetAlgo(const Algorithm &algorithm) {
    this->algo = algorithm;
}
cv::Mat display;
int spatialSigmaTrack = 500;
int spectralSigmaTrack=20;
int ws;
float spatialSigmaEffective=0.5;
char wName[] = "InteractiveGUI";
char spectralSigmaName[] = "spectralSigma";
char spatialSigmaName[] = "spatialSigma(in thousands)";
cv::Mat JointBilateralUpsampling, JointBilateralUpsamplingDiff;
cv::Mat source;
cv::Mat depth;
cv::Mat depthCmp;
auto lastExecution = std::chrono::high_resolution_clock::now();

void Evaluator::startGUI(){
    source = cv::imread("./data/img/books.png",cv::IMREAD_GRAYSCALE);
    depth = cv::imread("./data/gt/books.png",cv::IMREAD_GRAYSCALE);
    display = cv::Mat::ones(500,500,CV_8UC1);
    ws= this->windowSize;
    cv::resize(depth,depthCmp,cv::Size(depth.size[1]*0.65,depth.size[0]*0.65));
    filters::JointBilateralUpSampling(source, depthCmp, JointBilateralUpsampling,windowSize, spatialSigmaEffective, spectralSigmaTrack, 0.5);
    JointBilateralUpsamplingDiff = JointBilateralUpsampling-depth;
    display = JointBilateralUpsampling;
    lastExecution = std::chrono::high_resolution_clock::now();


    cv::namedWindow(wName,cv::WINDOW_AUTOSIZE);
    cv::createTrackbar(spectralSigmaName, wName, &spectralSigmaTrack, 100);
    cv::setTrackbarMin(spectralSigmaName, wName, 1);
    cv::setTrackbarPos(spectralSigmaName,wName,20);
    cv::createTrackbar(spatialSigmaName, wName, &spatialSigmaTrack, 2000);
    cv::setTrackbarMin(spatialSigmaName, wName, 5);
    cv::setTrackbarPos(spatialSigmaName,wName,500);
    while(true){
        spatialSigmaEffective = spatialSigmaTrack / 1000.0;
        auto now = std::chrono::high_resolution_clock::now();
        auto duration =std::chrono::duration_cast<std::chrono::milliseconds>(now-lastExecution);
        if( duration.count() >= 1000){
            std::cout << spectralSigmaTrack << " " << spatialSigmaTrack << std::endl;
            lastExecution = now;
            filters::JointBilateralUpSampling(source, depthCmp, JointBilateralUpsampling,ws,  spatialSigmaEffective, spectralSigmaTrack, 0.5);
            JointBilateralUpsamplingDiff = JointBilateralUpsampling-depth;
            display = JointBilateralUpsampling;
            cv::hconcat(display,JointBilateralUpsamplingDiff,display);
            cv::resize(display,display,cv::Size(display.size[1]*0.65, display.size[0]*0.65));
            std::cout<<"render complete"<<std::endl;
        }
        cv::imshow(wName,display);
        if(cv::waitKey(1000)==27)
            break;
    }
}
