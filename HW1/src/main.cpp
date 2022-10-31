#include <iostream>
#include <regex>
#include "StereoMatcher.h"

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        std::cerr << "Usage: " << argv[0] << " IMAGE1 IMAGE2 <OUTPUT_FILE> <ws=WINDOW_SIZE> <f=FOCAL_LENGTH> <bl=BASELINE> <l=LAMBDA> <dm=DMIN> <gt=GROUND_TRUTH_FILENAME>" << std::endl;
        return 1;
    }
    const std::string firstFileName = argv[1];
    const std::string secondFileName = argv[2];
    std::string outputFileName = "./data/output.png";
    std::string gtFileName = "";
    int window_size = 3;
    int dmin = 200;
    int focal = 3740;
    int baseline = 160;
    float lambda = 80;
    float scalingFactor = 1.0;
    // Uses regex to recognise parameters, missing error sensibility.
    static std::regex const sizeMatcher( "ws=*" );
    static std::regex const focalMatcher( "f=*" );
    static std::regex const dminMatcher( "dm=*" );
    static std::regex const lambdaMatcher( "l=*" );
    static std::regex const baselineMatcher( "bl=*" );
    static std::regex const groundMatcher( "gt=*" );
    static std::regex const scalingMatcher("s=*");
    for(int i = 3; i< argc; i++){
        std::string tmp = argv[i];
        if(std::regex_search(tmp, sizeMatcher)){
            window_size = std::stoi(tmp.substr(3, tmp.size() - 3));
        }
        else{
            if(std::regex_search(tmp, focalMatcher)){
                focal = std::stoi(tmp.substr(2, tmp.size() - 2));
            }
            else{
                if(std::regex_search(tmp, lambdaMatcher)){
                    lambda = std::stof(tmp.substr(3, tmp.size() - 3));
                }
                else{
                    if(std::regex_search(tmp, dminMatcher)){
                        dmin = std::stoi(tmp.substr(3, tmp.size() - 3));
                    }
                    else{
                        if(std::regex_search(tmp, baselineMatcher)){
                            baseline = std::stoi(tmp.substr(3, tmp.size() - 3));
                        }
                        else {
                            if (std::regex_search(tmp,groundMatcher)){
                                gtFileName = tmp.substr(3,tmp.size()-3);
                            }
                            else
                            if(std::regex_search(tmp,scalingMatcher)){
                                scalingFactor = std::stof(tmp.substr(2, tmp.size() - 2));
                            }
                            else{
                                outputFileName = tmp;
                            }
                        }
                    }
                }
            }
        }
    }
    std::ofstream metricsSerialiser(outputFileName.substr(0,outputFileName.length()-4)+".csv");
    metricsSerialiser<< StereoMatcher::CSVHeader();
    // user input version
    /*StereoMatcher stereo(firstFileName,secondFileName,gtFileName,window_size,dmin,focal,baseline,lambda,scalingFactor);
    PointCloud points = stereo.stereoMatch(outputFileName);
    stereo.showMetrics();
    stereo.saveCSVPerformances(metricsSerialiser);*/
    // Performance analysis version in scaling, window size for disparity matching and lambda value for
    for(float scaling=0.4; scaling <0.41; scaling += 0.2){
        StereoMatcher stereo(firstFileName,secondFileName,gtFileName,window_size,dmin,focal,baseline,lambda,scaling);
        for(int windowSize = 3; windowSize<4;windowSize+=2){
            for(int lambda= 80; lambda< 81; lambda+=20){
                std::cout<<"scaling :"<<scaling<<" ws:"<<windowSize<<" lambda:"<<lambda<<std::endl;
                stereo.updateParams(windowSize,lambda);
                //PointCloud points = stereo.stereoMatch(outputFileName,Naive);
                PointCloud points = stereo.stereoMatch(outputFileName,DP);
                stereo.showMetrics();
            }
        }
        stereo.saveCSVPerformances(metricsSerialiser);
        metricsSerialiser.flush();
    }
    metricsSerialiser.close();
    cv::waitKey(0);
}