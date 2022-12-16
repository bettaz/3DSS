//
// Created by Alessio Bettarello on 24/11/22.
//
#include <iostream>
#include "opencv2/imgproc.hpp"
#include "filters.h"
#include "Evaluator.h"

int main(int argc, char **argv){
    std::vector<std::pair<std::string,std::string>> files;
    std::string pathToData = "./data/img/";
    std::string pathToGt = "./data/gt/";
    std::string names[12] = {
            "aloe.png",
            "art.png",
            "baby.png",
            "books.png",
            "bowling.png",
            "cloth.png",
            "dolls.png",
            "flowerpots.png",
            "lampshade.png",
            "laundry.png",
            "moebius.png",
            "reindeer.png"
    };
    for ( const std::string &name : names) {
        files.emplace_back(pathToData+name,pathToGt+name);
    }
    Evaluator evaluator(B,files,7,0.5,1.5,0.25,4,20,2);
    //Cyclic algorithm
    Algorithm algos[4] = {
            B,
            BM,
            JB,
            JBU,
            //IterativeUpSampling
    };
    for(Algorithm algorithm : algos){
        evaluator.SetAlgo(algorithm);
        evaluator.Evaluate();
    }
    evaluator.Save("data.csv");
    // Live GUI for spatialSigma tuning on upsampling algorithms with only one image
    //evaluator.startGUI();
}