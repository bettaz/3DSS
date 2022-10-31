//
// Created by Alessio Bettarello on 21/10/22.
//
#include "StereoMatcher.h"

std::ostream& operator<<(std::ostream& stream, matchMode e) {
    switch(e) {
        case DP: stream << "Dynamic Programming"; break;
        case Naive: stream<< "Naive"; break;
    }
    return stream;
}

StereoMatcher::StereoMatcher(const std::string &leftFileName, const std::string &rightFileName, const std::string &gtFileName, const int windowSize, const double dmin, const double focal, const double baseline,const float & lambda, const float &resizeFactor) {
    this->lambda = lambda;
    this->baseLine = baseline;
    this->dMin = dmin/resizeFactor;
    this->focalLength = focal;
    this->windowSize = windowSize;
    this->scalingFactor = resizeFactor;
    this->history = std::list<historyEntry>();
    this->fileName= leftFileName;
    this->originalLeft = cv::imread(leftFileName,cv::IMREAD_COLOR);
    this->left = cv::imread(leftFileName, cv::IMREAD_GRAYSCALE);
    this->right = cv::imread(rightFileName, cv::IMREAD_GRAYSCALE);
    this->calculateGt(gtFileName);
    this->scaleImages();
    this->height = this->left.size().height;
    this->width = this->left.size().width;
}

void StereoMatcher::getNaiveDisparities() {
    this->matchingDuration.start= std::chrono::high_resolution_clock::now();
    int half_window_size = this->windowSize / 2;
    int progress = 0;
    this->disparities = cv::Mat::zeros(this->height,this->width,CV_8UC1);
#pragma omp parallel for
    for (int i = half_window_size; i < this->height - half_window_size; ++i) {
#pragma omp parallel for
        for (int j = half_window_size; j < this->width - half_window_size; ++j) {
            int min_ssd = INT_MAX;
            int disparity = 0;
            for (int d = -j + half_window_size; d < this->width - j - half_window_size; ++d) {
                uint ssd = 0;
                for (int s = 0; s < this->windowSize * this->windowSize; s++) {
                    int tmp = (int) this->left.at<uchar>(i - half_window_size + (s / this->windowSize),
                                                         j - half_window_size + (s % this->windowSize)) -
                              (int) this->right.at<uchar>(i - half_window_size + (s / this->windowSize),
                                                          j - half_window_size + (s % this->windowSize) + d);
                    ssd += tmp * tmp;
                }
                if (ssd < min_ssd) {
                    min_ssd = ssd;
                    disparity = d;
                }
            }

            this->disparities.at<uchar>(i - half_window_size, j - half_window_size) = std::abs(disparity);
        }
        progress++;
        std::cout
                << "Computing disparities with the naive approach... "
                << std::ceil(((progress - half_window_size + 1) /
                              static_cast<double>(this->height - this->windowSize - half_window_size + 1)) * 100) << "%\r"
                << std::flush;
    }
#pragma omp critical
    {
        std::cout << "Computing disparities with the naive approach... Done.\r" << std::flush;
        std::cout << std::endl;
        this->matchingDuration.end= std::chrono::high_resolution_clock::now();
    };
}

float StereoMatcher::costCalculation(const int i, const int j, const cv::Mat &dissimilarity, cv::Mat &C,cv::Mat &M, const float &lambda) {
    if(C.at<float>(i,j) != 0.0)
        return C.at<float>(i,j);
    float diagonal = this->costCalculation(i-1,j-1, dissimilarity,C,M, lambda) + dissimilarity.at<float>(i,j);
    float lx = this->costCalculation(i-1, j, dissimilarity,C,M, lambda) + lambda;
    float ux = this->costCalculation(i,j-1,dissimilarity,C,M, lambda) + lambda;
    if (diagonal< lx){
        if(diagonal < ux){
            C.at<float>(i,j) = diagonal;
            M.at<uchar>(i,j) = match;
            return diagonal;
        }
        else{
            C.at<float>(i,j) = ux;
            M.at<uchar>(i,j) = lOcc;
            return ux;
        }
    }
    else{
        if(ux >lx){
            C.at<float>(i,j) = lx;
            M.at<uchar>(i,j) = rOcc;
            return lx;
        }
        else{
            C.at<float>(i,j) = ux;
            M.at<uchar>(i,j) = lOcc;
            return ux;
        }
    }
}

PointCloud StereoMatcher::reconstruct3DPointCloud(const std::string &filename) {
    this->reconstructionDuration.start = std::chrono::high_resolution_clock::now();
    PointCloud pc = PointCloud(this->width);
    for (int i = 0; i < this->height - this->windowSize; ++i)
    {
        std::cout << "Reconstructing 3D point cloud from disparities... " << std::ceil(((i) / static_cast<double>(this->height - this->windowSize + 1)) * 100) << "%\r" << std::flush;
        for (int j = 0; j < this->width - this->windowSize; ++j)
        {
            if (this->disparities.at<uchar>(i, j) <= 0){
                pc.pushPoint(0.0, 0.0, 0.0,0,0,0);
                continue;
            }
            // Fetches the disparities and adds the dmin displacement due to cropping
            const double d = static_cast<double>(this->disparities.at<uchar>(i, j)) + this->dMin;
            // u1 is the current position on the source image while u2 becomes the diff with left disparities (d=u1-u2)
            const double u1 = j;
            const double u2 = u1 - d;
            // adjust the Y-axis reference
            const double v1 = this->height - i;

            // Z = fb/d, X=-b(u1+u2)/2d, Y=bv/d
            const double Z = this->focalLength * this->baseLine / d;
            const double X = 0.0 - (this->baseLine * (u1 + u2) / (2 * d));
            const double Y = this->baseLine * v1 / d;
            cv::Vec3b colors = this->originalLeft.at<cv::Vec3b>(i,j);
            pc.pushPoint(X, Y, Z, colors[2],colors[1],colors[0]);
        }
    }
    this->reconstructionDuration.end = std::chrono::high_resolution_clock::now();
    this->normalsDuration.start = std::chrono::high_resolution_clock::now();
    pc.calculateNormals();
    this->normalsDuration.end = std::chrono::high_resolution_clock::now();
    pc.savePLY(filename);
    std::cout << "Reconstructing 3D point cloud from disparities... Done.\r" << std::flush;
    std::cout << std::endl;
    return pc;
}

PointCloud StereoMatcher::stereoMatch(const std::string outputFileName, const matchMode &mode) {
    this->matchingDuration.start= std::chrono::high_resolution_clock::now();
    int half_window_size = this->windowSize / 2;
    this->disparities = cv::Mat::zeros(this->height,this->width,CV_8UC1);
    if (mode == DP){
        int progress = 0;
        //for each row
#pragma omp parallel for
        for (int y_0 = half_window_size; y_0 < this->height - half_window_size; ++y_0) {
            // dissimilarity(i,j) for each (i,j)
            cv::Mat dissimilarity = this->getDissimilarityMatrix(half_window_size, y_0);
            // allocate C and M
            cv::Mat C = cv::Mat::zeros(this->width, this->width, CV_32F);
            cv::Mat M = cv::Mat::zeros(this->width, this->width, CV_8UC1);
            // initialise C and M
            // (these nodes have no predecessor)
            // populate C and M (Dynamic programming approach, recursive function evaluation)
            // for(horizontally)
            //   for(vertically)
            // a lambda proportioned starting vector for C and their occlusion direction.
            for (int i = 1; i < this->width; ++i) {
                C.at<float>(0,i) = i*lambda;
                C.at<float>(i,0) = i*lambda;
                M.at<uchar>(0,i) = lOcc;
                M.at<uchar>(i,0) = rOcc;
            }
            C.at<float>(0,0) = 1;
            M.at<directions>(0,0) = match;

            // Calculates costs
            float cost = this->costCalculation(this->width-1,this->width-1,dissimilarity,C,M,lambda);
#pragma omp critical
            {
                // trace sink->source (from bottom-right to top-left of C/M)
                // fill i-th row of disparities
                progress++;
                std::cout
                        << "Calculating disparities with Dynamic Programming... "
                        << std::ceil(((progress - half_window_size + 1) /
                                      static_cast<double>(this->height - this->windowSize - half_window_size + 1)) * 100)
                        << "%\r"
                        << std::flush;

                int i=this->width-1;
                int j=this->width-1;
                while(i >= 0 && j >= 0) {
                    switch (M.at<uchar>(i,j)) {
                        case match:
                            this->disparities.at<uchar>(y_0,i) = std::abs(j-i);
                            i--;
                            j--;
                            break;
                        case lOcc:
                            this->disparities.at<uchar>(y_0,i) = 0;
                            j--;
                            break;
                        case rOcc:
                            i--;
                            break;
                    }
                }
            };
        }
    }
    else{
        this->getNaiveDisparities();
    }
    cv::normalize(this->disparities,this->disparities, 0, 255, cv::NORM_MINMAX);
    std::string formatFilename = outputFileName.substr(0,outputFileName.length()-4)+"_ws"+std::to_string(windowSize)+"_s"+std::to_string(this->scalingFactor)+"_l"+std::to_string(this->lambda)+".png";
    cv::imwrite(formatFilename,this->disparities);
    this->matchingDuration.end= std::chrono::high_resolution_clock::now();
    PointCloud pc = this->reconstruct3DPointCloud(outputFileName);
    this->computeMetrics();
    this->saveMetrics();
    return pc;
}

StereoMatcher::~StereoMatcher() {
    left.deallocate();
    right.deallocate();
}

cv::Mat StereoMatcher::getDissimilarityMatrix(const float &halfWindowSize, const int &row) {
    cv::Mat dissimilarity = cv::Mat::zeros(this->width, this->width,CV_32FC1);
    int iMax = this->width - halfWindowSize;
#pragma omp parallel for
    for (int i = halfWindowSize; i < iMax; ++i) {
#pragma omp parallel for
        for (int j = halfWindowSize; j < iMax; ++j) {
            float sum = 0.0;
            // Fetches and stores the abs difference of the luminance of a pixels' window in the two images.
            for (int u = -halfWindowSize; u <= halfWindowSize; ++u) {
                for (int v = -halfWindowSize; v <= halfWindowSize; ++v) {
                    float i1 = static_cast<float>(this->left.at<uchar>(row+v,i+u));
                    float i2 = static_cast<float>(this->right.at<uchar>(row+v,j+u));
                    sum+=std::abs(i1-i2);
                }
            }
            dissimilarity.at<float>(i,j) = sum;
        }
    }
    return dissimilarity;
}

void StereoMatcher::showMetrics() {
    // Prints previous run metrics
    std::cout<<"<------------------------------->"<<std::endl;
    std::cout << "The SSD is " << std::to_string((this->metrics)[0]) << std::endl;
    std::cout << "The SSIM is " << std::to_string(this->metrics[1]) << std::endl;
    std::cout << "The NCC is " << std::to_string(this->metrics[2]) << std::endl;
    std::cout<<"<------------------------------->"<<std::endl;
    return;
}

std::string StereoMatcher::durationFormatter(const struct duration &param) {
    //  A human-readable format for ms times
    std::chrono::high_resolution_clock::duration duration = param.end - param.start;
    return std::to_string(std::chrono::duration_cast<std::chrono::hours>(duration).count())
           +"h"+std::to_string(std::chrono::duration_cast<std::chrono::minutes>(duration).count() % 60)
           +"'"+std::to_string(std::chrono::duration_cast<std::chrono::seconds>(duration).count() % 60)
           +"\""+std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(duration).count() % 1000)+"ms";
}

void StereoMatcher::saveCSVPerformances(std::ofstream &filestream) {
    // Const iterator to pass and format in CSV the history
    for (std::list<historyEntry>::const_iterator it=this->history.begin();it!=this->history.end();it++) {
        filestream<<"\""<<this->fileName.substr(0, this->fileName.length()-4)<<"\";"
                  <<it->lambda<<";"
                  <<it->windowSize<<";"
                  <<it->scalingFactor<<";"
                  <<it->mode<<";"
                  <<this->durationFormatter(it->matchingDuration)<<";"
                  <<std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(it->matchingDuration.end-it->matchingDuration.start).count())<<";"
                  <<this->durationFormatter(it->reconstructionDuration)<<";"
                  <<std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(it->reconstructionDuration.end-it->reconstructionDuration.start).count())<<";"
                  <<it->metrics[0]<<";"
                  <<it->metrics[1]<<";"
                  <<it->metrics[2]<<std::endl;
    }
}

void StereoMatcher::computeMetrics() {
    cv::Mat diff = cv::Mat::zeros(this->disparities.size[0],this->disparities.size[1],CV_8UC1);
    cv::Mat first(this->gt);
    first.convertTo(first,CV_8UC1);
    cv::Mat second(this->disparities);
    // SSD calculation
    cv::absdiff(first, second, diff);
    cv::Scalar sqrsum = diff.dot(diff);
    this->metrics[0]=*sqrsum.val;
    // SSIM
    cv::Scalar ssim = cv::quality::QualitySSIM::compute(first, second, diff);
    this->metrics[1]=*ssim.val;
    // Normalised cross correlation
    cv::matchTemplate(first, second, diff, cv::TM_CCORR_NORMED);
    this->metrics[2]=*(cv::sum(diff).val);
    //cv::imshow("NCC", diff);
}

std::string StereoMatcher::CSVHeader() {
    // Standard csv heading file for matching performances
    return "\"fileName\",\"lambda\";\"window size\";\"scaling\";\"matching mode\";\"matching execution time\";\"in ms\";\"3D reconstruction execution time\";\"in ms\";\"gtSSD\";\"gtSSIM\";\"gtNCC\"\n";
}

void StereoMatcher::saveMetrics() {
    // For further in-code performance analysis.
    struct historyEntry his = {.windowSize =this->windowSize,.lambda=this->lambda,.scalingFactor=this->scalingFactor,.matchingDuration=this->matchingDuration,.gtDuration=this->gtDuration,.reconstructionDuration=this->reconstructionDuration,.normalsDuration=this->normalsDuration,.mode=this->matchingMode};
    his.metrics[0] = this->metrics[0];
    his.metrics[1] = this->metrics[1];
    his.metrics[2] = this->metrics[2];
    this->history.push_back(his);
}

void
StereoMatcher::updateParams(const int windowSize, const float &lambda) {
    this->windowSize = windowSize;
    this->lambda = lambda;
}

void StereoMatcher::scaleImages() {
    cv::resize(this->left,this->left,cv::Size(), this->scalingFactor,this->scalingFactor, cv::INTER_LINEAR);
    cv::resize(this->right,this->right,cv::Size(), this->scalingFactor,this->scalingFactor, cv::INTER_LINEAR);
    cv::resize(this->originalLeft, this->originalLeft, cv::Size(), this->scalingFactor,this->scalingFactor,cv::INTER_LINEAR);
    cv::resize(this->gt,this->gt,cv::Size(), this->scalingFactor,this->scalingFactor, cv::INTER_LINEAR);
}

void StereoMatcher::calculateGt(const std::string &fileName) {
    // If ground truth is provided, uses it O(1)
    if (fileName!=""){
        this->gt = cv::imread(fileName,cv::IMREAD_GRAYSCALE);
        this->gtDuration.start=std::chrono::high_resolution_clock::now();
        this->gtDuration.end = std::chrono::high_resolution_clock::now();
    }
    else{
        // If no ground truth is given, calculates a Stereo Estimation using openCV
        this->gtDuration.start = std::chrono::high_resolution_clock::now();
        auto gtCalculator = cv::StereoSGBM::create(0,16,3,0,0,0,0,0,0,0,cv::StereoSGBM::MODE_HH4);
        gtCalculator->compute(this->left,this->right, this->gt);
        // As alternative, you can use the naive ones (less accurate)
        //this->getNaiveDisparities(this->gt);
        cv::normalize(this->gt,this->gt,0,255,cv::NORM_MINMAX);
        this->gt.convertTo(this->gt,CV_8UC1);
        this->gtDuration.end = std::chrono::high_resolution_clock::now();
        cv::imwrite("./data/computed_gt"+std::to_string(this->scalingFactor)+".png",this->gt);
        //cv::imshow("Ground Truth", this->gt);
        //cv::waitKey(0);
    }
}
