//
// Created by Alessio Bettarello on 08/12/22.
//

#include "filters.h"

void getGaussianMask(const int & kernelSize, cv::Mat &out, double &sum, float sigma = 0.0) {
    if(sigma == 0.0)
        sigma = 0.3*((kernelSize-1)*0.5 - 1) + 0.8;
    out = cv::Mat(cv::Size(kernelSize, kernelSize), CV_8UC1);

    //const double rmax = kernelSize / 2;
    const double rmax = 2.5 * sigma;

    // most of the bell of the distribution will fit in the window
    // if too small -> to narrow bell, lose information
    // if too big   -> to wide bell, get noise

    //mask = cv::getGaussianKernel(kernelSize,-1,0); // Gaussian mask

    // MASK and WEIGHT CREATION in order to normalize filtering

    for (int r = 0; r < kernelSize; ++r) {
        for (int c = 0; c < kernelSize; ++c) {
            // Box Mask
            //mask.at<uchar>(r, c) = 1;

            // Gaussian Mask
            double r2 = (r - rmax) * (r - rmax) + (c - rmax) * (c - rmax);
            out.at<uchar>(r, c) = 255 * std::exp(-r2 / (2 * sigma * sigma));
            // so center pixel will be 255 we dont lose info
            // 1.0  -> 1
            // 0.98 -> 1
            // normalize between [0,255] since it is the maximum resolution

            sum += out.at<uchar>(r, c);
        }
    }
}

void
filters::BilateralFilter(const cv::Mat &source, cv::Mat &output,const int  &kernelSize,const float &sigmaSpace, const float & sigmaSpectral) {
    const int width = source.cols;
    const int height = source.rows;

    output = cv::Mat::zeros(height, width, CV_8UC1);

    cv::Mat gaussianKernel;
    double maskSum = 0;

    getGaussianMask(kernelSize, gaussianKernel, maskSum, sigmaSpace);

    const float sigmaRangeSq = sigmaSpectral * sigmaSpectral;

    // compute the range kernel's value
    float histogram[256];
    for (int level = 0; level < 256; ++level)
    {
        histogram[level] = std::exp(-level * level / (2 * sigmaRangeSq));
        // we did it outside because we knew the values of difference
        // from [0,255]
    }

    // must center the window
    for (int r = kernelSize / 2; r < height - kernelSize / 2; ++r) {
        for (int c = kernelSize / 2; c < width - kernelSize / 2; ++c) {

            // get center intensity
            int centralPixel = static_cast<int>(source.at<uchar>(r, c));

            int sum = 0;
            float bilateralSum = 0; // to aggregate weights in this gaussianKernel
            for (int i = -kernelSize / 2; i <= kernelSize / 2; ++i) {
                for (int j = -kernelSize / 2; j <= kernelSize / 2; ++j) {

                    int intensity = static_cast<int>(source.at<uchar>(r + i, c + j));
                    // weight based on the distance from the middle pixel
                    // compute range difference to center pixel value
                    int diff = cv::abs(centralPixel - intensity); // [255,0]
                    // get range kernel's value
                    float weight_range = histogram[diff]; // boundary kept
                    // combine the range kernel's value with weights
                    int weight_spatial = static_cast<int>(gaussianKernel.at<uchar>(i + kernelSize / 2,
                                                                         j + kernelSize / 2));
                    // distribution to be applied, will be cut by the weight_range
                    // (look for immages in the presentations)
                    // combined weight
                    float weight = weight_range * weight_spatial;
                    sum += intensity * weight; // convolution
                    bilateralSum += weight;
                }
            }
            output.at<uchar>(r, c) = sum / bilateralSum;
        }
    }
}

void
filters::JointBilateralFilter(const cv::Mat &source,const cv::Mat &joint, cv::Mat &output, int kernelSize,
                              const float &sigmaSpace, const float &sigmaSpectral) {
    const int width = source.cols;
    const int height = source.rows;
    output = cv::Mat::zeros(source.size(),CV_8UC1);

    cv::Mat gaussianKernel;
    double sum;
    getGaussianMask(kernelSize,gaussianKernel, sum, sigmaSpace);

    for (int r = kernelSize / 2; r < height - kernelSize / 2; ++r) {
        for (int c = kernelSize / 2; c < width - kernelSize / 2; ++c) {

            double sum_w = 0.0;
            double sum = 0.0;
            for (int i = -kernelSize / 2; i <= kernelSize / 2; ++i) {
                for (int j = -kernelSize / 2; j <= kernelSize / 2; ++j) {
                    double rangeSpread
                            = std::abs(source.at<uchar>(r, c) - source.at<uchar>(r + i, c + j));
                    const float sigmaSq = sigmaSpectral * sigmaSpectral;
                    const float normCoefficient = std::sqrt(2 * M_PI) * sigmaSpectral;

                    double w
                            = (1 / normCoefficient) * std::exp(-rangeSpread / (2 * sigmaSq))
                              * gaussianKernel.at<float>(i + kernelSize / 2, j + kernelSize / 2);

                    sum
                            += joint.at<uchar>(r + i, c + j) * w;
                    sum_w
                            += w;
                }
            }

            output.at<uchar>(r, c) = sum / sum_w;

        }
    }
}


void filters::BilateralMedianFilter(const cv::Mat &source, cv::Mat &output,const int &kernelSize, const float &sigmaSpace, const float &sigmaSpectral) {
    const int width = source.cols;
    const int height = source.rows;

    output = cv::Mat::zeros(height, width, CV_8UC1);

    cv::Mat gaussianKernel;
    double maskSum = 0;

    getGaussianMask(kernelSize, gaussianKernel, maskSum, sigmaSpace);

    const float sigmaRangeSq = sigmaSpectral * sigmaSpectral;

    // compute the range kernel's value
    float histogram[256];
    for (int level = 0; level < 256; ++level)
    {
        histogram[level] = std::exp(-level * level / (2 * sigmaRangeSq));
    }

    // must center the window
    for (int r = kernelSize / 2; r < height - kernelSize / 2; ++r) {
        for (int c = kernelSize / 2; c < width - kernelSize / 2; ++c) {

            // get center intensity
            int intensity_center = static_cast<int>(source.at<uchar>(r, c));

            int sum = 0;
            float sum_Bilateral_mask = 0; // to aggregate weights in this gaussianKernel

            for (int i = -kernelSize / 2; i <= kernelSize / 2; ++i) {
                for (int j = -kernelSize / 2; j <= kernelSize / 2; ++j) {

                    int intensity = static_cast<int>(source.at<uchar>(r + i, c + j));
                    // weight based on the distance from the middle pixel
                    // compute range difference to center pixel value
                    int diff = cv::abs(intensity_center - intensity); // [255,0]
                    // get range kernel's value
                    float weight_range = histogram[diff]; // boundary kept
                    // combine the range kernel's value with weights
                    int weight_spatial = static_cast<int>(gaussianKernel.at<uchar>(i + kernelSize / 2,
                                                                                   j + kernelSize / 2));
                    // distribution to be applied, will be cut by the weight_range
                    // (look for immages in the presentations)
                    // combined weight
                    float weight = weight_range * weight_spatial;
                    sum += intensity * weight; // convolution
                    sum_Bilateral_mask += weight;
                }
            }
            output.at<uchar>(r, c) = sum / sum_Bilateral_mask;
        }
    }
}

void filters::JointBilateralUpSampling(const cv::Mat &source, const cv::Mat &joint, cv::Mat &output, int kernelSize,
                              const float &sigmaSpace, const float &sigmaSpectral, const float &coefficient) {
    output = cv::Mat::ones(cv::Size(source.size().width,source.size().height),CV_8UC1);
    output *=255;
    const int width = output.cols;
    const int height = output.rows;

    cv::Mat gaussianKernel;
    double sum;
    getGaussianMask(kernelSize,gaussianKernel, sum, sigmaSpace);

    for (int r = kernelSize / 2; r < height - kernelSize / 2; ++r) {
        for (int c = kernelSize / 2; c < width - kernelSize / 2; ++c) {

            double sum_w = 0.0;
            double sum = 0.0;
            for (int i = -kernelSize / 2; i <= kernelSize / 2; ++i) {
                for (int j = -kernelSize / 2; j <= kernelSize / 2; ++j) {
                    double rangeSpread
                            = std::abs(source.at<uchar>(r, c) - source.at<uchar>((r + i), (c + j)));
                    const float sigmaSq = sigmaSpectral * sigmaSpectral;
                    const float normalization = std::sqrt(2 * M_PI) * sigmaSpectral;

                    double w
                            = (1 / normalization) * std::exp(-rangeSpread / (2 * sigmaSq))
                              * gaussianKernel.at<float>(i + kernelSize / 2, j + kernelSize / 2);

                    sum
                            += joint.at<uchar>((r + i)*coefficient, (c + j)*coefficient) * w;
                    sum_w
                            += w;
                }
            }

            output.at<uchar>(r, c) = sum / sum_w;

        }
    }
}

void filters::IterativeUpSampling(const cv::Mat &source, const cv::Mat &joint, cv::Mat &output, const int &windowSize,
                                  const float &sigmaSpace, const float & sigmaSpectral) {
    int upsampleFactor = log2(joint.rows / source.rows);
    output = source;

    for (int i = 1; i < upsampleFactor; ++i) {
        // Doubling the size of the image
        resize(output, output, cv::Size(), 2, 2);

        // Downscaling guide image
        cv::Mat downscaledGuideImage;
        cv::resize(joint, downscaledGuideImage, output.size());

        // Filtering upscaled image
        filters::JointBilateralFilter(output, downscaledGuideImage, output , windowSize, sigmaSpace , sigmaSpectral);
    }
    resize(output, output, joint.size());
    filters::JointBilateralFilter(output, joint, output , windowSize, sigmaSpace, sigmaSpectral);
}