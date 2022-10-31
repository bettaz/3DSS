//
// Created by bettaz on 26/10/22.
//

#ifndef FIRST_ASSIGNMENT_POINTCLOUD_H
#define FIRST_ASSIGNMENT_POINTCLOUD_H
#include "vector"
#include "string"
#include "sstream"
#include "fstream"
#include <list>
#include "opencv2/core.hpp"


class PointCloud{
private:
    class point{
    public:
        float x_val,y_val,z_val;
        double x_norm,y_norm,z_norm;
        int R,G,B;
        point(const float &x, const float &y, const float &z): point(x,y,z,255,255,255) {};
        point(const float &x, const float &y, const float &z, const int &R, const int &G, const int &B): x_val(x),y_val(y),z_val(z), R(R), G(G), B(B) {}
        std::string formatPNY();
    };
    int width;
    int outDim;
    std::vector<point> matrix;
public:
    PointCloud(const int &width);
    PointCloud(PointCloud const &src): matrix(src.matrix), width(src.width), outDim(src.outDim){};
    void pushPoint(const float &x, const float &y, const float &z, const int &R, const int &G, const int &B);
    void reconstructSurfaces(const std::string &outputFileName);
    void calculateNormals();
    void saveXYZ(const std::string & outputFileName);
    void reconstructSurfaces();
    void savePLY(const std::string & fileName);
};


#endif //FIRST_ASSIGNMENT_POINTCLOUD_H
