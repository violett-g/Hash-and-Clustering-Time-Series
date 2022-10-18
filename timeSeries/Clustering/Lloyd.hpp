#pragma once
#include <vector>
#include "../../NN/LSH/LSH.hpp"
#include "../../NN/Proj/Proj.hpp"
#include "../../Common/Common.hpp"
struct cluster{
    float *centroidCoords;
    std::vector<point *> points;
};
struct clusterOfCurves{
    discrete2DCurve *meanCurve;
    std::vector<discrete2DCurve *> curves;
};
//exact kmeans 
std::vector<cluster> KMeans_Exact(std::vector<point *> points,int k);
//two operator overloaded function for performing kmeans 
//using aproximate range search in the assignment step
std::vector<cluster> KMeans_RangeSearch(LSH &lsh,std::vector<point *> points,int k);
std::vector<cluster> KMeans_RangeSearch(HyperCube &hq,std::vector<point *> points,int k);

std::vector<float> getSilhouettes(std::vector<cluster> clusters,std::vector<point *>points);


std::vector<clusterOfCurves> KMeans_Exact(std::vector<discrete2DCurve *> curves,int k);
std::vector<float>  getSilhouettes(std::vector<clusterOfCurves> clusters,std::vector<discrete2DCurve *>curves);
std::vector<clusterOfCurves> KMeans_RangeSearch(
    LSH **lshs,grid2D *grids,std::vector<discrete2DCurve *> curves,
    int k,int vectorSize,int L,float padValue);