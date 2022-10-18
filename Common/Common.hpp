#pragma once
#include <string>
#include <vector>

struct point;

struct discrete2DCurve
{
    point * vertices; //vertex is a 2d point
    int d; //the number of vertices 
    std::string id;
};
struct grid2D
{
    std::vector<float> t; //t vector(uniformly chosen)
    float delta; //delta
    int d; //size of the greed
};
// struct continious1Dcurve{
//     float *coords;
//     int d;
//     std::string id;
// };
//this is the point data structure used by all algorithms and structures
struct point{
    int d;
    float *coords;
    std::string id;
    discrete2DCurve *associatedCurve;//NULL by default, used only for time series
    float *curve; //used only in continious freschet can be NULL
    int curve_d; //used only in continious freschet 
 
};


//all 3 compute the euclidean distance 
//as you can see our implemenation is generic and you can use our code 
//by switching to a different distance e.g  L1
float dist(point a,float * b);
float dist(point a,point b);
float dist(float *a,float * b,int d);

//converts a line of id coord1 coord2 coord3.....
//to a point. It also changes the id and make sure all the points in the file
//have the same dimension
float * strToCoords(const char *buf,int &d,std::string &id);
float * strToTS(const char *buf,int &d,std::string &id);




//finds the exact knn  brute force
std::vector<std::pair<point*, float>> findNNBruteForce(point p,std::vector<point *> data,int k=2,float (*distance_funct)(point a,point b)=dist);
std::vector<std::pair<discrete2DCurve*, float>> findNNBruteForce(
    discrete2DCurve p,std::vector<discrete2DCurve *> data,int k,float (*distance_funct)(discrete2DCurve a,discrete2DCurve b));
std::vector<std::pair<point*, float>> FrechetNNBruteForce(point p,std::vector<point *> data,int k);
