#pragma once 
#include "../../Common/Common.hpp"
#include "../LSH/LSH.hpp"
#include <unordered_map>

struct  f_funct
{
    hash_funct h; //taken from LSH (file LSH.cpp/hpp)
    /*
    if hashed_points[i]==false -> return 0
    if hashed_points[i]==True ->return 1
    else toss a coin and store the outcome
    */
    std::unordered_map<int,bool>hashed_points;

};
struct hq_bucket//aka a vertex
{
    std::vector<point *> points;
};

class HyperCube{

    private:
        int k; //d' a.k.a number of bits per vertex, max=32
        int M ;//max number of points allowed to check
        int probes; //max number of vertices allowed to check
        hq_bucket * table; //a collection of vertices , each vertice has (a lot of) points
        f_funct * coin_tossers; //total k coin_tossers
        int tableSize; //2^k, max=2^32 
    public:
        //just like LSH, the constructor only deals with hyperparameters
        //it also seeds the required generators
        HyperCube(int k,int M,int probes);
        //just like LSH, this function allocates space and store points
        void fit(int n,point **points,int d);

        //it only checks for points in the same vertex or in vertices nearby(according
        //to probes), it only checks up to M points. Hence it might return less than
        //the desired number N. Especially if k is very big compared to the number of points
        std::vector<std::pair<point*, float>>  approxKNN(
            point p,int N=2);

        //it only checks for points in the same vertex or in vertices nearby(according
        //to probes), it only checks up to M points. If you use it for clustering it is 
        //in your best intereset to have a small k and relatively large M and probes
        std::vector<point *> approxRangeSearch(point p,float r);
        //you can put your debuging code here
        void printTables();

        //destructor
        ~HyperCube();
};

//this function just like fileToLSH, reads a file and stores all the points to hypercube
//it makes sure all points are of the same dimension and returns a vector of them.
//it does not need to know how hypercube works, and it does not access its private parts
std::vector<point *> fileToHyperCube(HyperCube &hq,const char *fileName,int &d);
