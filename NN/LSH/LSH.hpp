#pragma once
#include <vector>
#include "../../Common/Common.hpp"
#include <string>
#include <sstream>
#include <fstream>
struct bucket_data{
    point *p;
    //hash_id is used to avoid too many computations of euclidean distances
    unsigned int hash_id;
};
struct bucket{
    std::vector<bucket_data >points;
};


struct hash_funct
{
    int d ; //the dimension of the point vectors and of vector v
    float t; //random uniformly in [0,w)
    float * v; //d dimensional vector whoose coords are N(0,1)
    int w; //provided by the user , usually in range [2,6]

};

struct  amplified_hash_fucnt
{
    int d; //dimension of the points 
    int k; //usually k is 4 but can go up to 10 #of weak hfs
    int *rs; // k coefficents for linear combination
    hash_funct *hs; //k weak hash functions
    int w; //provided by the user , usually in range [2,6]
  
};

class LSH
{
    private:
        int d;  //the dimension of the points
        int w;  //user specified in range [2,6] usually
        int k; //number of weak hash functions used by each amplified hf usually 4
        int L;  //number of amplified hash functions
        int tableSize; //usually n/4 where n=number of data points
        int n; //number of points
        bucket **tables; //all L hash tables each for an ampl hash funct
        amplified_hash_fucnt *gs;
        float (*distance_funct)(point a,point b);

    public:
        //constructior only deals with hyperparams
        LSH(int d,int w,int k,int L, float (*distance_funct)(point a,point b)=dist);
        //each point is an array of floats
        //LSH gets an array of arrays of floats <-> array of points 
        //this function constructs LSH
        //if table size=-1 default tableSize =n/8
        //else the user provides the desired table size
        void fit(int n,point **points,int d,int tableSize=-1);

        //it only checks points in the SAME bucket with the SAME ID
        //if it finds less than N points, it returns only those found
        std::vector<std::pair<point*, float>>  approxKNN(point p,int N=2);
        std::vector<std::pair<point*, float>> frechetKNN(point p,int N);
        void fitFrechet(int n,std::vector<point*> data,int d,int tableSize, float delta,float e);

        //it only checks points in the SAME bucket
        std::vector<point *> approxRangeSearch(point p,float r);
        //debug
        void printTables();

        //destructor
        ~LSH();
};



//methods for creating and using hash functions
hash_funct create_hash_funct(int d,int w);
amplified_hash_fucnt create_ampl_hash_funct(int k,int d,int w);
int get_weak_hf(hash_funct h,point p);
unsigned int hash_data_point(amplified_hash_fucnt g,point p);

//this function reads a space seperated file line by line and puts it in lsh
//it does not access private parts of lsh
//it also checks if the file exists, if there is a problem with oppening the file
//it terminates
std::vector<point *> fileToLSH(LSH &lsh,const char *fileName,int &d);
