#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <cstring>
#include <random>
#include <assert.h>
#include <math.h>
#include <unordered_map>
#include "../../Common/Common.hpp"
#include "../../NN/LSH/LSH.hpp"


extern std::default_random_engine generator;


//generates a grid according to slides
grid2D createGrid(float delta,int d)
{
    grid2D n;
    n.delta=delta;
    n.d=d;
    std::uniform_real_distribution<float> uni_distribution(0.0,delta);
    //first of all pick the coordinates of t
    for(int i=0;i<d;i++)
    {  
        n.t.push_back(uni_distribution(generator));
    }
    
    
    
    return n;
}
float myFloor(float a)
{
    int r=a;
    return float(r);
}
float *TS_to_Vector(discrete2DCurve curve,grid2D G,int vectorSize,float paddValue)
{
    int d=G.d;
    //snap all the coordinates of the curve onto the grid
    std::vector<point > unpadded;
    std::vector<point > unpaddedWithDuplicates;
    //for each vertex of the curve
    for(int i=0;i<curve.d;i++)
    {
        //snapped the coordinates according to the formulla 
        //discused in classroom and eclass forum
        point snappedCoord;
        snappedCoord.d=2;
        snappedCoord.coords=new float[2];
        snappedCoord.coords[0]=
            myFloor((curve.vertices[i].coords[0]-G.t[i])/G.delta + 1/2)*G.delta+G.t[i];
        snappedCoord.coords[1]=
            myFloor((curve.vertices[i].coords[1]-G.t[i])/G.delta + 1/2)*G.delta+G.t[i];

        unpaddedWithDuplicates.push_back(snappedCoord);
        

        
    }
    //go  through all the snapped coordinates and remove duplicates(in a row)
    unpadded.push_back(unpaddedWithDuplicates[0]);
    for(int i=1;i<unpaddedWithDuplicates.size();i++)
    {
        if(dist(unpaddedWithDuplicates[i-1],unpaddedWithDuplicates[i])==0)
            //indicesToRemove.push_back(i);
            continue;
        else
        {
            unpadded.push_back(unpaddedWithDuplicates[i]);
        }
    }
    
    //concatenate the coordinates to a single vector x
    std::vector<float > concatedUnpadded;
    for(int i=0;i<unpadded.size();i++)
    {
        concatedUnpadded.push_back(unpadded[i].coords[0]);
        concatedUnpadded.push_back(unpadded[i].coords[1]);

    }
    float *x=new float[vectorSize];
    for(int i=0;i<vectorSize;i++)
    {
        if(i<concatedUnpadded.size())
            x[i]=concatedUnpadded[i];
        else
            x[i]=paddValue;
    }
    for(int i=0;i<unpaddedWithDuplicates.size();i++)
        delete[] unpaddedWithDuplicates[i].coords;
    return x;
}
float max(float a,float b)
{
    if(a>b)
        return a;
    return b;
}
float min(float a,float b,float c)
{
    if(a<b && a<c)
        return a;
    if(b<a && b<c)
        return b;
    return c;
        
}

float getFrechet(discrete2DCurve A,discrete2DCurve B)
{
    float **distanceMatrix=new float *[A.d];
    for(int i=0;i<A.d;i++)
        distanceMatrix[i]=new float [B.d];
    distanceMatrix[0][0]=dist(A.vertices[0],B.vertices[0]);
    for(int i=1;i<A.d;i++)
    {
        distanceMatrix[i][0] = max(
            distanceMatrix[i-1][0], dist(A.vertices[i],B.vertices[0]));
    }
    for(int j=1;j<B.d;j++)
        distanceMatrix[0][j] = max(distanceMatrix[0][j-1], dist(A.vertices[0], B.vertices[j]));
    for(int i=1;i<A.d;i++)
    {
        for(int j=1;j<B.d;j++)
        {
            distanceMatrix[i][j] = 
            max
            (
                min(
                    distanceMatrix[i-1][j], distanceMatrix[i][j-1], 
                    distanceMatrix[i-1][j-1]
                    ),
                dist(A.vertices[i],B.vertices[j])
            );
        }

    }
    float res=distanceMatrix[A.d-1][B.d-1];
    for(int i=0;i<A.d;i++)
        delete[] distanceMatrix[i];
    delete[] distanceMatrix;
    return res;
}
float getFrechet(point a,point b)
{

    return (getFrechet(*(a.associatedCurve),*(b.associatedCurve)));
}
discrete2DCurve pointTo2dCurve(point *a)
{
    discrete2DCurve cv;
    cv.d=a->d;
    cv.id=a->id;
    cv.vertices=new point[a->d];
    for(int i=0;i<a->d;i++)
    {
        cv.vertices[i].coords=new float[2];
        cv.vertices[i].coords[0]=i+1;
        cv.vertices[i].coords[1]=a->coords[i];
        cv.vertices[i].d=2;
    }
    return cv;
}