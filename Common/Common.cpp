#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <vector>
#include <random>
#include "Common.hpp"
#include "../Fred/include/frechet.hpp"
#include "../timeSeries/Search/1dCurves.hpp"
#include <math.h>
#include <map>


/*overloaded euclidean distance functions */
//our implementation is generic if you change them lsh will still work
float dist(point a,point b)
{
    float sum=0;
    for(int i=0;i<a.d;i++)
        sum = sum + (a.coords[i]-b.coords[i])*(a.coords[i]-b.coords[i]);
    return sqrt(sum);
}
float dist(point a,float * b)
{
    float sum=0;
    for(int i=0;i<a.d;i++)
        sum = sum + (a.coords[i]-b[i])*(a.coords[i]-b[i]);
    return sqrt(sum);
}
float dist(float *a,float * b,int d)
{
    float sum=0;
    for(int i=0;i<d;i++)
        sum = sum + (a[i]-b[i])*(a[i]-b[i]);
    return sqrt(sum);
}

float * strToTS(const char *buf,int &d,std::string &id)
{
    char *str=strdup(buf);
    std::stringstream ss;
    ss.str(str);
    std::string t;
    std::vector<float >res ;
    int count=0;
    while ( ss >> t) 
    {
     if(count==0)
        id=t;
     else
        res.push_back(std::stof(t)); //fill your array one by one
     count++;
    }
    if(d==-1)
    {
        d=res.size();
        //printf("d=%d\n",d);
    }
    
    
    //printf("id=%d\n",*id);
    float *coords=new float[d];
    for(int i=0;i<res.size();i++)
        coords[i]=res[i];
    free(str);
    return coords;
}
//this function takes as input a string in the format id cord1 coord2 ...]
//and returns the coords, it also returns (through pointer) the point's id 
//finally, it makes sure all points have the same number of coordinates d
float * strToCoords(const char *buf,int &d,std::string &id)
{
    char *str=strdup(buf);
    std::stringstream ss;
    ss.str(str);
    std::string t;
    std::vector<float >res ;
    int count=0;
    while ( ss >> t) 
    {
     if(count==0)
        id=t;
     else
        res.push_back(std::stof(t)); //fill your array one by one
     count++;
    }
    if(d==-1)
    {
        d=res.size();
        //printf("d=%d\n",d);
    }
    else
    {
        if(d!=res.size())
        {
            printf("d=%d res.size()=%d\n",d,res.size());
            printf("WRONG FILE\n");
            exit(0);
        }
    }
    //printf("id=%d\n",*id);
    float *coords=new float[d];
    for(int i=0;i<res.size();i++)
        coords[i]=res[i];
    free(str);
    return coords;

}



std::vector<std::pair<point*, float>> findNNBruteForce(point p,std::vector<point *> data,int k,float (*distance_funct)(point a,point b))
{
    //printf("GOT TO HERE\n");

    point **b=new point *[k]; //k best candidates
    float *Db=new  float[k];//their distances
    for(int i=0;i<k;i++)
        Db[i]=1999999;
    for(int i=0;i<data.size();i++)
    {
        float tempDist=distance_funct(*data[i],p);        
        for(int j=0;j<k;j++)
        {
            if(Db[j] > tempDist)
            {
                //printf("fuckin hell\n");
                Db[j]=tempDist;
                b[j]=data[i];
                break;
            }
        }
    }
    std::vector<std::pair<point*, float>>result;
    for(int i=0;i<k;i++)
    {
        if(Db[i]!=1999999)
            result.push_back(std::make_pair(b[i],Db[i]));
    }   
    delete[] b;
    delete[] Db;
    //printf("GOT OUT OF HERE\n");
    return result;
}
std::vector<std::pair<discrete2DCurve*, float>> findNNBruteForce(
    discrete2DCurve p,std::vector<discrete2DCurve *> data,int k,float (*distance_funct)(discrete2DCurve a,discrete2DCurve b))
{

    discrete2DCurve **b=new discrete2DCurve *[k]; //k best candidates
    float *Db=new  float[k];//their distances
    for(int i=0;i<k;i++)
        Db[i]=1999999;
    for(int i=0;i<data.size();i++)
    {
        float tempDist=distance_funct(*data[i],p);        
        for(int j=0;j<k;j++)
        {
            if(Db[j] > tempDist || Db[j]==1999999)
            {
                //printf("fuckin hell\n");
                Db[j]=tempDist;
                b[j]=data[i];
                break;
            }
        }
    }
    std::vector<std::pair<discrete2DCurve*, float>>result;
    for(int i=0;i<k;i++)
    {
        if(Db[i]!=1999999)
            result.push_back(std::make_pair(b[i],Db[i]));
    }   
    delete[] b;
    delete[] Db;
    //printf("GOT OUT OF HERE\n");
    return result;
}


std::vector<std::pair<point*, float>> FrechetNNBruteForce(point p,std::vector<point *> data,int k)
{
    //printf("GOT TO HERE\n");

    point **b=new point *[k]; //k best candidates
    float *Db=new  float[k];//their distances
    for(int i=0;i<k;i++)
        Db[i]=1999999;
    for(int i=0;i<data.size();i++)
    {
        Curve query_c = pointToCurve(&p);
        point d = *data[i];
        Curve data_c = pointToCurve(&d);
        float tempDist =  Frechet::Continuous::distance(query_c,data_c).value;
        //float tempDist=distance_funct(*data[i],p);        
        for(int j=0;j<k;j++)
        {
            if(Db[j] > tempDist)
            {
                //printf("fuckin hell\n");
                Db[j]=tempDist;
                b[j]=data[i];
                break;
            }
        }
    }
    std::vector<std::pair<point*, float>>result;
    for(int i=0;i<k;i++)
    {
        if(Db[i]!=1999999)
            result.push_back(std::make_pair(b[i],Db[i]));
    }   
    delete[] b;
    delete[] Db;
    //printf("GOT OUT OF HERE\n");
    return result;
}