#include "Proj.hpp"
#include <iostream>
#include <random>
#include <math.h>       /* floor */
#include <algorithm>
#include <functional>
#include <stdlib.h>
#include <cmath>
#include <bits/stdc++.h>
#include <bitset>


//this variables are kept for debugging/testing purposes
//to see if hash function is fair
int numZero=0;
int numOne=0;
//uniform generator for his from LSH.cpp
extern std::default_random_engine generator;

bool coin_toss(f_funct f,point p)
{
    int hi=get_weak_hf(f.h,p);
    //hi not found
    if(f.hashed_points.find(hi)==f.hashed_points.end())
    {
        unsigned int r=rand()%2;
        //map hi to 1
        if(r==1)
        {
            f.hashed_points[hi]=true;
            numOne+=1;
        }
        //map hi to 0
        else //r=0
        {
            f.hashed_points[hi]=false;
            numZero+=1;
        }
        return f.hashed_points[hi];

    }
    //hi is allready stored
    else{
        return f.hashed_points[hi];
    }
}
f_funct create_f_funct(int d,int w)
{
    f_funct f;
    hash_funct h=create_hash_funct(d,w);
    f.h=h;
    return f;
}
//taken from stack overflow
unsigned int bitArrayToInt32(bool arr[], int count)
{
    unsigned int ret = 0;
    int tmp;
    //we only care about the count(=k) first bits
    for (int i = 0; i < count; i++) {
        tmp = arr[i];
        ret |= tmp << (count - i - 1);
    }
    return ret;
}
unsigned int hash_data_point_hypercube(point p,f_funct *coin_tossers,int k)
{
    /*
    IMPORTANT NOTE : 
        An integer has 32 bits, we only care about the first k.
        We initiallise all the bits to 0 and then determine the first k using coin tossers
        Every function in this file that operates on bits, only modifies the first k bits
    */
    bool * res=new bool[32];
    for(int  i=0;i<32;i++)
        res[i]=0;
    for(int i=0;i<k;i++)
    {
        bool cur=coin_toss(coin_tossers[i],p);
        res[i]=cur;
    }
    //return (unsigned int) res;
    int index=bitArrayToInt32(res,k);
    delete[] res;
    return index;
    
}

//taken from https://stackoverflow.com/questions/1505675/power-of-an-integer-in-c
int myPow(int x, unsigned int p)
{
  if (p == 0) return 1;
  if (p == 1) return x;
  
  int tmp = myPow(x, p/2);
  if (p%2 == 0) return tmp * tmp;
  else return x * tmp * tmp;
}
/*
    HYPERCUBE CODE BELLOW
*/
HyperCube::HyperCube(int k,int M,int probes)
{
    this->k=k;
    this->M=M;
    this->probes=probes;
    if(k <= 32 )
        this->tableSize=myPow(2,k);
    else 
    {
        this->tableSize=myPow(2,32);
        this->k=32;
    }
    this->table=NULL;
    //srand(185);
    //generator.seed(185);

}

void HyperCube:: fit(int n,point **points,int d)
{
    //first step ) create hash functions
    this->coin_tossers=new f_funct[this->k];
    for(int i=0;i<this->k;i++)
        this->coin_tossers[i]=create_f_funct(d,300);
    //second step) create hash table
    //printf("this->tableSize=%d\n",tableSize);
    this->table=new hq_bucket[tableSize];
    //thirst step) hash all the points
    for(int i=0;i<n;i++)
    {
        unsigned int index=hash_data_point_hypercube(*(points[i]),this->coin_tossers,this->k);
        //printf("index = %d\n",index);
        this->table[index].points.push_back(points[i]);
    }
    

}
 


unsigned int toggleKthBit(unsigned int n, int k)
{
    return (n ^ (1 << (k-1)));
}

std::vector<int> magic(int bits,int i,int changesLeft)
{
    if(changesLeft==0)
    {
        std::vector<int> res;
        res.push_back(bits);
        return res;
    }
    else if(i<1)
    {
        std::vector<int> res;
        return res;
    }
    std::vector<int>res;
    //togle ith bit
    bits=toggleKthBit(bits,i);
    
    std::vector<int> res2=magic(bits, i-1, changesLeft-1);
    //untoggle it
    bits=toggleKthBit(bits,i);
    std::vector<int> res3 = magic(bits, i-1, changesLeft);
    if(res2.size()>0)
        res.insert(res.end(),res2.begin(),res2.end());
    if(res3.size()>0)
        res.insert(res.end(),res3.begin(),res3.end());
    return res;

    
}
std::vector<int > getAdjacentVertices(int index,int k,int probes,int tableSize)
{
    std::vector<int> res;
    res.push_back(index);
    probes=probes-1;
    if(probes==0)
        return res;
    int maxDistance=50; //very large
    for(int dist=1;dist<=maxDistance;dist++)
    {
        std::vector<int> res2 = magic(index,k,dist);//it only modifies k first bits
        for(int i=0;i<res2.size();i++)
        {
            if(res2[i]<0 || res2[i] >=tableSize)
                continue;
            res.push_back(res2[i]);
            probes=probes-1;
            if(probes==0)
                return res;
        }
    }
    return res;
}
std::vector<std::pair<point*, float>>  HyperCube:: approxKNN(point p,int N)
{
    unsigned int index=hash_data_point_hypercube(p,this->coin_tossers,this->k);
    std::vector<int> vertices=getAdjacentVertices(index,this->k,this->probes,this->tableSize);
    int totalVisited=0;
    //// store info about the candidates
    point **b=new point *[N]; //k best candidates
    float *Db=new  float[N];//their distances
    for(int i=0;i<N;i++)
        Db[i]=1999999;
    for(int i=0;i<vertices.size();i++)
    {
        for(int j=0;j<table[vertices[i]].points.size();j++)
        {
            if(totalVisited > this->M)
            {
                std::vector<std::pair<point*, float>>result;
                for(int i=0;i<N;i++)
                {
                    if(Db[i]!=1999999)
                        result.push_back(std::make_pair(b[i],Db[i]));
                }
                delete[] Db;
                delete[] b;   
                return result; 
                
            }
            float tempDist=dist(p,*(table[vertices[i]].points[j]));
            for(int x=0;x<N;x++)
            {
                //points[j] must be kth best and it must not be allready inserted
                if(Db[x] > tempDist)
                {
                    Db[x]=tempDist;
                    b[x]=this->table[vertices[i]].points[j];
                    break;
                }
            }

            totalVisited+=1;
        }
    }
    //done, pass data to a vector and return
    std::vector<std::pair<point*, float>>result;
    for(int i=0;i<N;i++)
    {
        if(Db[i]!=1999999)
            result.push_back(std::make_pair(b[i],Db[i]));
    }
    delete[] Db;
    delete[] b;   
    return result; 
}

std::vector<point *> HyperCube:: approxRangeSearch(point p,float r)
{
    //first find the vertice of the data point p
    unsigned int index=hash_data_point_hypercube(p,this->coin_tossers,this->k);
    //then find the vertices you can check according to probes argument
    std::vector<int> vertices=getAdjacentVertices(index,this->k,this->probes,this->tableSize);
    int totalVisited=0;
    std::vector<point *> res;
    for(int i=0;i<vertices.size();i++)
    {
        for(int j=0;j<table[vertices[i]].points.size();j++)
        {

            if(totalVisited >=this->M)
                break;
            float tempDist=dist(p,*(table[vertices[i]].points[j]));
            if(tempDist<=r)
                res.push_back(table[vertices[i]].points[j]);
            totalVisited+=1;
        }
    }
    return res;
}
//debug
void HyperCube:: printTables()
{
    printf("numZero=%d , numOne= %d\n",numZero,numOne);

}

//destructor
HyperCube:: ~HyperCube()
{
    delete[] this->table;
    for(int i=0;i<this->k;i++)
        delete[] this->coin_tossers[i].h.v;
    delete[] this->coin_tossers;
}



std::vector<point *> fileToHyperCube(HyperCube &hq,const char *fileName,int &d)
{

    std::ifstream infile(fileName);
    if(infile.fail())
    {
        printf("ERROR : PROBLEM WITH INPUT FILE\n");
        exit(-1);
    }
    std::string line;
    //will store all the points here
    std::vector<point *>points;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        //std::cout<<"LINE = "<<line <<"\n";
        if(line.length() <2)
            break;
        point *temp=new point();
        temp->coords=strToCoords(line.data(),d,temp->id);
        temp->d=d;
        //std::cout<<"temp coords ="<<temp->coords[0] << "," << temp->coords[1]<<"\n";
        points.push_back(temp);
    }

    hq.fit(points.size(),points.data(),d);
    return points;
}