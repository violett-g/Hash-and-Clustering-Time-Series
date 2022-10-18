#include <iostream>
#include <random>
#include <math.h>       /* floor */
#include <algorithm>
#include <functional>

#include "LSH.hpp"
#include "../../timeSeries/Search/1dCurves.hpp"
#include <unordered_map>

#include <map>
#include <set>

unsigned int M=4294967291;

//choose t uniformly in [0,w)
//choose a d dim vector v whoose coords are N(0,1)
std::default_random_engine generator;

hash_funct create_hash_funct(int d,int w)
{   
    hash_funct ht;
    ht.d=d;
    ht.w=w;
    //pic v from gaussian
    std::normal_distribution<float> distribution(0.0,1.0);
    ht.v=new float[d];
    for(int i=0;i<d;i++)
        ht.v[i]= distribution(generator);

    //pic t from uniform
    std::uniform_real_distribution<float> uni_distribution(0.0,w);
    ht.t=uni_distribution(generator);

    return ht;
} 
//create k "weak" has functions  and store them in hs
//create k integers and store them in rs
amplified_hash_fucnt create_ampl_hash_funct(int k,int d,int w)
{
    amplified_hash_fucnt ahf;
    ahf.d=d;
    ahf.k=k;
    ahf.w=w;
    ahf.rs=new int[k];
    ahf.hs=new hash_funct[k];
    std::uniform_int_distribution<> distrib(-1000,1000);
    for(int i=0;i<k;i++)
    {
        ahf.hs[i]=create_hash_funct(d,w);
        //random int in [1,50]
        ahf.rs[i]=distrib(generator);

    }
    return ahf;
}
/*
for a given point p , h(p) = floor((vp+t/w))
use a linear combination of k h's allready stored
*/
int get_weak_hf(hash_funct h,point p)
{
    float dotProd=0;
    for(int i=0;i<p.d;i++)
    {
        dotProd=dotProd+h.v[i]*p.coords[i];
    }
    dotProd+=h.t;
    return floor(dotProd/h.w);
}
unsigned int myMod(int a,unsigned int b)
{
    return (a%b+b)%b;
}
unsigned int hash_data_point(amplified_hash_fucnt g,point p)
{
    unsigned int linearComb=0;
    for(int i=0;i<g.k;i++)
    {
        linearComb=myMod( 
                    myMod(linearComb,M) 
                        + 
                    myMod(
                        myMod(get_weak_hf(g.hs[i],p),M)
                        *
                        myMod(g.rs[i],M)
                    ,M)
                ,M);
    }
    
    return linearComb;
}
/*
BELLOW IS CODE FOR LSH
*/

//constructor only deals with hyperparms
LSH::LSH(int d,int w,int k,int L,float (*distance_funct)(point a,point b))
{
    this->d=d;
    this->w=w;
    this->k=k;
    this->L=L;
    this->tables=NULL;
    this->gs=NULL;
    this->n=0;
    this->distance_funct=distance_funct;
    
}
//each point is an array of floats
//LSH gets an array of arrays of floats <-> array of points 
//this function constructs LSH
void LSH::fit(int n,point **points,int d,int tableSize)
{
    //this are non hyperparam variables
    // how many points and which points
    this->d=d;
    this->n=n;
    //default is n/8 but if data too big or if want to perform
    //accurate range search for clustering we might need a smaller table size
    if(tableSize==-1)
    {
        this->tableSize=n/8;
        if(tableSize==0) //for applications with tiny dataset
            tableSize=n;
    }
    else
        this->tableSize=tableSize;
    //create hash functions
    this->gs=new amplified_hash_fucnt[this->L];
    for(int i=0;i<this->L;i++) //L amplified hash Functs
    {
        this->gs[i]=create_ampl_hash_funct(this->k,this->d,this->w);
    }
    //create the  tables
    this->tables=new bucket *[L]; //L hash tables
    for(int i=0;i<this->L;i++)
    {
        //each table has tableSize buckets and each bucket can have a lot of points
        this->tables[i]=new bucket[this->tableSize];
    }
    

    for(int i=0;i<this->L;i++)
    {
        //go through all the points(pointers)
        for(int j=0;j<n;j++)
        {
            //hash data point according to the hash function corresp
            //to the i-th table and put the point in the bucket
          
            //unsigned int index=hash_data_point(gs[i],*points[j])%this->tableSize;
            bucket_data bd;
            bd.p=points[j];
            bd.hash_id=hash_data_point(gs[i],*points[j]);
            unsigned int index=bd.hash_id%this->tableSize;
           
            this->tables[i][index].points.push_back(bd);
            //this->tables[i][index].p.points.push_back(points[j]);
            
        }
    }

}

std::vector<std::pair<point*, float>> LSH:: approxKNN(point p,int N)
{
    //printf("WILL FIND %d NEAREST NEIGHBORS\n",N);
    point **b=new point *[N]; //k best candidates
    float *Db=new  float[N];//their distances
    //we use an unordered map to quickly check if a point has allready been considered
    std::unordered_map <point *,bool>allreadyInserted;
    for(int i=0;i<N;i++)
        Db[i]=1999999;
    //for each table of buckets 
    for(int i=0;i< this -> L;i++)
    {   
        //find the coresponding bucket
        unsigned int index=hash_data_point(gs[i],p)%this->tableSize;
        //printf("index=%d,hash_funct=%d, size=%d\n",index,hash_data_point(gs[i],p), this->tables[i][index].points.size());
        //std::cout<<"BUCKET "<<i<<" INDEX="<<index<<" num points= "<<this->tables[i][index].points.size()<<"\n";
        //for each point in coresponding bucket
        for(int j=0;j< this->tables[i][index].points.size();j++)
        {
            unsigned int id1=this->tables[i][index].points[j].hash_id;
            unsigned int id2=hash_data_point(this->gs[i],p);
            if(id1!=id2)
            {
                //printf("DIFFERENT FUCKIN IDS %u,%u\n",id1,id2);
               continue;
            }
            float tempDist=this->distance_funct(p,*(this->tables[i][index].points[j].p));
            //std::cout<<"tempDist="<<tempDist<<"best="<<Db[i]<<"\n";
            //go through all the distances and check if it is k-th best dist
            //if so replace the current k-th neighbor with points[j]
            for(int x=0;x<N;x++)
            {
                //points[j] must be kth best and it must not be allready inserted
                if((Db[x] > tempDist || Db[x]==1999999) && 
                (allreadyInserted.find(this->tables[i][index].points[j].p) == allreadyInserted.end()))
                {
                    Db[x]=tempDist;
                    b[x]=this->tables[i][index].points[j].p;
                    allreadyInserted[this->tables[i][index].points[j].p]=true;
                    break;
                }
            }
        }
    }
    std::vector<std::pair<point*, float>>result;
    int fuckinSize=0;
    for(int i=0;i<N;i++)
    {
        if(Db[i]!=1999999)
        {
            result.push_back(std::make_pair(b[i],Db[i]));
            fuckinSize+=1;
        }
    }
    delete[] Db;
    delete[] b;   
    //printf("FOUND %d FUCKIN NEIGHBORS \n",fuckinSize)
    return result;   
}
 std::vector<point *> LSH:: approxRangeSearch(point p,float r)
{
    std::vector<point *> pointsInRange;
    //for each table of buckets
    for(int i=0;i<this->L;i++)
    {
        //find the coresponding bucket
        unsigned int index=hash_data_point(gs[i],p)%this->tableSize;
        //std::cout<<"BUCKET "<<i<<" INDEX="<<index<<" num points= "<<this->tables[i][index].points.size()<<"\n";
        //for each point in coresponding bucket
        for(int j=0;j< this->tables[i][index].points.size();j++)
        {
            point *currentPoint=this->tables[i][index].points[j].p;
            //make sure we haven't allready found it in some other bucket
            if
            (
                std::find(pointsInRange.begin(), pointsInRange.end(),currentPoint) 
                != pointsInRange.end()
            )
            {
                continue;
            }    
            float tempDist=this->distance_funct(p,*(this->tables[i][index].points[j].p));
            //printf("tempDist=%f r=%f\n",tempDist,r);
            if (tempDist <= r)
                pointsInRange.push_back(this->tables[i][index].points[j].p);
            //optional early stoppage condition up to the user to choose
            //if((pointsInRange.size() > 20 *this->L ) && (earlyStop==true))
            //       return pointsInRange;
         
        }
        //printf("points in Range size= %d, 20*L =%d\n",pointsInRange.size(),20*this->L);
        //return pointsInRange;
    }

    return pointsInRange;
}

//destructor
LSH::~LSH()
{
    for(int i=0;i<this->L;i++)
    {
        delete[] this->gs[i].rs;
        for(int j=0;j<this->k;j++)
            delete[] this->gs[i].hs[j].v;
        delete[] this->gs[i].hs;
        delete[] this->tables[i];
        
    }
    delete[] this->gs;
    delete[] this->tables;

}
void LSH::printTables()
{
    for(int i=0;i<this->L;i++)
    {
        std::cout<<"TABLE " << i <<": \n" ;
        for(int j=0;j<this->tableSize;j++)
            for(int a=0;a <this->tables[i][j].points.size();a++)
                for(int c=0; c<this->tables[i][j].points[a].p->d; c++){
                    std::cout<<this->tables[i][j].points[a].p->coords[c]<<" ";
                }
                std::cout << std::endl;

    }
}
std::vector<point *> fileToLSH(LSH &lsh,const char *fileName,int &d)
{

    std::ifstream infile(fileName);
    if(infile.fail())
    {
        printf("ERROR : PROBLEM WITH INPUT FILE %s\n",fileName);
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
        points.push_back(temp);
    }

    lsh.fit(points.size(),points.data(),d,points.size()/2);


    return points;
}




std::vector<std::pair<point*, float>> LSH:: frechetKNN(point p,int N)
{
    //printf("WILL FIND %d NEAREST NEIGHBORS\n",N);
    point **b=new point *[N]; //k best candidates
    float *Db=new  float[N];//their distances
    //we use an unordered map to quickly check if a point has allready been considered
    std::unordered_map <point *,bool>allreadyInserted;
    for(int i=0;i<N;i++)
        Db[i]=1999999;
    //for each table of buckets 
    for(int i=0;i< this -> L;i++)
    {   
       // std::cout << "L " << i << std::endl;
        //find the coresponding bucket
        unsigned int index=hash_data_point(gs[i],p)%this->tableSize;
        //printf("index=%d,hash_funct=%d, size=%d\n",index,hash_data_point(gs[i],p), this->tables[i][index].points.size());
        //std::cout<<"BUCKET "<<i<<" INDEX="<<index<<" num points= "<<this->tables[i][index].points.size()<<"\n";
        //for each point in coresponding bucket

        for(int j=0;j< this->tables[i][index].points.size();j++)
        {
            int id1=this->tables[i][index].points[j].hash_id;
            int id2=hash_data_point(this->gs[i],p);
            if(id1!=id2)
            {
                //printf("DIFFERENT FUCKIN IDS %d,%d\n",id1,id2);
               continue;
            }
            Curve query_c = pointToCurve(&p);
            point d = *(this->tables[i][index].points[j].p);
            Curve data_c = pointToCurve(&d);

            float tempDist =  Frechet::Continuous::distance(query_c,data_c).value;
            // float tempDist=this->distance_funct(p,*(this->tables[i][index].points[j].p));
            //std::cout<<"tempDist="<<tempDist<<"best="<<Db[i]<<"\n";
            //go through all the distances and check if it is k-th best dist
            //if so replace the current k-th neighbor with points[j]
            for(int x=0;x<N;x++)
            {
                //points[j] must be kth best and it must not be allready inserted
                if((Db[x] > tempDist || Db[x]==1999999) && 
                (allreadyInserted.find(this->tables[i][index].points[j].p) == allreadyInserted.end()))
                {
                    Db[x]=tempDist;
                    b[x]=this->tables[i][index].points[j].p;
                    allreadyInserted[this->tables[i][index].points[j].p]=true;
                    break;
                }
            }
        }
    }
    std::vector<std::pair<point*, float>>result;
    int fuckinSize=0;
    for(int i=0;i<N;i++)
    {
        if(Db[i]!=1999999)
        {
            result.push_back(std::make_pair(b[i],Db[i]));
            fuckinSize+=1;
        }
    }
    delete[] Db;
    delete[] b;   
    //printf("FOUND %d FUCKIN NEIGHBORS \n",fuckinSize);
    return result;   
}

