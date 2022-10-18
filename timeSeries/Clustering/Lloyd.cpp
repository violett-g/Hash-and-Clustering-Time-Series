#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <random>
#include <map>
#include <cstring>
#include "Lloyd.hpp"
#include "../../Common/Common.hpp"
#include <algorithm>    // std::min_element, std::max_element


bool AllreadyChoosen(int * indexesChoosen,int index,int k)
{
    for(int i=0;i<k;i++)
        if(indexesChoosen[i]==index)
            return true;
    return false;
}
//this function picks k centroids randomly 
//it is left for comparison
//on the input_small_id file it performs slightly worse , but not too bad for small k
float ** getCentroids(std::vector<point *> points,int k)
{
    float ** centroids = new float*[k];
    int *indexesChoosen=new int[k];//to make sure we don't choose the same twice

    for(int i=0;i<k;i++)
        indexesChoosen[i]=-1;

    int numCentroidsChoosen=0;
    while(numCentroidsChoosen<k)
    {
        int index=std::rand() % points.size();
        while(AllreadyChoosen(indexesChoosen,index,k)==true)
            index=std::rand() % points.size();
        centroids[numCentroidsChoosen]=new float[points[0]->d];
        for(int a=0;a<points[0]->d;a++)
            centroids[numCentroidsChoosen][a]=points[index]->coords[a];
        //store the choosen index
        indexesChoosen[numCentroidsChoosen]=index;
        numCentroidsChoosen+=1;
    }
    return centroids;
}



float **getCentroidsPlusPlus(std::vector<point *> points,int k)
{
    //FIRST SET UP GENERATORS
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    gen.seed(185);
    std::uniform_int_distribution<> distrib(0,points.size()-1);

    //data structures we will use
    std::vector <float *> centroids(k);
    int *indexesChoosen=new int[k];
    for(int i=0;i<k;i++)
        indexesChoosen[i]=-1;
    std::vector<point *> nonCentroids(points);


    //generate the first centroid randomly
    int index=distrib(gen);
    centroids[0]=new float[points[index]->d];
    for(int i=0;i<points[index]->d;i++)
        centroids[0][i]=points[index]->coords[i];
    indexesChoosen[0]=index;
    nonCentroids.erase(nonCentroids.begin()+index);
    int t=1;

    while(t<k)
    {
        //step 2)Compute Ds  for all non centroids
        std::vector<float> D(nonCentroids.size());
        for(int i=0;i<nonCentroids.size();i++)
        {
            float bestDist=dist(*points[i],centroids[0]);

            //go through all the centroids and find the min dist
            for(float *c :centroids )
            {      
                if(c==nullptr)
                    continue;
                float tempDist=dist(*points[i],c);
                if(tempDist <bestDist)
                    bestDist=tempDist;
            }
            D[i]=bestDist;
        }
        //optional) normalize all the Ds, we didn't do it

        //step3)compute Ps according to Ds
        std::vector<float> P(nonCentroids.size());
        for(int i=0;i<nonCentroids.size();i++)
        {
            P[i]=0.0;
            for(int j=0;j<=i;j++)
            {
                P[i]+=D[j]*D[j];
            }
        }
        //pick x uniformly
        std::default_random_engine generator;
        float maxP= *std::max_element(P.begin(), P.end());
        float minP= *std::min_element(P.begin(),P.end());
        std::uniform_real_distribution<double> distribution(minP,maxP);
        float x=distribution(generator);
        //now choose point i from non centroids such that
        // P[i-1]<x<=P[i]
        for(int i=1;i<nonCentroids.size();i++)
        {
            if(P[i-1]<x<=P[i]) //nonCentroids[i] will now be a centroid
            {
                centroids[t]=new float[nonCentroids[i]->d];
                for(int j=0;j<nonCentroids[i]->d;j++)
                    centroids[t][j]=nonCentroids[i]->coords[j];
                indexesChoosen[t]=i;
                nonCentroids.erase(nonCentroids.begin()+i);
                break;
            }
        }
        t++;

    }

    float **result=new float*[k];
    for(int i=0;i<k;i++)
        result[i]=centroids[i];
    delete[] indexesChoosen;
    return result;


    
}

 
std::vector<cluster> KMeans_Exact(std::vector<point *> points,int k)
{
    float **centroids=getCentroidsPlusPlus(points,k);
     
    //if dependencies[i]=j that means the i-th point is in cluster j
    std::vector<int > dependencies(points.size());
    
    //now that centroids are choosen procceed to the algorithm
    int numEpochs=100;
    while(numEpochs)
    {
        //printf("epoch = %d\n",100-numEpochs+1);
        //assignment step 
        /*
        assign each point to it's nearest center(Brute force)
        */

       //find best centroid for each point
       int numDifferenet=0;
       for(int i=0;i<points.size();i++)
       {
           float *bestCentroid=centroids[0];
           int bestCentroidIndex=0;
           float bestDist=dist(*points[i],centroids[0]);
           for(int j=1;j<k;j++)
           {
               float tempDist=dist(*points[i],centroids[j]);
               if(tempDist <bestDist)
               {
                   bestDist=tempDist;
                   bestCentroid=centroids[j];
                   bestCentroidIndex=j;
               }
           }
           int prevDep=dependencies[i];
           dependencies[i]=bestCentroidIndex;
           if(prevDep!=dependencies[i])
            numDifferenet+=1;
       }



        //update
        /*
        compute the new centers of each cluster
        */ 

        for(int i=0;i<k;i++)
        {   //for each cluster
            //go through all the points , and find the ones dependent on it
            int *meanCoords=new int[points[0]->d];
            int totalPointsInCluster=0;
            for(int j=0;j<points[0]->d;j++)
                meanCoords[j]=0.0;
            //meancoords[cluster i]=(0,0,0,....0) initially
            //go through all the points
            for(int j=0;j<points.size();j++)
            {
                //found point j that bellongs to cluster i
                if(dependencies[j]==i)
                {
                    totalPointsInCluster+=1;
                    for(int a=0;a<points[0]->d;a++)
                        meanCoords[a]+=points[j]->coords[a];
                }

            }
            //delete[] centroids[i];
            //centroids[i]=new float[points[0]->d];
            for(int a=0;a<points[0]->d;a++)
                centroids[i][a]=meanCoords[a]/totalPointsInCluster;
            delete[] meanCoords;
        }
        
        
        numEpochs--;
        /*algorithm is done store the result */
        if(numEpochs==0 || numDifferenet==0)
        {
            std::vector<cluster> result;
            for(int i=0;i<k;i++)
            {
                cluster temp;
                temp.centroidCoords=centroids[i];
                for(int j=0;j<points.size();j++)
                {    
                    if(dependencies[j]==i)
                        temp.points.push_back(points[j]);
                }
                result.push_back(temp);
            }
            
            delete[] centroids;
            return result;

        }
    }
} 



  
//returns the sillhouete index for a given point
float getSilhouette(std::vector<cluster> clusters,point p,float *centroidCoords)
{
    float bestDist=1999999;
    int bestIndex=-1;
    float a=0;
    for(int i=0;i<clusters.size();i++)
    {
        if(dist(clusters[i].centroidCoords,centroidCoords,p.d)==0)
        {
            for(int j=0;j<clusters[i].points.size();j++)
                a+=dist(*clusters[i].points[j],p);
            a=a/clusters[i].points.size();
            continue;
        }
        float tempDist=dist(p,clusters[i].centroidCoords);
        if(tempDist<bestDist || bestDist==1999999)
        {
            bestDist=tempDist;
            bestIndex=i;
        }
    }
    float b=0;
    for(int i=0;i<clusters[bestIndex].points.size();i++)
    {
        b+=dist(*clusters[bestIndex].points[i],p);
    }
    b=b/clusters[bestIndex].points.size();
    float max;
    if(b > a)
        max=b;
    else
        max=a;
    return (b-a)/max;
}   
std::vector<float>  getSilhouettes(std::vector<cluster> clusters,std::vector<point *>points)
{
    float totalSilhouette=0;
    std::vector<float> allShillouetes;
    int numPoints=0;
    for(int i=0;i<clusters.size();i++)
    {
        float sk=0;
        for(int j=0;j<clusters[i].points.size();j++)
        {   
            float sil=getSilhouette(clusters,*clusters[i].points[j],clusters[i].centroidCoords);
            sk+=sil;
            totalSilhouette+=sil;
        }
        sk=sk/clusters[i].points.size();
        allShillouetes.push_back(sk);
    }
    totalSilhouette=totalSilhouette/points.size();

    allShillouetes.push_back(totalSilhouette); //totalSilhouetete is at the end
   
    return allShillouetes;

}
//0(k^2 ) complexity , returns the minimum distance between centroids
float getMinDist(float **centroids,int k,int d)
{
    float minDist=1999999;
    for(int i=0;i<k;i++)
    {
        for(int j=0;j<k;j++)
        {
            if(j==i)
                continue;
            float tempDist=dist(centroids[i],centroids[j],d);
            if(tempDist < minDist)
                minDist=tempDist;

        }
    }
    return minDist;
}
