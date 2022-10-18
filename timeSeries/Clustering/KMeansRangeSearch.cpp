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
float **getCentroidsPlusPlus(std::vector<point *> points,int k);
float getMinDist(float **centroids,int k,int d);


std::vector<cluster> KMeans_RangeSearch(LSH &lsh,std::vector<point *> points,int k)
{
    float **centroids=getCentroidsPlusPlus(points,k);

    
    //if dependencies[i]=j that means the i-th point is in cluster j
    std::vector<int > dependencies(points.size());
    
    
    //now that centroids are choosen procceed to the algorithm
    int numEpochs=100;
    while(numEpochs)
    {
       float r=getMinDist(centroids,k,points[0]->d)/2.0;
       float threshold=10.0*r;
       std::unordered_map<point *,std::vector<int>> points_in_balls;
       std::unordered_map<point *,int>assigned;
       int *numPointsPerBall=new int[k];
       for(int i=0;i<k;i++)
            numPointsPerBall[i]=0;
        int iters=0;
        int numPointsAssigned=0;

        //FIRST CLUSTER AS MANY POINTS AS POSSIBLE WITH RANGE SEARCH
        while(r<threshold)
        {
            int numUpdatedBalls=0;
            //GO THROUGH EACH CENTROID AND FIND A BALL OF RANGE R
            for(int i=0;i<k;i++)
            {
                //first range search to get the points in a ball near a centroid
                point temp;
                temp.coords=centroids[i];
                temp.d=points[0]->d;
                temp.id=-1;
                std::vector<point *> ballPoints=lsh.approxRangeSearch(temp,r);

                 //check if ball is updated
                if(ballPoints.size() > numPointsPerBall[i])
                {
                    numPointsPerBall[i]=ballPoints.size();
                    numUpdatedBalls+=1;
                }

                //FOR EACH POINT IN THAT BALL THAT HAS NOT BEEN ASSIGNED
                //MARK IT FOR ASSIGNMENT
                for(int j=0;j<ballPoints.size();j++)
                {
                    //if it is assigned in a previous step don't touch it
                    if(assigned.find(ballPoints[j])!=assigned.end())
                        continue;
                    if(points_in_balls.find(ballPoints[j])==points_in_balls.end())
                        numPointsAssigned+=1;
                    points_in_balls[ballPoints[j]].emplace_back(i);
                }
            }
            //FOR EVERY POINT IN >= 1 BALL
            //SEE ALL THE POSSIBLE CENTROIDS AND ASSIGN IT TO THE CLOSEST
            for(auto & it:points_in_balls)
            {
                //it might have been assigned in a previous step
                if(assigned.find(it.first)!=assigned.end())
                        continue;
                point *p=it.first;
                float bestDist=dist(*p,centroids[it.second[0]]);
                int bestCentroidIndex=it.second[0];
                for(int j=1;j<it.second.size();j++)
                {
                    float tempDist=dist(*p,centroids[it.second[j]]);
                    if(tempDist<bestDist)
                    {
                        bestDist=tempDist;
                        bestCentroidIndex=it.second[j];
                    }
                }
                assigned[it.first]=bestCentroidIndex;
                numPointsAssigned+=1;
                //result[bestCentroidIndex].points.emplace_back(p);

            }
            
            iters+=1;
            r=r*2;
            if(iters>=4 && numUpdatedBalls <0.5*k )
                break;
        }

        //NOW ASSIGN THE POINTS IN O(N)
        //THOSE THAT HAVEN'T BEEN ASSIGNED WILL BE ASSIGNED USING TRADITIONAL LLOYD'S ALGO
        int numDifferent=0;
        for(int i=0;i<points.size();i++)
        {
            //ALLREADY ASSIGNED WE DON'T NEED TO COMPUTE EUCLIDEAN DISTANCE
            if(assigned.find(points[i])!=assigned.end())
            {
                int prevDep=dependencies[i];
                dependencies[i]=assigned[points[i]];
                if(prevDep!=dependencies[i])
                    numDifferent+=1;

            }
            //NOT ALLREADY ASSIGNED
            else
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
                    numDifferent+=1;
            }
        }
        numEpochs--;
        delete[] numPointsPerBall;
         /*algorithm is done store the result */
        if(numEpochs==0 || numDifferent==0)
        {
            std::vector<cluster> result;
            for(int i=0;i<k;i++)
            {
                cluster temp;
                temp.centroidCoords=centroids[i];
                for(int j=0;j<points.size();j++)
                {    
                    if(dependencies[j]==i)
                        temp.points.emplace_back(points[j]);
                }
                result.emplace_back(temp);
            }
            //for(int i=0;i<centroids.size();i++)
            
            delete[] centroids;
            return result;

        }
    }
}
        
       
std::vector<cluster> KMeans_RangeSearch(HyperCube &hq,std::vector<point *> points,int k)
{
    float **centroids=getCentroidsPlusPlus(points,k);

    
    //if dependencies[i]=j that means the i-th point is in cluster j
    std::vector<int > dependencies(points.size());
    
   
    //now that centroids are choosen procceed to the algorithm
    int numEpochs=100;
    int tolerance=0.02*points.size();
    while(numEpochs)
    {
       float r=getMinDist(centroids,k,points[0]->d)/2.0;
       float threshold=10.0*r;
       std::unordered_map<point *,std::vector<int>> points_in_balls;
       std::unordered_map<point *,int>assigned;
       int *numPointsPerBall=new int[k];
       for(int i=0;i<k;i++)
            numPointsPerBall[i]=0;
        int iters=0;
        int numPointsAssigned=0;

        //FIRST CLUSTER AS MANY POINTS AS POSSIBLE WITH RANGE SEARCH
        while(r<threshold)
        {
            int numUpdatedBalls=0;
            //GO THROUGH EACH CENTROID AND FIND A BALL OF RANGE R
            for(int i=0;i<k;i++)
            {
                //first range search to get the points in a ball near a centroid
                point temp;
                temp.coords=centroids[i];
                temp.d=points[0]->d;
                temp.id=-1;
                std::vector<point *> ballPoints=hq.approxRangeSearch(temp,r);

                 //check if ball is updated
                if(ballPoints.size() > numPointsPerBall[i])
                {
                    numPointsPerBall[i]=ballPoints.size();
                    numUpdatedBalls+=1;
                }

                //FOR EACH POINT IN THAT BALL THAT HAS NOT BEEN ASSIGNED
                //MARK IT FOR ASSIGNMENT
                for(int j=0;j<ballPoints.size();j++)
                {
                    //if it is assigned in a previous step don't touch it
                    if(assigned.find(ballPoints[j])!=assigned.end())
                        continue;
                    if(points_in_balls.find(ballPoints[j])==points_in_balls.end())
                        numPointsAssigned+=1;
                    points_in_balls[ballPoints[j]].emplace_back(i);
                }
            }
            //FOR EVERY POINT IN >= 1 BALL
            //SEE ALL THE POSSIBLE CENTROIDS AND ASSIGN IT TO THE CLOSEST
            for(auto & it:points_in_balls)
            {
                //it might have been assigned in a previous step
                if(assigned.find(it.first)!=assigned.end())
                        continue;
                point *p=it.first;
                float bestDist=dist(*p,centroids[it.second[0]]);
                int bestCentroidIndex=it.second[0];
                for(int j=1;j<it.second.size();j++)
                {
                    float tempDist=dist(*p,centroids[it.second[j]]);
                    if(tempDist<bestDist)
                    {
                        bestDist=tempDist;
                        bestCentroidIndex=it.second[j];
                    }
                }
                assigned[it.first]=bestCentroidIndex;
                numPointsAssigned+=1;
                //result[bestCentroidIndex].points.emplace_back(p);

            }
            
            iters+=1;
            r=r*2;
            if(iters>=4 && numUpdatedBalls <0.5*k )
                break;
        }
        //NOW ASSIGN THE POINTS IN O(N)
        //THOSE THAT HAVEN'T BEEN ASSIGNED WILL BE ASSIGNED USING TRADITIONAL LLOYD'S ALGO
        int numDifferent=0;
        for(int i=0;i<points.size();i++)
        {
            //ALLREADY ASSIGNED
            if(assigned.find(points[i])!=assigned.end())
            {
                int prevDep=dependencies[i];
                dependencies[i]=assigned[points[i]];
                if(prevDep!=dependencies[i])
                    numDifferent+=1;

            }
            //NOT ALLREADY ASSIGNED
            else
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
                    numDifferent+=1;
            }
        }
 
        numEpochs--;
        delete[] numPointsPerBall;
         /*algorithm is done store the result */
        if(numEpochs==0 || numDifferent<=tolerance)
        {
            std::vector<cluster> result;
            for(int i=0;i<k;i++)
            {
                cluster temp;
                temp.centroidCoords=centroids[i];
                for(int j=0;j<points.size();j++)
                {    
                    if(dependencies[j]==i)
                        temp.points.emplace_back(points[j]);
                }
                result.emplace_back(temp);
            }
            
            delete[] centroids;
            return result;

        }
    }
}
