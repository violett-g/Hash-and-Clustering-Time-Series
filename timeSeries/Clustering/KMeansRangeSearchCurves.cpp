#include "meanCurve.hpp"
#include <random>
#include <algorithm>    // std::min_element, std::max_element
discrete2DCurve **getCentroidsPlusPlus(std::vector<discrete2DCurve *> curves,int k);
float getMinDistFrechet(discrete2DCurve **centroids,int k);
bool special_compare(const discrete2DCurve *&a, const discrete2DCurve *&b)
{
    // match with wild card
    if(a->id==b->id)
        return true;
    return false;
}
std::vector<clusterOfCurves> KMeans_RangeSearch(
    LSH **lshs,grid2D *grids,std::vector<discrete2DCurve *> curves,
    int k,int vectorSize,int L,float padValue)
{
    //kmeans++ initialisation
    discrete2DCurve **centroids=getCentroidsPlusPlus(curves,k);

    //if dependencies[i]=j that means the i-th point is in cluster j
    std::vector<int > dependencies(curves.size());
    
    //now that centroids are choosen procceed to the algorithm
    int numEpochs=100;//for testing
    float tol=curves.size()*0.15;
    while(numEpochs)
    {
        /*
        1st step assignment, range search, points that didn't get assigned, 
        get assigned in the traditional way
        */

        float r=getMinDistFrechet(centroids,k)/2;
;
        float threshold=10.0*r;

        std::unordered_map<discrete2DCurve *,std::vector<int>> curves_in_balls;
        std::unordered_map<discrete2DCurve *,int>assigned;
        int *numCurvesPerBall=new int[k];
        for(int i=0;i<k;i++)
            numCurvesPerBall[i]=0;
        int iters=0;
        int numCurvesAssigned=0;
        while(r<threshold)
        {
            int numUpdatedBalls=0;
            //GO THROUGH EACH CENTROID AND FIND A BALL OF RANGE r
            for(int i=0;i<k;i++)
            {   std::vector<discrete2DCurve *> ballCurves;
                //first range search to get the points in a ball near a centroid
                for(int j=0;j<L;j++)
                {
                    point *x=new point();
                    x->d=vectorSize;
                    
                    x->coords=TS_to_Vector(*centroids[i],grids[j],vectorSize,padValue);
                    x->associatedCurve=centroids[i];
                    //printf("starting range search\n");
                    std::vector<point *> psTemp=lshs[j]->approxRangeSearch(*x,r);
                    for(auto & t: psTemp)
                    {
                        if(std::find(ballCurves.begin(),ballCurves.end(),t->associatedCurve)!=ballCurves.end())
                            continue;
                        else 
                            ballCurves.push_back(t->associatedCurve);
                    }
            
                    
                    delete[] x->coords;
                    delete x;
                    
                }
                 //check if ball is updated
                if(ballCurves.size() > numCurvesPerBall[i])
                {
                        numCurvesPerBall[i]=ballCurves.size();
                        numUpdatedBalls+=1;
                }
                //FOR EACH CURVE IN THAT BALL THAT HAS NOT BEEN ASSIGNED
                //MARK IT FOR ASSIGNMENT
                for(int j=0;j<ballCurves.size();j++)
                {
                    //if it is assigned in a previous step don't touch it
                    if(assigned.find(ballCurves[j])!=assigned.end())
                        continue;
                    if(curves_in_balls.find(ballCurves[j])==curves_in_balls.end())
                        numCurvesAssigned+=1;
                    curves_in_balls[ballCurves[j]].emplace_back(i);
                }
            }
            
            //FOR EVERY POINT IN >= 1 BALL
            //SEE ALL THE POSSIBLE CENTROIDS AND ASSIGN IT TO THE CLOSEST
            for(auto & it:curves_in_balls)
            {
                //it might have been assigned in a previous step
                if(assigned.find(it.first)!=assigned.end())
                        continue;
                discrete2DCurve *p=it.first;
                float bestDist=getFrechet(*p,*centroids[it.second[0]]);
                int bestCentroidIndex=it.second[0];
                for(int j=1;j<it.second.size();j++)
                {
                    float tempDist=getFrechet(*p,*centroids[it.second[j]]);
                    if(tempDist<bestDist)
                    {
                        bestDist=tempDist;
                        bestCentroidIndex=it.second[j];
                    }
                }
                assigned[it.first]=bestCentroidIndex;
                numCurvesAssigned+=1;
                //result[bestCentroidIndex].points.emplace_back(p);

            }
            if(assigned.size()>=0.7*curves.size())
                break;
             if(numUpdatedBalls <0.5*k )
                break;
            iters+=1;
            r=r*2;
        }

        //NOW ASSIGN THE POINTS IN O(N) (WE DONT CHECK THEIR DISTANCE TO CENTROIDS)
        //THOSE THAT HAVEN'T BEEN ASSIGNED WILL BE ASSIGNED USING TRADITIONAL LLOYD'S ALGO
        int numDifferent=0;
        int numEscaped=0;
        for(int i=0;i<curves.size();i++)
        {
            //ALLREADY ASSIGNED WE DON'T NEED TO COMPUTE FRECHET DISTANCE
            if(assigned.find(curves[i])!=assigned.end())
            {
                int prevDep=dependencies[i];
                dependencies[i]=assigned[curves[i]];
                if(prevDep!=dependencies[i])
                    numDifferent+=1;
                numEscaped+=1;

            }
            else
            {
                discrete2DCurve *bestCentroid=centroids[0];
                int bestCentroidIndex=0;
                
                float bestDist=getFrechet(*centroids[0],*curves[i]);
                //float bestDist=dist(*points[i],centroids[0]);
                for(int j=1;j<k;j++)
                {
                    float tempDist=getFrechet(*centroids[j],*curves[i]);
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

        /*
        
       2nd step update with mean  curve (frechet)
       */ 
      for(int i=0;i<k;i++)
      {
          //for each cluster, store all the curves in a vector
          std::vector<discrete2DCurve *> clustersCurves;
          for(int j=0;j<curves.size();j++)
          {
            if(dependencies[j]==i)
            {
                clustersCurves.push_back(curves[j]);
            }  
          }
          if(clustersCurves.size()>=2)
          {
            int num_of_curves=clustersCurves.size();  
            for(int j=0;j<centroids[i]->d;j++)
            {
                delete[] centroids[i]->vertices[j].coords;
            }
            delete[] centroids[i]->vertices;
            delete centroids[i];

            treeNode *t=initialize_tree(clustersCurves);
            centroids[i]=meanOfNcurves(*t);  
            deallocateTree(t,num_of_curves);

          }  

      }
      numEpochs--;
      delete[] numCurvesPerBall;
      /*algorithm is done store the result */
      if(numEpochs==0 || numDifferent<=tol)
      {
        std::vector<clusterOfCurves> result;
        for(int i=0;i<k;i++)
        {
            clusterOfCurves temp;
            temp.meanCurve=new discrete2DCurve();
            temp.meanCurve->d=centroids[i]->d;
            temp.meanCurve->vertices=new point[temp.meanCurve->d];
            for(int j=0;j<temp.meanCurve->d;j++)
            {
                temp.meanCurve->vertices[j].d=2;
                temp.meanCurve->vertices[j].coords=new float[2];
                temp.meanCurve->vertices[j].coords[0]=centroids[i]->vertices[j].coords[0];
                temp.meanCurve->vertices[j].coords[1]=centroids[i]->vertices[j].coords[1];

            }
            //centroids[i];
            for(int j=0;j<curves.size();j++)
            {    
                if(dependencies[j]==i)
                {
                    temp.curves.push_back(curves[j]);
                }
            }
            result.push_back(temp);
        }
        //CLEAN UP MEMORY
        for(int i=0;i<k;i++)
        {
            for(int j=0;j<centroids[i]->d;j++)
                delete[] centroids[i]->vertices[j].coords;
            delete[] centroids[i]->vertices;
            delete centroids[i];
        }
        delete[] centroids;
        return result;

        }

    }
    
}
//0(k^2 ) complexity , returns the minimum distance between centroids
float getMinDistFrechet(discrete2DCurve **centroids,int k)
{
    float minDist=1999999;
    for(int i=0;i<k;i++)
    {
        for(int j=0;j<k;j++)
        {
            if(j==i)
                continue;
            float tempDist=getFrechet(*centroids[i],*centroids[j]);
            if(tempDist < minDist || minDist==1999999)
                minDist=tempDist;

        }
    }
    return minDist;
}