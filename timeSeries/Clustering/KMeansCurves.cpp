#include "meanCurve.hpp"
#include <random>
#include <algorithm>    // std::min_element, std::max_element

bool AllreadyChoosen(int * indexesChoosen,int index,int k);
//random choice of centroids
discrete2DCurve ** getCentroids(std::vector<discrete2DCurve *> curves,int k)
{
    discrete2DCurve ** centroids = new discrete2DCurve*[k];
    int *indexesChoosen=new int[k];//to make sure we don't choose the same twice
    int d=curves[0]->d;
    for(int i=0;i<k;i++)
        indexesChoosen[i]=-1;

    int numCentroidsChoosen=0;
    while(numCentroidsChoosen<k)
    {
        int index=std::rand() % curves.size();
        while(AllreadyChoosen(indexesChoosen,index,k)==true)
            index=std::rand() % curves.size();
        centroids[numCentroidsChoosen]=new discrete2DCurve;
        centroids[numCentroidsChoosen]->vertices=new point[d];
        centroids[numCentroidsChoosen]->d=d;
        for(int i=0;i<d;i++)
        {   
            centroids[numCentroidsChoosen]->vertices[i].d=2;
            centroids[numCentroidsChoosen]->vertices[i].coords=new float[2];
            centroids[numCentroidsChoosen]->vertices[i].coords[0]=curves[index]->vertices[i].coords[0];
            centroids[numCentroidsChoosen]->vertices[i].coords[1]=curves[index]->vertices[i].coords[1];

        }
        indexesChoosen[numCentroidsChoosen]=index;
        numCentroidsChoosen+=1;
       /*centroids[numCentroidsChoosen]=new float[curves[0]->d];
        for(int a=0;a<curves[0]->d;a++)
            centroids[numCentroidsChoosen][a]=points[index]->coords[a];
        //store the choosen index
        indexesChoosen[numCentroidsChoosen]=index;
        numCentroidsChoosen+=1;*/
    }
    return centroids;
}
//kmeans plus plus initialisation , the same as assigment one but moddified for curves(using frechet)
discrete2DCurve **getCentroidsPlusPlus(std::vector<discrete2DCurve *> curves,int k)
{
    //FIRST SET UP GENERATORS
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    gen.seed(0);

    std::uniform_int_distribution<> distrib(0,curves.size()-1);

    //data structures we will use
    std::vector <discrete2DCurve *> centroids(k);
    int *indexesChoosen=new int[k];
    for(int i=0;i<k;i++)
        indexesChoosen[i]=-1;
    std::vector<discrete2DCurve *> nonCentroids(curves);
    //generate the first centroid randomly
    int index=distrib(gen);
    centroids[0]=new discrete2DCurve();
    centroids[0]->vertices=new point[curves[index]->d];
    centroids[0]->d=curves[index]->d;

    for(int i=0;i<curves[index]->d;i++)
    {
        centroids[0]->vertices[i].d=2;
        centroids[0]->vertices[i].coords=new float[2];
        centroids[0]->vertices[i].coords[0]=curves[index]->vertices[i].coords[0];
        centroids[0]->vertices[i].coords[1]=curves[index]->vertices[i].coords[1];

    }
    indexesChoosen[0]=index;
    nonCentroids.erase(nonCentroids.begin()+index);
    int t=1;

    while(t<k)
    {
        //step 2)Compute Ds  for all non centroids
        std::vector<float> D(nonCentroids.size());
        for(int i=0;i<nonCentroids.size();i++)
        {
            float bestDist=getFrechet(*curves[i],*centroids[0]);

            //go through all the centroids and find the min dist
            for(discrete2DCurve *c :centroids )
            {      
                if(c==nullptr)
                    continue;
                float tempDist=getFrechet(*curves[i],*c);
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
                centroids[t]=new discrete2DCurve();
                centroids[t]->d=nonCentroids[i]->d;
                centroids[t]->vertices=new point[nonCentroids[i]->d];
                for(int j=0;j<nonCentroids[i]->d;j++)
                {
                    centroids[t]->vertices[j].coords=new float[2];
                    centroids[t]->vertices[j].d=2;
                    centroids[t]->vertices[j].coords[0]=nonCentroids[i]->vertices[j].coords[0];
                    centroids[t]->vertices[j].coords[1]=nonCentroids[i]->vertices[j].coords[1];

                }
                indexesChoosen[t]=i;
                nonCentroids.erase(nonCentroids.begin()+i);
                break;
            }
        }
        t++;

    }

    discrete2DCurve **result=new discrete2DCurve*[k];
    for(int i=0;i<k;i++)
        result[i]=centroids[i];
    delete[] indexesChoosen;
    return result;


    
}
//same as assignment one only two differences
//1)frechet distance
//2)mean curve
std::vector<clusterOfCurves> KMeans_Exact(std::vector<discrete2DCurve *> curves,int k)
{
    //discrete2DCurve **centroids=getCentroids(curves,k);
    discrete2DCurve **centroids=getCentroidsPlusPlus(curves,k);

  
    //if dependencies[i]=j that means the i-th point is in cluster j
    std::vector<int > dependencies(curves.size());
    
    //now that centroids are choosen procceed to the algorithm
    int numEpochs=100;
    int tol=0.1*curves.size();
    while(numEpochs)
    {
        //assignment step 
        /*
        assign each curve to it's nearest center(Brute force)
        we use frechet distance
        */

       //find best centroid for each point
       int numDifferenet=0;
       for(int i=0;i<curves.size();i++)
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
           {
               numDifferenet+=1;
           }
       }
       /*
            now time for the update step 
            will use mean curves 

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
      /*algorithm is done store the result */
      if(numEpochs==0 || numDifferenet<tol)
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
                    temp.curves.push_back(curves[j]);
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
//returns the sillhouete index for a given point
float getSilhouette(std::vector<clusterOfCurves> clusters,discrete2DCurve p,discrete2DCurve meanCurve,int clusterIndex)
{
    float bestDist=1999999;
    int bestIndex=-1;
    float a=0;
    for(int i=0;i<clusters.size();i++)
    {
        if(i==clusterIndex)
        {  
            for(int j=0;j<clusters[i].curves.size();j++)
            {
                float frech=getFrechet(*clusters[i].curves[j],p);
                a+=frech;
            }
            a=a/clusters[i].curves.size();
            continue;
        }
        float tempDist=getFrechet(p,*clusters[i].meanCurve);
        if(tempDist<bestDist || bestDist==1999999)
        {
            bestDist=tempDist;
            bestIndex=i;
        }
    }
    float b=0;
    for(int i=0;i<clusters[bestIndex].curves.size();i++)
    {
        b+=getFrechet(*clusters[bestIndex].curves[i],p);
    }
    b=b/clusters[bestIndex].curves.size();
    float max;
    if(b > a)
        max=b;
    else
        max=a;
    return (b-a)/max;
}   
std::vector<float>  getSilhouettes(std::vector<clusterOfCurves> clusters,std::vector<discrete2DCurve *>curves)
{
    float totalSilhouette=0;
    std::vector<float> allShillouetes;
    int numPoints=0;
    for(int i=0;i<clusters.size();i++)
    {
        float sk=0;
        for(int j=0;j<clusters[i].curves.size();j++)
        {   
            //printf("TESTING %f\n",getFrechet(*clusters[i].meanCurve,*clusters[i].meanCurve));
            float sil=getSilhouette(clusters,*clusters[i].curves[j],*clusters[i].meanCurve,i);
            sk+=sil;
            //if(sk==1)
            //    printf("silhouette=1 and cluster size=%d\n",clusters[i].curves.size());
            totalSilhouette+=sil;
        }
        sk=sk/clusters[i].curves.size();
        allShillouetes.push_back(sk);
    }
    totalSilhouette=totalSilhouette/curves.size();

    allShillouetes.push_back(totalSilhouette); //totalSilhouetete is at the end
   
    return allShillouetes;

}

