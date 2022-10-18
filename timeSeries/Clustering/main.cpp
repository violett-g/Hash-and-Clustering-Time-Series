#include "Lloyd.hpp"
#include "../../NN/LSH/LSH.hpp"
#include "../../NN/Proj/Proj.hpp"
#include "../../Common/Common.hpp"
#include "../Search/2dCurves.hpp"
#include <string>
#include <iostream>
#include <cstring>
#include <sstream>
#include <fstream>
#include <ctime>
#include <limits.h>
#include <random>
extern std::default_random_engine generator;

int main(int argc,char **argv)
{
    //for reproducability
    generator.seed(185);
    std::uniform_real_distribution<float> uni_distribution(1000000,INT_MAX);
    srand(185);
    float paddValue=uni_distribution(generator);

    //default variables
    int k=-1,L=3,k_LSH=4,M_H=10,k_H=3,p_H=2;
    bool complete=false;
    bool silhouete=false;
    /*
    DEAL WITH INPUT FROM ARGC,ARGV
    */
    std::string assignment;
    std::string update;
    if(argc!=11 && argc!=12 && argc!=13)
    {
        printf("NOT ENOUGH ARGUMENTS\n");
        exit(-1);
    }
    std::string fname;
    std::string config_fname;
    std::string out_fname;
    for(int i=0;i<argc;i++)
    {
        if(strcmp(argv[i],"-i")==0)
            fname=argv[i+1];
        if(strcmp(argv[i],"-c")==0)
            config_fname=argv[i+1];
        if(strcmp(argv[i],"-o")==0)
            out_fname=argv[i+1];   
        if(strcmp(argv[i],"-update")==0)
        {
            update=argv[i+1];
            update=update+" ";
            update=update+argv[i+2];
        }
        if(strcmp(argv[i],"-assignment")==0)
            assignment=argv[i+1];       
        if(strcmp(argv[i],"-complete")==0)
            complete=true;
        if(strcmp(argv[i],"-silhouette")==0)
            silhouete=true;
    }
    //NOW READ CONFIG FILE
    std::ifstream confile(config_fname);
    if(confile.fail())
    {
        printf("SOMETHING IS WRONG WITH CONFIG FILE \n");
        exit(-1);
    }
    std::string line;
    while (std::getline(confile, line))
    {
        char *ptr=strtok((char *)line.data(),":");
        if(ptr==NULL)
            continue;
        if(strcmp(ptr,"number_of_clusters")==0)
            k=atoi(strtok(NULL,":"));
        if(strcmp(ptr,"number_of_vector_hash_tables")==0)
            L=atoi(strtok(NULL,":"));
        if(strcmp(ptr,"number_of_vector_hash_functions")==0)
            k_LSH=atoi(strtok(NULL,":"));
        //fill the rest here for variables ,M_H=10,k_H=3,p_H=2
        if(strcmp(ptr,"max_number_M_hypercube")==0)
            M_H=atoi(strtok(NULL,":"));
        if(strcmp(ptr,"number_of_hypercube_dimensions")==0)
            k_H=atoi(strtok(NULL,":"));
        if(strcmp(ptr,"number_of_probes")==0)
            p_H=atoi(strtok(NULL,":"));


    }
    if(k<=1 || L<=0 || k_LSH <=0 || M_H <=0 ||k_H<=0 || p_H<=0)
    {
        printf("MISSING OR INVALID  PARAMETERS \n");
        exit(-1);
    }
    /*
        NOW MOVING ON TO THE MAIN PART , CLUSTERING
    */
    int d=-1;
    std::vector<point *> points;
    std::vector<cluster> resultPoints;
    std::vector<discrete2DCurve *> curves;
    std::vector<clusterOfCurves >resultCurves;
    std::vector<float > sils;
    double duration;

    //open file for output
    std::ofstream outFile(out_fname);
    /*
    CLASSIC K-MEANS
    */
    if(assignment=="Classic")
    {
        //SUCH AS ASSIGNMENT 1
        if(update=="Mean Vector")  
        {
            outFile<<"Assignment: Lloyd's, Update: Mean Vector\n";
            std::clock_t start;
            int d=-1;
            LSH lsh(d,4,k_LSH,L);
            //SOS REMEMBER TO CHANGE IT TO FILETOVECTOR
            points=fileToLSH(lsh,fname.c_str(),d);
            //printf("Starting Kmeans Exact\n");
            start = std::clock();
            resultPoints=KMeans_Exact(points,k);
            duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
            //printf("duration=%lf\n",duration);
            if(silhouete==true)
                sils = getSilhouettes(resultPoints,points);
        }
        //CLASSIC K MEANS ,CURVES AND NOT VECTORS
        else if(update=="Mean Frechet")
        {
            outFile<<"Assignment: LLoyd's, Update:Mean Frechet\n";
            std::ifstream infile(fname.c_str());
            if(infile.fail())
            {
                printf("ERROR: FAIL IN INPUT FILE\n");
                exit(-1);
            }
            std:: string line;
            int d=-1;
            //GO THROUGH THE INPUT FILE LINE BY LINE AND SAVE THE CURVES
            while (std::getline(infile, line))
            {
        
                discrete2DCurve * cv=new discrete2DCurve();
                float *temp=strToTS(line.data(),d,cv->id);
                //printf("d=%d\n",d);
                cv->vertices=new point[d];
                cv->d=d;
                for(int i=0;i<d;i++)
                {
                    cv->vertices[i].coords=new float[2];
                    cv->vertices[i].coords[0]=i+1;
                    cv->vertices[i].coords[1]=temp[i];
                    cv->vertices[i].d=2; 
                } 
                //test->vertices=strToCoords(line.data(),d,temp->id);
                //test->d=
                delete[] temp;
                curves.push_back(cv);
            }
            std::clock_t start;
            start = std::clock();
            resultCurves=KMeans_Exact(curves,k);
            duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
            if(silhouete==true)
                sils = getSilhouettes(resultCurves,curves);
            
             
        }
        else{
            printf("INVALID ASSIGNMENT FOR KMEANS\n");
            exit(1);
        }
        
    }
    /*
    K-MEANS USING LSH  
    CURVES ARE VIEWED AS VECTORS (SAME AS ASSIGNMENT 1)
    */
    else if(assignment=="LSH")
    {
        if(update!="Mean Vector")  
        {
            printf("INVALID ASSIGNMENT \n");
            exit(1);
        }

        outFile<<"Assingmnet:  Range Search LSH, Update: Mean Vector\n";
        int d=-1;
        LSH lsh(d,300,k_LSH,L);
        //make LSH
        points=fileToLSH(lsh,fname.c_str(),d);
          /* Perform  Kmeans using LSH*/
        std::clock_t start;
        start = std::clock();

        resultPoints=KMeans_RangeSearch(lsh,points,k);

        duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

        if(silhouete==true)
            sils= getSilhouettes(resultPoints,points);
        
    }
    /*
    K MEANS WITH LSH RANGE SEARCH
    (CURVES ARE VIEWED AS 2D DIRCRETE CURVES AND NOT VECTORS)
    */
    else if(assignment=="LSH_Frechet")
    {
         
        
        if(update!="Mean Frechet")
        {
            printf("INVALID ASSIGNMENT\n");
            exit(1);
        }
 
        /*
        put code here
        */
       outFile<<"Assignmnet:  Range Search LSH Frechet, Update: Mean Frechet\n";
       float delta=0.1; //try smaller delta
       LSH **lsh=new LSH*[L];
       grid2D *grids=new grid2D[L];
       int vectorSize=1500;
       //OPEN THE INPUT FILE
       std::ifstream infile(fname);
        if(infile.fail())
        {
            printf("ERROR: FAIL IN INPUT FILE\n");
            exit(-1);
        }
        std:: string line;
        int d=-1;
        //GO THROUGH THE INPUT FILE LINE BY LINE AND SAVE THE CURVES
        while (std::getline(infile, line))
        {
        
            discrete2DCurve * cv=new discrete2DCurve();
            float *temp=strToTS(line.data(),d,cv->id);
            //printf("d=%d\n",d);
            cv->vertices=new point[d];
            cv->d=d;
            for(int i=0;i<d;i++)
            {
                cv->vertices[i].coords=new float[2];
                cv->vertices[i].coords[0]=i+1;
                cv->vertices[i].coords[1]=temp[i];
                cv->vertices[i].d=2; 
            } 
                //test->vertices=strToCoords(line.data(),d,temp->id);
                //test->d=d;
            delete[] temp;
            curves.push_back(cv);

        }
        
        //CREATE DATA STRUCTURES
        for(int i=0;i<L;i++)
        {
            lsh[i]=new LSH(vectorSize,150,k,1,getFrechet);
            grids[i]=createGrid(delta,d);
        }
            //exit(1);
        std::unordered_map<int,std::vector<point *> > pointsPerGrid;
        //FOR EACH CURVE
        for(int i=0;i<curves.size();i++)
        {   
                
            //FOR EACH OF THE L GRIDS
            for(int j=0;j<L;j++)
            {
                    //find vector x for the curve
                point *x=new point();
                x->d=vectorSize;
                x->coords=TS_to_Vector(*curves[i],grids[j],vectorSize,paddValue);

                    //connect x with current curve through pointer
                x->associatedCurve=curves[i]; 

                    //hash x on lsh[i]
                pointsPerGrid[j].push_back(x);
            }
        }
            //exit(1);
        for(int i=0;i<L;i++)
        {
            //lsh[i]->fit(pointsPerGrid[i].size(),pointsPerGrid[i].data(),vectorSize,pointsPerGrid[i].size()/4);

             lsh[i]->fit(pointsPerGrid[i].size(),pointsPerGrid[i].data(),vectorSize,pointsPerGrid[i].size()/8);
        }
        std::clock_t start;
        start = std::clock();
        resultCurves=KMeans_RangeSearch(lsh,grids,curves,k,vectorSize,L,paddValue);
        duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        if(silhouete==true)
            sils = getSilhouettes(resultCurves,curves);
        for(int i=0;i<L;i++)
            delete lsh[i];
        delete[] lsh;
        
    }
    /*
    K MEANS USING RANGE SEARCH WITH HYPERCUBE
    CURVES ARE VIEWED AS VECTORS
    */
    else if(assignment=="Hypercube")
    {
        if(update!="Mean Vector")  
        {
            printf("INVALID OPTION\n");
            exit(1);
        }
        outFile<<"Assingment: Range Search Hypercube, Update: Mean Vector\n";
        int d=-1;
        
        HyperCube hq(k_H,M_H,p_H);
        points=fileToHyperCube(hq,fname.c_str(),d);
        std::clock_t start;

        start = std::clock();
    
        resultPoints=KMeans_RangeSearch(hq,points,k);
        duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        if(silhouete==true)
            sils= getSilhouettes(resultPoints,points);
       
    }
    else{
        printf("invalid assignment\n");
        exit(-1);
        
    }
    
    /*
    CLUSTERING DONE TIME TO OUTPUT RESULTS
    */
    if(update=="Mean Vector")
    {
        //output clusters
        for(int i=0;i<resultPoints.size();i++)
        {
            

            outFile<<"CLUSTER-"<<i+1<<" { size: "<<resultPoints[i].points.size()<<", centroid :\n[";
            for(int j=0;j<points[0]->d-1;j++)
                outFile<<resultPoints[i].centroidCoords[j]<<", ";
            outFile<<resultPoints[i].centroidCoords[points[0]->d-1]<<" ] }\n";
        }
        outFile<<"clustering_time: "<<duration<<" sec \n";
        //output silhouetes
        if(silhouete==true)
        {
            outFile<<"Silhouettes : [";
            for(int i=0; i<sils.size()-1;i++)
                outFile<<sils[i]<<", ";
            outFile<<sils[sils.size()-1]<<" ]\n";
            
        }

        /* IF COMPLETE==TRUE*/
        if(complete==true)
        {
            for(int i=0;i<resultPoints.size();i++)
            {
                outFile<<"CLUSTER-"<<i+1<<"{ [";
                for(int j=0;j<points[0]->d-1;j++)
                    outFile<<resultPoints[i].centroidCoords[j]<<", ";
                outFile<<resultPoints[i].centroidCoords[points[0]->d-1]<<" ], ";
                for(int j=0;j<resultPoints[i].points.size()-1;j++)
                    outFile<<resultPoints[i].points[j]->id<<", ";
                outFile<<resultPoints[i].points[resultPoints[i].points.size()-1]->id<<" }\n";

            }
        }
        //CLEAN UP MEMORY
        
        for(int i=0;i<points.size();i++)
        {
            delete[] points[i]->coords;
            delete points[i];
        }
        for(int i=0;i<resultPoints.size();i++)
        {
            delete[] resultPoints[i].centroidCoords;
        }
    }
    else if(update=="Mean Frechet")
    {
        //output clisters
        for(int i=0;i<resultCurves.size();i++)
        {
            
            outFile<<"CLUSTER-"<<i+1<<" { size: "<<resultCurves[i].curves.size()<<", centroid :\n[";
            for(int j=0;j<curves[0]->d-1;j++)
                outFile<<resultCurves[i].meanCurve->vertices[j].coords[1]<<", ";
            outFile<<resultCurves[i].meanCurve->vertices[curves[0]->d-1].coords[1]<<" ] }\n";
        }
        outFile<<"clustering_time: "<<duration<<" sec \n";
        //output silhouetes
        if(silhouete==true)
        {
            outFile<<"Silhouettes : [";
            for(int i=0; i<sils.size()-1;i++)
                outFile<<sils[i]<<", ";
            outFile<<sils[sils.size()-1]<<" ]\n";
        }
        /* IF COMPLETE==TRUE*/
        if(complete==true)
        {
            for(int i=0;i<resultCurves.size();i++)
            {
                outFile<<"CLUSTER-"<<i+1<<"{ [";
                for(int j=0;j<curves[0]->d-1;j++)
                    outFile<<resultCurves[i].meanCurve->vertices[j].coords[1]<<", ";
                outFile<<resultCurves[i].meanCurve->vertices[curves[0]->d-1].coords[1]<<" ] \n";
                for(int j=0;j<resultCurves[i].curves.size()-1;j++)
                    outFile<<resultCurves[i].curves[j]->id<<", ";
                outFile<<resultCurves[i].curves[resultCurves[i].curves.size()-1]->id<<" }\n";

            }
        }
        //clean up memory
        for(int i=0;i<curves.size();i++)
        {
            for(int j=0;j<curves[i]->d;j++)
            {
                delete[] curves[i]->vertices[j].coords;
            }
            delete[] curves[i]->vertices;
            delete curves[i];
        }
        for(int i=0;i<resultCurves.size();i++)
        {
            for(int j=0;j<resultCurves[i].meanCurve->d;j++)
                delete[] resultCurves[i].meanCurve->vertices[j].coords;
            delete[] resultCurves[i].meanCurve->vertices;
            delete resultCurves[i].meanCurve;
        }
    }
}
