#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <cstring>
#include "../../Common/Common.hpp"
#include "../../NN/LSH/LSH.hpp"
#include "../../NN/Proj/Proj.hpp"
#include <ctime>
#include "2dCurves.hpp"

//THIS FUNCTION READS THE QUERRY FILE LINE BY LINE
//FOR EACH CURVE (TREATED AS A VECTOR), IT QUERIES THE APPROXIMATE NN THROUGH LSH
//AND THE BRUTE FORCE NN , IT OUTPUTS THE REQUIRED FORMAT IN FILE outFname
void dealWithTestSet(
    LSH &lsh,const char *testFName,const char *outFName,
    std::vector<point *> &data,
    int d=2,int k=4,int L=5,int N=1
    )
{

    //these are some performance metrics
    float maxError=-1.0;
    float avgError=0.0;
    int numErrors=0;

    /*
    FIRST OF ALL OPEN THE FILES AND CHECK IF THERE IS SOMETHING WRONG
    */
    std::ifstream infile(testFName);
    if(infile.fail())
    {
        printf("ERROR: FAIL IN TEST FILE\n");
        exit(-1);
    }
    std::string line;
    std::ofstream outFile(outFName);
    if(outFile.fail())
    {
        printf("ERROR: FAIL IN OUTPUT FILE \n");
        exit(-1);
    }
    float timeLSH=0.0;
    float timeTrue=0.0;
    //the files are all open go through the query file line by line 
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);

        if(line.length() <2)
            break;
        //extract point from the test set
        point *temp=new point();
        temp->d=d;
        temp->coords=strToTS(line.data(),d,temp->id);
        std::clock_t start;
        float durationLSH,durationTrue;
        //first LSH querry for APROX-NN
        start = std::clock();
        std::vector<std::pair<point*, float>> neighbors=lsh.approxKNN(*temp,N);
        durationLSH = ( std::clock() - start ) / (float) CLOCKS_PER_SEC;

        //then brute force query for NN
        start = std::clock();
        std::vector<std::pair<point*, float>> bruteForceNeighbors=findNNBruteForce(*temp,data,N,dist);
        durationTrue = ( std::clock() - start ) / (float) CLOCKS_PER_SEC;
        
        if(neighbors.size()==0)
        {
            delete[] temp->coords;
            delete temp;
            continue;
        }
        //convert the points to curves to calculate freched distances
        discrete2DCurve original=pointTo2dCurve(temp);
        discrete2DCurve aprox=pointTo2dCurve(neighbors[0].first);
        discrete2DCurve tr=pointTo2dCurve(bruteForceNeighbors[0].first);
        /* calculate the current error and see if its worse than max error*/
        float curDist=getFrechet(original,aprox)/getFrechet(original,tr);
        if(curDist > maxError)
            maxError=curDist;
        avgError+=curDist;
        numErrors+=1;
        
        
        //now we got the k best neighbors (both LSH and brute force) , print then on outfile
        outFile<<"Querry :"<<temp->id <<"\n";
        
        outFile<<"Algorithm: {LSH-Vector} \n";
       

        
        outFile<<"Approximate Nearest neighbor:"<<neighbors[0].first->id<<"\n";
        outFile<<"True Nearest neighbor:"<<bruteForceNeighbors[0].first->id<<"\n";
        outFile<<"distance Approximate: "<<getFrechet(original,aprox)<<"\n";
        outFile<<"distance True: "<<getFrechet(original,tr)<<"\n";

        //outFile<<"tLSH: "<<durationLSH<<"\n"; 
        //outFile<<"tTrue: "<<durationTrue<<"\n";
        timeLSH+=durationLSH;
        timeTrue+=durationTrue;
        outFile<<"-----------------------------------------------------------\n";
        //clean up space and move to the next line in document
       delete[] temp->coords;
       delete temp;

        for(int i=0;i<original.d;i++)
            delete[] original.vertices[i].coords;
        delete[] original.vertices;

        for(int i=0;i<aprox.d;i++)
            delete[] aprox.vertices[i].coords;
        delete[] aprox.vertices;

         for(int i=0;i<tr.d;i++)
            delete[] tr.vertices[i].coords;
        delete[] tr.vertices;        

    }
    outFile<<"tApproximateAverage: "<<timeLSH/numErrors<<"\n";
    outFile<<"tTrueAverage: "<<timeTrue/numErrors<<"\n";
    outFile<<"MAF: "<<maxError<<"\n";
    outFile.close();
    printf("MAXIMUM ERROR distLSH/distTrue =%f\n",maxError);
    printf("AVERAGE ERROR distLSH/distTrue =%f\n",avgError/numErrors);

}

/*

    SAME AS BEFORE BUT WITH HYPERCUBE INSTEAD OF LSH
*/
void dealWithTestSet(
    HyperCube &hq,const char *testFName,const char *outFName,
    std::vector<point *> &data,
    int d=2,int k=4,int L=5,int N=1
    )
{

    //these are some performance metrics
    float maxError=-1.0;
    float avgError=0.0;
    int numErrors=0;

    /*
    FIRST OF ALL OPEN THE FILES AND MAKE SURE EVERYTHING IS OK
    */
    std::ifstream infile(testFName);
    if(infile.fail())
    {
        printf("ERROR: FAIL IN TEST FILE\n");
        exit(-1);
    }
    std::string line;
    std::ofstream outFile(outFName);
    if(outFile.fail())
    {
        printf("ERROR: FAIL IN OUTPUT FILE \n");
        exit(-1);
    }
    float timeLSH=0.0;
    float timeTrue=0.0;
    //the files are all open go through the query file line by line 
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);

        if(line.length() <2)
            break;
        //extract point from the test set
        point *temp=new point();
        temp->d=d;
        temp->coords=strToTS(line.data(),d,temp->id);
       
        std::clock_t start;
        float durationLSH,durationTrue;
        //first LSH querry for aprroximate NN
        start = std::clock();
        std::vector<std::pair<point*, float>> neighbors=hq.approxKNN(*temp,N);
        if(neighbors.size()==0)
        {
            delete[] temp->coords;
            delete temp;
            continue;
        }
        durationLSH = ( std::clock() - start ) / (float) CLOCKS_PER_SEC;

        //then brute force query for kNN
        start = std::clock();
        std::vector<std::pair<point*, float>> bruteForceNeighbors=findNNBruteForce(*temp,data,N,dist);
        durationTrue = ( std::clock() - start ) / (float) CLOCKS_PER_SEC;
        //convert the points to curves to calculate freched distances
        discrete2DCurve original=pointTo2dCurve(temp);
        discrete2DCurve aprox=pointTo2dCurve(neighbors[0].first);
        discrete2DCurve tr=pointTo2dCurve(bruteForceNeighbors[0].first);
        /* calculate the current error and see if its worse than max error*/
        if(neighbors.size()>0)
        {
            float curDist=getFrechet(original,aprox)/getFrechet(original,tr);
            if(curDist > maxError)
                maxError=curDist;
            avgError+=curDist;
            numErrors+=1;
        }
        //now we got the k best neighbors (both LSH and brute force) , print then on outfile
        outFile<<"Querry :"<<temp->id <<"\n";
        
        outFile<<"Algorithm: {Hypercube-Vector} \n";
       

        
        outFile<<"Approximate Nearest neighbor:"<<neighbors[0].first->id<<"\n";
        outFile<<"True Nearest neighbor:"<<bruteForceNeighbors[0].first->id<<"\n";
        outFile<<"distance Approximate: "<<neighbors[0].second<<"\n";
        outFile<<"distance True: "<<bruteForceNeighbors[0].second<<"\n";
        
        //outFile<<"tLSH: "<<durationLSH<<"\n";
        //outFile<<"tTrue: "<<durationTrue<<"\n";
        timeLSH+=durationLSH;
        timeTrue+=durationTrue;
        outFile<<"-----------------------------------------------------------\n";
        //clean up space and move to the next line in document
        delete[] temp->coords;
        delete temp;
        for(int i=0;i<original.d;i++)
            delete[] original.vertices[i].coords;
        delete[] original.vertices;
        for(int i=0;i<aprox.d;i++)
            delete[] aprox.vertices[i].coords;
        delete[] aprox.vertices;

        for(int i=0;i<tr.d;i++)
            delete[] tr.vertices[i].coords;
        delete[] tr.vertices;        

       
    }
    outFile<<"tApproximateAverage: "<<timeLSH/numErrors<<"\n";
    outFile<<"tTrueAverage: "<<timeTrue/numErrors<<"\n";
    outFile<<"MAF: "<<maxError<<"\n";
    outFile.close();
    printf("MAXIMUM ERROR distLSH/distTrue =%f\n",maxError);
    printf("AVERAGE ERROR distLSH/distTrue =%f\n",avgError/numErrors);

}

//THIS FUNCTION IS USED IN AII)
//IT TAKES AN ARRAY OF LSH TABLES AND THEIR GRIDS, AND FOR EACH QUERY CURVE
//IT SEARCHES FOR THE APPROXIMATE NEAREST CURVE AND THE TRUE NEAREST CURVE
//OUTPUTING THE REQUIRED FORMAT. FRECHET DISTANCE IS USED FOR EVERYTHING
void dealWithTestSet
(
    LSH **lshs,grid2D *grids,const char *testFName,const char *outFName,
    std::vector<discrete2DCurve *>data,int L,int d,int vectorSize,float paddValue
)
{
    //these are some performance metrics
    float maxError=-1.0;
    float avgError=0.0;
    int numErrors=0;


    std::ifstream infile(testFName);
    if(infile.fail())
    {
        printf("ERROR: FAIL IN TEST FILE\n");
        exit(-1);
    }
    std::string line;
    std::ofstream outFile(outFName);
    if(outFile.fail())
    {
        printf("ERROR: FAIL IN OUTPUT FILE \n");
        exit(-1);
    }
    float timeLSH=0.0;
    float timeTrue=0.0;
    //the files are all open go through the query file line by line 
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);

        if(line.length() <2)
            break;
        //first of all extract curve from line
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
        //now find the approximate nearest neighbors
        std::clock_t start;
        float durationLSH,durationTrue;
        start=std::clock();
        //for each of the grids, 
        float bestDist=1999999;
        std::string bestId="";
        discrete2DCurve *aprox;
        for(int i=0;i<L;i++)
        {   
            //map the curve to a point x 
            point *x=new point();
            x->d=vectorSize;
            x->coords=TS_to_Vector(*cv,grids[i],vectorSize,paddValue);
            
            //connect x with current curve through pointer
            x->associatedCurve=cv;    
            //search for x in LSH
            std::vector<std::pair<point *,float>> approxNN=lshs[i]->approxKNN(*x,1);
            delete[] x->coords;
            delete x;
            if(approxNN.size()==0)
                continue;
            if(approxNN[0].second<bestDist || bestDist==1999999)
            {
                bestDist=approxNN[0].second;
                bestId=approxNN[0].first->associatedCurve->id;
                aprox=approxNN[0].first->associatedCurve;
        
            }
            
        }  
        durationLSH = ( std::clock() - start ) / (float) CLOCKS_PER_SEC; 
        //brute force
        start = std::clock();
        std::vector<std::pair<discrete2DCurve *, float>> bruteForceNeighbors=findNNBruteForce(*cv,data,1,getFrechet);
        durationTrue = ( std::clock() - start ) / (float) CLOCKS_PER_SEC;
        //LSH MIGHT FAIL IN THAT CASE, ASSIGN IT THE SAME AS BRUTE FORCE 
        if(bestDist==1999999)
        {
            printf("LSH FAIL\n");
            bestDist=bruteForceNeighbors[0].second;
            bestId=bruteForceNeighbors[0].first->id;
            durationLSH+=durationTrue;
        }
        /* calculate the current error and see if its worse than max error*/
        float curDist=bestDist/bruteForceNeighbors[0].second;
        if(curDist > maxError)
            maxError=curDist;
        avgError+=curDist;
        numErrors+=1;
        //could not locate any neihgourhs
        
        //now we got the k best neighbors (both LSH and brute force) , print then on outfile
        outFile<<"Querry :"<<cv->id <<"\n";
        
        outFile<<"Algorithm: {LSH_Frechet_Discrete} \n";
       

        
        outFile<<"Approximate Nearest neighbor:"<<bestId<<"\n";
        outFile<<"True Nearest neighbor:"<<bruteForceNeighbors[0].first->id<<"\n";
        outFile<<"distance Approximate: "<<bestDist<<"\n";
        outFile<<"distance True: "<<bruteForceNeighbors[0].second<<"\n";
        
        //outFile<<"tLSH: "<<durationLSH<<"\n";
        //outFile<<"tTrue: "<<durationTrue<<"\n";
        timeLSH+=durationLSH;
        timeTrue+=durationTrue;
        outFile<<"-----------------------------------------------------------\n";
        //clean up space and move to the next line in document
        delete[] temp;
        for(int i=0;i<cv->d;i++)
            delete[] cv->vertices[i].coords;
        delete[] cv->vertices;
        delete cv;
    }
    outFile<<"tApproximateAverage: "<<timeLSH/numErrors<<"\n";
    outFile<<"tTrueAverage: "<<timeTrue/numErrors<<"\n";
    outFile<<"MAF: "<<maxError<<"\n";
    outFile.close();
    
}

//QUESTION AIII
//THIS FUNCTION IS USED IN AIII)
//IT TAKES AN ARRAY OF LSH TABLES AND THEIR GRIDS, PREPROCCESSED DATA AND QUERIES 
//AND FOR EACH QUERY CURVE IT SEARCHES FOR THE APPROXIMATE NEAREST CURVE AND THE TRUE NEAREST CURVE
//OUTPUTING THE REQUIRED FORMAT. CONTINIOUS FRECHET DISTANCE IS USED FOR EVERYTHING
void dealWithTestSet
(LSH **lsh,std::vector<point*>data, std::vector<std::vector<point*>>queries_vector,std::vector<point*>queries,
  std::string outFName, int L, float duration_preprocces
)
{
    // std::cout << "DEAAAAAL" << std::endl;
    //these are some performance metrics
    float maxError=-1.0;
    float avgError=0.0;
    int numErrors=0;
    
    float timeLSH=0.0;
    float timeTrue=0.0;

    //open output file
    std::ofstream outFile(outFName);
    if(outFile.fail())
    {
        printf("ERROR: FAIL IN OUTPUT FILE \n");
        exit(-1);
    }
    
    // point *b=new point; //k best candidates
    // b->id = "dummy";
    // float Db= 1999999; 
        
    // printf("fook %d\n",queries_vector[0].size()); 
    for(int q=0; q<queries_vector[0].size(); q++)
    {
        //dummy point
        float best_dist = 1999999;
        std::string best_id ="dummy";
        
        clock_t start;
        start = std::clock();
        for(int l=0; l<L; l++){
            std::vector<std::pair<point*, float>> nearest = lsh[l]->frechetKNN(*queries_vector[l][q],1);//for each query
            if(nearest.size()==0)
                continue;
            // std::cout << "frechet returned  " << nearest[0].first-> id << "  " << nearest[0].second << std::endl;
            if(nearest[0].second<best_dist || best_dist==1999999)
            {
                best_dist=nearest[0].second;
                best_id=nearest[0].first->id;
            }
        }
        float durationLSH = ( std::clock() - start ) / (float) CLOCKS_PER_SEC;

        //then brute force query for kNN
        // printf("GOING BRUTE FORCE\n");
        start = std::clock();
        std::vector<std::pair<point*, float>> bruteForceNeighbors = FrechetNNBruteForce(*queries[q],data,1);
       // std:: cout << "True NEAREST of " << queries[q]->id << "---> " << bruteForceNeighbors[0].first->id << " , distance " <<bruteForceNeighbors[0].second<<std::endl;
        float durationTrue = ( std::clock() - start ) / (float) CLOCKS_PER_SEC;
        //std::cout << "duration true " << durationTrue << std::endl;

        //in case lsh failed to return any result than use brute force
        if(best_dist==1999999)
        {
            //printf("LSH FAIL\n");
            best_dist=bruteForceNeighbors[0].second;
            best_id=bruteForceNeighbors[0].first->id;
            durationLSH+=durationTrue;
        }

        // durationLSH+=duration_preprocces;

        //std::cout << " lsh NEAREST of " << queries_vector[0][q]->id << "---> " << best_id << " , distance " << best_dist << std::endl;
        //std::cout << "duration lsh " << durationLSH << std::endl;


        // calculate the current error and see if its worse than max error
        float cur_dist=best_dist/bruteForceNeighbors[0].second;
        if(cur_dist > maxError)
            maxError=cur_dist;
        avgError+=cur_dist;
        numErrors+=1;


        //now we got the k best neighbors (both LSH and brute force) , print then on outfile
        outFile<<"Querry :"<<queries[q]->id <<"\n";
        
        outFile<<"Algorithm: {LSH_Frechet_Continous} \n";
       

        
        outFile<<"Approximate Nearest neighbor:"<<best_id<<"\n";
        outFile<<"True Nearest neighbor:"<<bruteForceNeighbors[0].first->id<<"\n";
        outFile<<"distance Approximate: "<<best_dist<<"\n";
        outFile<<"distance True: "<<bruteForceNeighbors[0].second<<"\n";
        
        //outFile<<"tLSH: "<<durationLSH<<"\n";
        //outFile<<"tTrue: "<<durationTrue<<"\n";
        timeLSH+=durationLSH;
        timeTrue+=durationTrue;
        outFile<<"-----------------------------------------------------------\n";


        //clean up space and move to the next line in document

    }//for each query
    timeLSH+=duration_preprocces;
    outFile<<"tApproximateAverage: "<<timeLSH/numErrors<<"\n";
    outFile<<"tTrueAverage: "<<timeTrue/numErrors<<"\n";
    outFile<<"MAF: "<<maxError<<"\n";
    outFile.close();
   

}
