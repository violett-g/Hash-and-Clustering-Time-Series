#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <cstring>
#include "Proj.hpp"
#include "../../Common/Common.hpp"
#include <ctime>






//this function will parse the test file line by line 
//for each vector it will output(on the file with name outFname):
// ----- for vector (..,...,...) k nearset ---
// ----- id 1 , id 2 , ... id k

void dealWithTestSet(
    HyperCube &hq,const char *testFName,const char *outFName,std::vector<point *> &data,
    int d=2,int k=14,int N=1,float R=100000 
    )
{
    float maxError=-1;
    float avgError=0.0;
    int numVeryBad=0;
    int numErrors=0;
    int numVeryVeryBad=0;
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
    //will store all the points here
    std::vector<point *>points;
    double timeHyperCube=0.0;double timeTrue=0.0;
    int numQueries=0;
    while (std::getline(infile, line))
    {
        numQueries+=1;
        std::istringstream iss(line);
        //std::cout<<"LINE = "<<line <<"\n";
        if(line.length() <2)
            break;
        //extract point fromthe test set
        point *temp=new point();
        temp->d=d;
        temp->coords=strToCoords(line.data(),d,temp->id);
        //std::cout<<"temp coords ="<<temp->coords[0] << "," << temp->coords[1]<<"\n";
        /*outFile<<"For point: [";
        for(int i=0;i<d;i++)
            outFile<<temp->coords[i]<<",";
        outFile<<"] :\n";*/
        //query in hq
        std::clock_t start;
        double durationHyperCube,durationTrue;

        start = std::clock();
        std::vector<std::pair<point*, float>> neighbors=hq.approxKNN(*temp,N);
        durationHyperCube = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        timeHyperCube+=durationHyperCube;

        start = std::clock();
        std::vector<std::pair<point*, float>> bruteForceNeighbors=findNNBruteForce(*temp,data,N);
        durationTrue = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        timeTrue+=durationTrue;
        /* calculate the current error and see if its worse than max error*/
        if(neighbors.size()>0)
        {
            float curDist=neighbors[0].second/dist(*(bruteForceNeighbors[0].first),*temp);
            if(curDist > maxError)
                maxError=curDist;
            avgError+=curDist;
            numErrors+=1;
            if(curDist>2)
                numVeryBad+=1;
            if(curDist>3)
                numVeryVeryBad+=1;
        }
        //now we got the k best neighbors (one for each array , print then in outfile)
        outFile<<"Querry :"<<temp->id <<"\n";
        for(int i=0;i<neighbors.size();i++)
        {  
            outFile<<"Nearest Neighbor-"<<i+1<<" : "<<neighbors[i].first->id<<"\n";
            outFile<<"distance HyperCube: "<<neighbors[i].second<<"\n";
            outFile<<"distance True: "<<dist(*(bruteForceNeighbors[i].first),*temp)<<"\n";
        }
        outFile<<"tHyperCube: "<<durationHyperCube<<"\n";
        outFile<<"tTrue: "<<durationTrue<<"\n";
    
        std::vector<point *>pointsInRange=hq.approxRangeSearch(*temp,R);
        outFile<<"\n R-near neighbors (" << pointsInRange.size() <<") :\n";
        for(int i=0;i<pointsInRange.size();i++)
        {
            outFile<<pointsInRange[i]->id <<"\n";
        }
        outFile<<"-----------------------------------------------------------\n";
        //clean up space and move to the next line in document
       delete[] temp->coords;
       delete temp;
    }
    //outFile<<"tHyperCube : "<<timeHyperCube<<"\n";
    //outFile<<"tTrueL "<<timeTrue<<"\n";
    outFile.close();
    printf("MAXIMUM ERROR distHyperCube/distTrue =%f\n",maxError);
    printf("AVERAGE ERROR distHyperCube/distTrue =%f\n",avgError/numErrors);
    printf("DIST ERRORS OVER 2.0: %d/%d = %f \n",numVeryBad,numQueries,(float)numVeryBad/(float)numQueries);
    printf("DIST ERRORS OVER 3.0: %d/%d = %f \n",numVeryVeryBad,numQueries,(float)numVeryVeryBad/(float)numQueries);

}

int main(int argc,char **argv)
{
    //srand(33);
    //default values , will later deal with command line arguments
    int k=14,M=10,probes=2,N=1;
    int d=-1;
    float R=10000;
    std::string in_fname="",query_fname="",out_fname="";
    std::vector<point *>points;
    HyperCube *hq;
    if(argc==1)
    {
        std::cout<<"PLEASE GIVE ME THE PATH TO THE INPUT FILE : "<<"\n";
    
        std::cin>>in_fname;
        std::cout<<"\n FILENAME = "<<in_fname <<"\n";
        
        hq=new HyperCube(k,M,probes);
        points=fileToHyperCube(*hq,in_fname.c_str(),d);

        std::cout<<"PLEASE GIVE ME THE PATH TO THE QUERRY FILE \n";
        std::cin>>query_fname;
        std::cout <<"PLEASE GIVE ME THE PATH TO THE OUTPUT FILE\n";
        std::cin>>out_fname;
        //std::cout<<" in main  d= "<<d<<'\n';
        dealWithTestSet(*hq,query_fname.c_str(),out_fname.c_str(),points,d,k,N,R);
    }
    else
    {
         for(int i=0;i<argc;i++)
        {
            if(strcmp(argv[i],"-i")==0)
                in_fname=argv[i+1];
            if(strcmp(argv[i],"-o")==0)
                out_fname=argv[i+1];
            if(strcmp(argv[i],"-q")==0)
                query_fname=argv[i+1];
            if(strcmp(argv[i],"-k")==0)
                k=atoi(argv[i+1]);   
            if(strcmp(argv[i],"-M")==0)
                M=atoi(argv[i+1]);
            if(strcmp(argv[i],"-probes")==0)
                probes=atoi(argv[i+1]);   
            if(strcmp(argv[i],"-R")==0)
                R=atoi(argv[i+1]);
            if(strcmp(argv[i],"-N")==0)
                N=atoi(argv[i+1]);            
            
        }
        if(R<=0 || k<=0 || M<=0 || N<=0 || probes<=0)
        {
            printf("INVALID ARGUMENTS\n");
            exit(-1);
        }
        hq=new HyperCube(k,M,probes);
        bool created=false;
        if(in_fname=="")
        {
            std::cout<<"PLEASE GIVE ME THE PATH TO THE INPUT FILE : "<<"\n";
    
            std::cin>>in_fname;

            points=fileToHyperCube(*hq,in_fname.c_str(),d);
            created=true;
        }
        if(query_fname=="")
        {
            std::cout<<"PLEASE GIVE ME THE PATH TO THE QUERRY FILE \n";
            std::cin>>query_fname;
        }
        if(out_fname=="")
        {
            std::cout <<"PLEASE GIVE ME THE PATH TO THE OUTPUT FILE\n";
            std::cin>>out_fname;
        }
        d=-1;
        if(created==false)
            points=fileToHyperCube(*hq,in_fname.c_str(),d);
        dealWithTestSet(*hq,query_fname.c_str(),out_fname.c_str(),points,d,k,N,R);

    }
    bool Continue=true;
    while(Continue)
    {

        printf("Would you like to continue for another query file (y/n) ?\n");
        std::string answer;
        std::cin>>answer;
        if(answer=="n")
            Continue=false;
        else if(answer=="y")
        {
            std::cout<<"PLEASE GIVE ME THE PATH TO THE QUERRY FILE \n";
            std::cin>>query_fname;
            std::cout<<"NEW QUERRY FILE : "<<query_fname<<"\n";
            dealWithTestSet(*hq,query_fname.c_str(),out_fname.c_str(),points,d,k,N,R);

        }
        else   
            printf("wrong answer please answer y or n \n");
       
    }

    
    

    delete hq;
    for(int i=0;i<points.size();i++)
    {
        delete[] points[i]->coords;
        delete points[i];
    }

      
}
