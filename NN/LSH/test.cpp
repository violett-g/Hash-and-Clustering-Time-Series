#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <cstring>
#include "LSH.hpp"
#include "../../Common/Common.hpp"
#include <ctime>






//this function will parse the test(query) file line by line 
//for each vector (point) it will output the k-nn and range search results
//according to the required format. If there is something wrong with
//the query file or the output file it will exit printing the appropriate message
void dealWithTestSet(
    LSH &lsh,const char *testFName,const char *outFName,std::vector<point *> &data,
    int d=2,int k=4,int L=5,int N=1,float R=100000.0 
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
    //the files are all open go through the query file line by line 
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);

        if(line.length() <2)
            break;
        //extract point fromthe test set
        point *temp=new point();
        temp->d=d;
        temp->coords=strToCoords(line.data(),d,temp->id);
       
        std::clock_t start;
        float durationLSH,durationTrue;
        //first LSH querry for kNN
        start = std::clock();
        std::vector<std::pair<point*, float>> neighbors=lsh.approxKNN(*temp,N);
        durationLSH = ( std::clock() - start ) / (float) CLOCKS_PER_SEC;

        //then brute force query for kNN
        start = std::clock();
        std::vector<std::pair<point*, float>> bruteForceNeighbors=findNNBruteForce(*temp,data,N);
        durationTrue = ( std::clock() - start ) / (float) CLOCKS_PER_SEC;
        /* calculate the current error and see if its worse than max error*/
        if(neighbors.size()>0)
        {
            float curDist=neighbors[0].second/dist(*(bruteForceNeighbors[0].first),*temp);
            if(curDist > maxError)
                maxError=curDist;
            avgError+=curDist;
            numErrors+=1;
        }
        //now we got the k best neighbors (both LSH and brute force) , print then on outfile
        outFile<<"Querry :"<<temp->id <<"\n";
        for(int i=0;i<neighbors.size();i++)
        {  
            outFile<<"Nearest Neighbor-"<<i+1<<" : "<<neighbors[i].first->id<<"\n";
            outFile<<"distance LSH: "<<neighbors[i].second<<"\n";
            outFile<<"distance True: "<<dist(*(bruteForceNeighbors[i].first),*temp)<<"\n";
        }
        outFile<<"tLSH: "<<durationLSH<<"\n";
        outFile<<"tTrue: "<<durationTrue<<"\n";
        outFile<<"\n R-near neighbors:\n";
        std::vector<point *>pointsInRange=lsh.approxRangeSearch(*temp,R);
        for(int i=0;i<pointsInRange.size();i++)
        {
            outFile<<pointsInRange[i]->id <<"\n";
        }
        //this line is to make the output file appear cleaner
        outFile<<"-----------------------------------------------------------\n";
        //clean up space and move to the next line in document
       delete[] temp->coords;
       delete temp;
    }
    //outFile<<"tLSH : "<<timeLSH<<"\n";
    //outFile<<"tTrueL "<<timeTrue<<"\n";
    outFile.close();
    printf("MAXIMUM ERROR distLSH/distTrue =%f\n",maxError);
    printf("AVERAGE ERROR distLSH/distTrue =%f\n",avgError/numErrors);

}

int main(int argc,char **argv)
{
    //default values , will later deal with command line arguments
    int k=4, N=1,L=5;
    float R=10000.0;
    int d=-1;
    std::string fname="",outFname="",testFname="";
    std::vector<point *>points;
    LSH *lsh;
    //read file name
    if(argc==1)
    {
        std::cout<<"PLEASE GIVE ME THE PATH TO THE INPUT FILE : "<<"\n";
    
        std::cin>>fname;
        std::cout<<"\n FILENAME = "<<fname <<"\n";
        
        lsh =new LSH(d,300,k,L);

        points=fileToLSH(*lsh,fname.c_str(),d);

        std::cout<<"PLEASE GIVE ME THE PATH TO THE QUERRY FILE \n";
        std::cin>>testFname;
        std::cout <<"PLEASE GIVE ME THE PATH TO THE OUTPUT FILE\n";
        std::cin>>outFname;
        std::cout<<" in main  d= "<<d<<'\n';
        dealWithTestSet(*lsh,testFname.c_str(),outFname.c_str(),points,d);
    }
    else
    {
        //if a parameter is missing  it is filled with the default value
        for(int i=0;i<argc;i++)
        {
            if(strcmp(argv[i],"-i")==0)
                fname=argv[i+1];
            if(strcmp(argv[i],"-o")==0)
                outFname=argv[i+1];
            if(strcmp(argv[i],"-q")==0)
                testFname=argv[i+1];
            if(strcmp(argv[i],"-k")==0)
                k=atoi(argv[i+1]);   
            if(strcmp(argv[i],"-L")==0)
                L=atoi(argv[i+1]);
            if(strcmp(argv[i],"-N")==0)
                N=atoi(argv[i+1]);   
            if(strcmp(argv[i],"-R")==0)
                R=atoi(argv[i+1]);         
            
        }
        if(k<=0 || L<=0 || N<=0 || R<=0)
        {
            printf("INVALID COMMAND LINE ARGUMENTS \n");
            exit(1);
        }
        d=-1;
        lsh =new LSH(d,300,k,L);
        bool created=false;
        if(fname=="")
        {  
            std::cout<<"PLEASE GIVE ME THE PATH TO THE INPUT FILE \n";
            std::cin>>fname;
            points=fileToLSH(*lsh,fname.c_str(),d);
            created=true;

        }
        if(testFname=="")
        {
            std::cout<<"PLEASE GIVE ME THE PATH TO THE QUERRY FILE \n";
            std::cin>>testFname;
        }
        if(outFname=="")
        {
            std::cout <<"PLEASE GIVE ME THE PATH TO THE OUTPUT FILE\n";
            std::cin>>outFname;
        }
        if(created==false)
            points=fileToLSH(*lsh,fname.c_str(),d);
        dealWithTestSet(*lsh,testFname.c_str(),outFname.c_str(),points,d,k,L,N,R);

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
            std::cin>>testFname;
            std::cout<<"NEW QUERRY FILE : "<<testFname<<"\n";
            //points=fileToLSH(lsh,fname.c_str(),d);
            //dealWithTestSet(lsh,testFname.c_str(),outFname.c_str(),points,d,k,L,N,R);
            dealWithTestSet(*lsh,testFname.c_str(),outFname.c_str(),points,d,k,L,N,R);
        }
        else 
            printf("wrong answer ! \n");
       
    }
    delete lsh;
    for(int i=0;i<points.size();i++)
    {
        delete[] points[i]->coords;
        delete points[i];
    }

}
