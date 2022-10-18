#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <cstring>
#include <limits.h>
//#include "../../Common/Common.hpp"
#include "../../NN/LSH/LSH.hpp"
#include "../../NN/Proj/Proj.hpp"
#include "../../Common/input.hpp"
#include "1dCurves.hpp"
#include <random>

#include <ctime>
#include "2dCurves.hpp"


//A i) LSH part
void dealWithTestSet
(
    LSH &lsh,const char *testFName,const char *outFName,
    std::vector<point *> &data,
    int d=2,int k=4,int L=5,int N=1
);
//A i) Hypercube part
void dealWithTestSet(
    HyperCube &hq,const char *testFName,const char *outFName,
    std::vector<point *> &data,
    int d=2,int k=4,int L=5,int N=1
    );
// A ii)
void dealWithTestSet
(
    LSH **lshs,grid2D *grids,const char *testFName,const char *outFName,
    std::vector<discrete2DCurve *>data,int L,int d,int vectorSize,float paddValue
);
///test///
//A iii)
void dealWithTestSet
(LSH **lsh,std::vector<point*>data, std::vector<std::vector<point*>>queries_vector,std::vector<point*>queries,
  std::string outFName, int L, float duration_preprocces
);
//this generator is used in all the experiments
extern std::default_random_engine generator;

int main(int argc,char **argv)
{
    //seed it in main for reproducability
    generator.seed(185);
    srand(185);
    //read argc,argv
    Input myInput(argc,argv);
 


    LSH *lsh;
    HyperCube *hq;
    int d=-1;

    std::uniform_real_distribution<float> uni_distribution(1000000,INT_MAX);
    float paddValue=uni_distribution(generator);
        
            
    //first case, lsh, time series are viewed as vectors
    if(myInput.algorithm=="LSH")
    {
        lsh =new LSH(d,50,myInput.k,myInput.L,dist);
        myInput.data=fileToLSH(*lsh,myInput.inputPath.c_str(),d);

        dealWithTestSet(*lsh,myInput.queryPath.c_str(),myInput.outputPath.c_str(),myInput.data,d);
        bool cont=true;
        while(cont)
        {
            std::string answer;
            printf("do you want to try for an other query file? (y/n)\n");
            std::cin>>answer;
            if(answer=="y")
            {
                printf("Give me the path to the new query file\n");
                std::cin>>myInput.queryPath;
                dealWithTestSet(*lsh,myInput.queryPath.c_str(),myInput.outputPath.c_str(),myInput.data,d);


            }
            else {
                cont=false;
            }

        }
        delete lsh;
        for(int i=0;i<myInput.data.size();i++)
        {
            delete[] myInput.data[i]->coords;
            delete myInput.data[i];
        }
    }
    //second case, hypercube time series are viewed as vectors
    if(myInput.algorithm=="Hypercube")
    {
            
        hq=new HyperCube(myInput.k,myInput.M,myInput.probes);
        myInput.data=fileToHyperCube(*hq,myInput.inputPath.c_str(),d);
    
        dealWithTestSet(*hq,myInput.queryPath.c_str(),myInput.outputPath.c_str(),myInput.data,d);
        bool cont=true;
        while(cont)
        {
            std::string answer;
            printf("do you want to try for an other query file? (y/n)\n");
            std::cin>>answer;
            if(answer=="y")
            {
                printf("Give me the path to the new query file\n");
                std::cin>>myInput.queryPath;
                dealWithTestSet(*hq,myInput.queryPath.c_str(),myInput.outputPath.c_str(),myInput.data,d);


            }
            else {
                cont=false;
            }

        }
        delete hq;
        for(int i=0;i<myInput.data.size();i++)
            {
                delete[] myInput.data[i]->coords;
                delete myInput.data[i];
            }
    }
    //FRECHED DISTANCE
    else if(myInput.algorithm=="Frechet")
    {
        //discrete case , 2d cruves
        if(myInput.metric=="discrete")
        {
            grid2D *grids=new grid2D[myInput.L];
            int vectorSize=1500;
            LSH **lsh=new LSH*[myInput.L];
            //THERE ARE L LSH TABLES ONE FOR EACH GRID
                
            //OPEN THE INPUT FILE
            std::ifstream infile(myInput.inputPath);
            if(infile.fail())
            {
                printf("ERROR: FAIL IN INPUT FILE\n");
                exit(-1);
            }
            std:: string line;
            int d=-1;
            std::vector<discrete2DCurve *>curves;
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
                    
                curves.push_back(cv);
                delete[] temp;

            }
            //CREATE DATA STRUCTURES
            for(int i=0;i<myInput.L;i++)
            {
                lsh[i]=new LSH(vectorSize, 150,myInput.k,1,getFrechet);
                grids[i]=createGrid(myInput.delta,d);
            }
            std::unordered_map<int,std::vector<point *> > pointsPerGrid;
            //FOR EACH CURVE
            for(int i=0;i<curves.size();i++)
            {   
                //FOR EACH OF THE L GRIDS
                for(int j=0;j<myInput.L;j++)
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
            //store the created points to the appropriate LSH tables
            for(int i=0;i<myInput.L;i++)
            {
                lsh[i]->fit(pointsPerGrid[i].size(),pointsPerGrid[i].data(),vectorSize,pointsPerGrid[i].size()/2);
            }
            //now read the test file and output
            dealWithTestSet(
                lsh,grids,myInput.queryPath.c_str(),myInput.outputPath.c_str(),
                curves,myInput.L,d,vectorSize,paddValue); 
            bool cont=true;
            while(cont)
            {
                std::string answer;
                printf("do you want to try for an other query file? (y/n)\n");
                std::cin>>answer;
                if(answer=="y")
                {
                    printf("Give me the path to the new query file\n");
                    std::cin>>myInput.queryPath;
                    dealWithTestSet(
                        lsh,grids,myInput.queryPath.c_str(),myInput.outputPath.c_str(),
                        curves,myInput.L,d,vectorSize,paddValue); 

                }
                else {
                    cont=false;
                }

            }
            //CLEAN UP SPACE
            delete[] grids;
            for(int i=0;i<myInput.L;i++)
            {
                delete lsh[i];
                for(int j=0;j<pointsPerGrid[i].size();j++)
                {
                    delete[] pointsPerGrid[i][j]->coords;
                    delete pointsPerGrid[i][j];
                }
            }
            delete[] lsh;
            for(int i=0;i<curves.size();i++)
            {
                for(int j=0;j<curves[i]->d;j++)
                {
                    delete[] curves[i]->vertices[j].coords;
                }
                delete[] curves[i]->vertices;
                delete curves[i];
            }
        }
        //QUESTION AIII CONTINOUS FRECHET
        else
        {
            //std::cout<<"CONTINOUS"<<std::endl;
                
            //read data and queries
            vector<point*> data;
            vector<point*> queries;
            data=myInput.read_data(data,myInput.inputPath );
            queries= myInput.read_data(queries,myInput.queryPath );

            float e = 0.01;//used on filtering
            int w=40;

                //filter data and queries
            vector<point*> filtered_data_curves = filter_data(data,e); //filter data curves
                
            clock_t start;
            start = std::clock();//calculate time nedded for filtering queries
            vector<point*> filtered_query_curves = filter_data(queries,e); //filter queries curves
            float duration_preprocces = ( std::clock() - start ) / (float) CLOCKS_PER_SEC;

            //create L different tables and L different grids, snap and pad data and queries in grids
            std::default_random_engine generator;
            std::uniform_real_distribution<> distribution(0, myInput.delta-1);
            vector<vector<point*>> queries_vector; //vector of vectors with snaped data curves in vectors
            vector<vector<point*>> data_vector;//vector of vectors with snaped queries curves in vectors
                // vector<float> ts;
            LSH **lsh=new LSH*[myInput.L];
            for(int l=0; l<myInput.L; l++){
                lsh[l] = new LSH(data.size(),w,myInput.k, 1,dist);
                float t = distribution(generator); //generate L diferent grids
                    // ts.push_back(t);
                vector<point*> snap_data = snap_pad(filtered_data_curves,data,myInput.delta,t, paddValue); //snap data curve in the grid
                data_vector.push_back(snap_data);

                lsh[l]->fit(snap_data.size(),snap_data.data(), data[0]->d, snap_data.size()/16 );//store in lsh tables

                start = std::clock();
                vector<point*> snap_q = snap_pad(filtered_query_curves,queries,myInput.delta,t, paddValue); //snap queries curve in the grid
                duration_preprocces +=( std::clock() - start ) / (float) CLOCKS_PER_SEC;
                queries_vector.push_back(snap_q);
            }
                // for(int i=0; i<myInput.L; i++){
                //     lsh[i]->printTables();
                // }
            dealWithTestSet(lsh,data,queries_vector,queries, myInput.outputPath,myInput.L,duration_preprocces);
            
            while(true)
            {
                std::string answer;
                printf("do you want to try for an other query file? (y/n)\n");
                std::cin>>answer;
                if(answer=="y")
                {
                    printf("Give me the path to the new query file\n");
                    std::cin>>myInput.queryPath;
                    queries= myInput.read_data(queries,myInput.queryPath );

                    dealWithTestSet(lsh,data,queries_vector,queries, myInput.outputPath,myInput.L,duration_preprocces);
                    for (int i = 0; i < queries.size(); i++)
                    {
                        delete[] queries[i]->coords;
                        delete[] queries[i]->curve;
                        delete queries[i];
                    }   

                }
                else {
                    break;
                }

            }
            //clean up space
            for(int i=0;i<filtered_data_curves.size(); i++){
                delete[] filtered_data_curves[i]->coords;
                delete filtered_data_curves[i];
            }
            for(int i=0;i<filtered_query_curves.size(); i++){
                delete[] filtered_query_curves[i]->coords;
                delete filtered_query_curves[i];
            }
            for (int i = 0; i < data.size(); i++)
            {
                delete[] data[i]->coords;
                delete[] data[i]->curve;
                delete data[i];
            }
            // for (int i = 0; i < queries.size(); i++)
            // {
            //     delete[] queries[i]->coords;
            //     delete[] queries[i]->curve;
            //     delete queries[i];
            // }
            for(int i=0; i<myInput.L; i++){
                delete lsh[i];
            }
            delete[] lsh;


        }
    }
    /*string more;
    string new_queryPath;
    cout << "Would you like to proccess a new queryset? (yes/no)" << endl;
        cin >> more;
        if(more == "yes"){
            cout << "Type the new file path for the new queryset" << endl;
            cin >> new_queryPath;
            myInput.queryPath = new_queryPath;
        }
        else {
            break;
        }*/
    

}
