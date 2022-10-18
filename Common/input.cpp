#include "input.hpp"
//#include "Common.hpp" //only for testing lsh
//#include "../NN/LSH/LSH.hpp"
#include <string>
#include <random>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cstring>
#include <ctype.h>
#include <cmath>
using namespace std;


Input::Input(int argc, char *argv[]){
    //default values
    k=4;
    L=5;
    delta=0.05;
    probes=2;
    M=10;
    k_HQ=14;
    //parse the arguments to see if some parameter is given by the user
    parse_arguments(argc,argv);
   
}

Input::~Input(){
  
}

void Input::parse_arguments(int argc, char *argv[]){
    for (int i = 0; i < argc; i++){
        if (strcmp(argv[i], "-i") == 0) {inputPath = argv[i + 1]; cout << inputPath << endl;}
        else if (strcmp(argv[i], "-q") == 0)    queryPath = argv[i + 1];
        else if (strcmp(argv[i], "-k") == 0) 
        {
            k = atoi(argv[i + 1]);
            k_HQ=k; //k_HQ is the same as k of lsh, only difference in default values
        }
        else if (strcmp(argv[i], "-L") == 0)    L = atoi(argv[i + 1]);
        else if (strcmp(argv[i], "-M") == 0)    M = atoi(argv[i + 1]);
        else if (strcmp(argv[i], "-probes") == 0)   probes = atoi(argv[i + 1]);
        else if (strcmp(argv[i], "-o") == 0)    outputPath = argv[i + 1];
        else if (strcmp(argv[i],"-algorithm") == 0) algorithm = argv[i+1];
        else if (strcmp(argv[i],"-metric") == 0) metric = argv[i+1];
        else if(strcmp(argv[i],"-delta") == 0) delta = atof(argv[i+1]);
    }
    
    //checks for missing obligatory parameters -> input file path and query file path
    if (inputPath.empty()) {
        cout << "Please provide dataset path" << endl;
        cin>>this->inputPath;
    }
    if (outputPath.empty()){
        cout << "Please provide output path" << endl;
        cin>>this->outputPath;
    }
    if (queryPath.empty()){
        cout << "Please provide query set path" << endl;
        cin>>this->queryPath;
    }
}
/*function to read data from csv file and convert and store to datastructures*/
vector<point*> Input::read_data(vector<point*> data,string path){
    //cout << path << endl;
    ifstream dataFile;
    string line;
    int d = -1;
    dataFile.open(path); //open data file

    if (!dataFile.is_open()){ //check if file is opened
        perror("Error open");
        exit(EXIT_FAILURE);
    }
    while (getline(dataFile, line)){
        istringstream iss(line);

        if(line.length() <2)
            break;
        //create point struct
        point *curr_point=new point();
        curr_point->curve=strToCoords(line.data(),d,curr_point->id);
        curr_point->curve_d=d;
        curr_point->coords=new float[730];//same as padding size
        data.push_back(curr_point); //push point to data structure
    }
    cout << endl;
    dataFile.close();
    //if something is wrong uncomment the following to debug
    // Displaying the 2D vector ----TESTING------
        // for (int i = 0; i < data.size(); i++) {
        //     float* temp_vec = data[i]->curve;
        //     for (int j = 0; j < data[i]->curve_d; j++)
        //         cout << temp_vec[j] << " ";
        //     cout << endl;
        // }
    return data;
}
