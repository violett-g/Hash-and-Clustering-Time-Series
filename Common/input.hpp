#ifndef INPUT_HPP
#define INPUT_HPP
#include <vector>
#include <string>
#include "Common.hpp"


using namespace std;

class Input{
    public:
        //input parameters
        int k; //number of hi functions
        int k_HQ; //k of hypercube
        int L; //number of hashtables
        int N; //number of closest neighbours
        int M; //max number of data to be checked
        int probes; //max number of edges on the hypercube
        string algorithm; //
        string metric; //
        float delta; // 

        string queryPath; 
        string outputPath;
        string inputPath;

        //structures
        vector<point*> data; //structure to store datas
        vector<point*> queries; //structure to store query points
        
        //functions
        void parse_arguments(int argc, char* argv[]);
        vector<point*> read_data(vector<point*> data,string path); // function to read input data from file create item and store to heap
        // void find_neighbors(Distance_Function dist_func);
        //void read_data(); // function to read input data from file create item and store to heap
        // void read_query(); // function to read queries create query item and store to  heap
        //vector<point*> filter_data(vector<point*> raw_data, int e);
        //vector<point*> snap_pad(vector<point*> filtered_curves, vector<point*>raw_data,int delta, float t);
        //void test(int e,int L,int delta,int w,int k);

    
    Input(int argc, char* argv[]);
    ~Input();
    void more_queries(string new_queryPath);

};

#endif


