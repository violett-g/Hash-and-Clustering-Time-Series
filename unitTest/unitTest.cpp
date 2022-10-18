#include "../timeSeries/Search/1dCurves.hpp"
#include "../timeSeries/Search/2dCurves.hpp"

#include "../../../../CUnit-2.1-2/CUnit/Headers/Basic.h"
#include "../Common/Common.hpp"
#include "../timeSeries/Clustering/meanCurve.hpp"

#include <string>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cstring>
using namespace std;

//Documents/Project/prj2/TimeSeriesPRJ/Data/nasd_query.csv

/*utility function to test snaping and padding procedure*/
vector<point*> read_data(string path){
    vector<point*> data;
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
        data.push_back(curr_point); //push point to data structure
    }
    dataFile.close();
    return data;
}

/*function to test snaping and padding procces of a curve
snap_pad function performs snaping for all the curves of a file and than pads them 
for a given padding value = 730 and returns a vector  of points where the point stores the padded
curve = coords and the size = d 
we are going to test if the returned vectors that corespond to curves has a size same as the give padding size =730*/

void testSNAP_PAD(){
    cout << "GIVE FILE PATH in order to get some real curves to snap in a random grid and pad" << endl;
    string path ;
    cin >> path;
    vector<point*> queries = read_data(path);

    float delta = 3;
    float t = 1;
    float paddValue = 0;
    
    vector<point*> ret = snap_pad(queries,queries,delta,t,paddValue);
    // int returned_size = ret[0]->d; 
    int returned_size = ret[0]->d;
    CU_ASSERT(730 == returned_size);
}


/*utility functions to test mean curves*/
vector<discrete2DCurve*> fileTo2D(string path)
{
    std::ifstream infile(path);
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
    }
    return curves;
}

void testMEAN_CURVE2(){
    cout << "GIVE FILE PATH in order to get some real curve and calculate the mean" << endl;
    string path ;
    cin >> path;
    vector<discrete2DCurve*> curves = fileTo2D(path);

    int d = 10; //we check only the 10 first curves
    discrete2DCurve* ret_mean = meanOfTwoCurves(curves[0],curves[1]); //test the first two curves of the file
    
    
    float calc_mean[] ={34.965, 34.515, 33.7 , 34.69 , 69.63 , 34.445 , 34.055 , 34.005 , 34.875 , 334.1 };
    float ret_arr[10];
    for(int i=0; i<10; i++){
        ret_arr[i] = ret_mean->vertices[i].coords[1]; //y from curve verticies
    }
    float calc_mean_x[] ={1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
    float ret_arr_x[10];
    for(int i=0; i<10; i++){
        ret_arr_x[i] = ret_mean->vertices[i].coords[0]; //x from curve verticies
    }
    CU_ASSERT_EQUAL(0, memcmp(calc_mean, ret_arr, 2));
}

/*function to test getFrechet for discrete distance
we create two small dummy curves and calculated the Frechet distance 
and compared it with the one that our function returned*/
void testFRECHET(){
    //create two dummy curves
    int curve_d = 3;
    int point_d = 2;
    discrete2DCurve* c1 = new discrete2DCurve;
    discrete2DCurve* c2 = new discrete2DCurve;
   
    c1->vertices = new point[curve_d];
    c2->vertices = new point[curve_d];

    c1->d = curve_d;
    c2->d = curve_d;

    int A[3][2] = {{4,4}, { 2 , 6 } , {8,2}};
    int B[3][2] = {{6,8}, { 2 , 2 } , {4,2}};
    // int A[3][2] = {{1,4}, { 2 , 6 } , {5,5}};
    // int B[3][2] = {{1,3}, { 2 , 2 } , {4,3}};
    for(int i=0; i<curve_d; i++){
        
        point* p1 = new point;
        p1->d = point_d;
        p1->coords = new float[point_d];
        p1->coords[0]=A[i][0];
        p1->coords[1]=A[i][1];
        c1->vertices[i] = *p1;

        point* p2 = new point;
        p2->d = point_d;
        p2->coords = new float[point_d];
        p2->coords[0]=B[i][0];
        p2->coords[1]=B[i][1];
        c2->vertices[i] = *p2;
    } 
    
    //calling getFrechet
    // float exp = 4.47214;
    // float ret = getFrechet(*c1,*c2);
    // CU_ASSERT(0 == memcmp(&exp, &ret, sizeof(float)));

    int ret = getFrechet(*c1,*c2);
    int exp = 4;
    CU_ASSERT(exp == ret );
}


int main(){
    CU_pSuite pSuite = NULL;

    /* initialize the CUnit test registry */
    if (CUE_SUCCESS != CU_initialize_registry())
        return CU_get_error();

    /* add a suite to the registry */
    pSuite = CU_add_suite("Suite_1", NULL, NULL);
    if (NULL == pSuite) {
        CU_cleanup_registry();
        return CU_get_error();
    }

    /* add the tests to the suite */
    if ((NULL == CU_add_test(pSuite, "test of snaping and padding process", testSNAP_PAD))
        || (NULL == CU_add_test(pSuite, "test of mean curve for two curves", testMEAN_CURVE2))
        || (NULL == CU_add_test(pSuite, "test of discrete frechet distance ", testFRECHET))
        )
    {
        CU_cleanup_registry();
        return CU_get_error();
    }

    /* Run all tests using the CUnit Basic interface */
    CU_basic_set_mode(CU_BRM_VERBOSE);
    CU_basic_run_tests();
    CU_cleanup_registry();
    return CU_get_error();
}