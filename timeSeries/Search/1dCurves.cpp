#pragma once
#include "1dCurves.hpp"
#include "../../Common/input.hpp"
#include "../../Common/Common.hpp" //only for testing lsh
#include "../../NN/LSH/LSH.hpp"
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
#include <ctime>
using namespace std;



/*this function filters data from input. If we consider a curve with point [a,b,c,d,e] then,
for each three points a,b,c on the curve tests if |a-b|<= e and |b-c|<=e removes b, 
(that means add point a  on the filtered curve and test  c,d,e and test b,c,d)
otherwise keeps all three points (that means add point a  on the filtered curve)
and finally returns a vector of vectors with the filtered curves. */
vector<point*> filter_data(vector<point*> raw_data, float e){
    vector<point*> filtered_curves; //e vector of curves with filtered coordinates
    for(int data_point=0; data_point<raw_data.size(); data_point++ ){ //iterate curves on raw data
        point *curr_point=new point(); //create a new curve to store filtered curve
        float* curr_point_coords = raw_data[data_point]->curve; //get coordinates
        curr_point->id = raw_data[data_point]->id;
        vector<float> filtered_coords; //temporary vector to store filtered coordinates
        //for each three coordinates on the current raw curve test the condition
       
        for(int a=0; a<raw_data[data_point]->curve_d-2; a++){
            int b = a+1;
            int c = a+2;
            if(abs(curr_point_coords[a]-curr_point_coords[b])<=e && abs(curr_point_coords[b]-curr_point_coords[c]) <=e ){
                filtered_coords.push_back(curr_point_coords[a]);
                a++;
            }
            else{
                filtered_coords.push_back(curr_point_coords[a]);
            }
        }
        curr_point->coords=new float[filtered_coords.size()];
        for(int i=0; i<filtered_coords.size(); i++){//copy filtered coordinates of the curve to the curve to be returned by the function
            curr_point->coords[i]=filtered_coords[i];
        }
        curr_point->d=filtered_coords.size(); //store dimension
        filtered_curves.push_back(curr_point);//push curve with filtered coordinates
    }
    return filtered_curves;
}

/* function that snap curves on a grid and than perform pading, new coordinates are stored on a point struct
which stores the default curve and the new one after filtering/snaping/padding
function is written in such way to snap in different G_delta_t grids, if you want to snap a curve in a G_delta
grid pass parameter t=0*/

vector<point*> snap_pad(vector<point*> filtered_curves, vector<point*>raw_data,float delta,float t,float paddValue){
    vector<point*> snaped_curves;
    int max_d = 0;
    for(int i=0; i<filtered_curves.size(); i++){
        point *curr_point=new point(); //create a new curve to store snaped curve
        curr_point->id = raw_data[i]->id;

        float* curr_point_coords = raw_data[i]->coords; //get coordinates
        vector<float> snaped_coords; //temporary vector to store snaped coordinates
        //snaping
        for(int x=0;x<filtered_curves[i]->d; x++){
            float snaped_x = floor((filtered_curves[i]->coords[x]+0.5)/delta)*delta; //snaping in G_delta
            snaped_x = snaped_x+t; //snaping in G_delta_t
            //check for dublicates
            if (!count(snaped_coords.begin(), snaped_coords.end(), snaped_x)) {
                snaped_coords.push_back(snaped_x); //push snaped cordinate
            }
        }
        curr_point->coords=new float[snaped_coords.size()];
        for(int i=0; i<snaped_coords.size(); i++){//copy filtered coordinates of the curve to the curve to be returned by the function
            //minima-maxima----------->copy only after check
            curr_point->coords[i]=snaped_coords[i];
        }
        if(max_d < snaped_coords.size()) max_d = snaped_coords.size(); //keep the longest curve dimensions
        curr_point->d = snaped_coords.size();
        // curr_point->coords = coords;
        snaped_curves.push_back(curr_point);
    }

    //padding
    for(int i=0; i<snaped_curves.size(); i++){
        //cout << "PADDING" << endl;
        vector<float> pad_coords;
        //float *pad_coords=new float[max_d];
        // vector<float> pad;
        for(int j=0; j<snaped_curves[i]->d; j++){ //copy snaped coords
            //pad_coords[j] = snaped_curves[i]->coords[j];
            pad_coords.push_back(snaped_curves[i]->coords[j]);
        }
        for(int j=snaped_curves[i]->d; j<730; j++){
            //pad_coords[j]=0; //pad 0 in the end
            // pad_coords.push_back(0);
            pad_coords.push_back(paddValue);
        }
        // raw_data[i]->coords=new float[pad_coords.size()];
        for(int c=0; c<pad_coords.size(); c++){
            raw_data[i]->coords[c]=pad_coords[c];
        }
        // raw_data[i]->coords = new_coords;
        raw_data[i]->d = pad_coords.size();
    }
    //delete snaped 
    for (int i = 0; i < filtered_curves.size(); i++){
        delete[] snaped_curves[i]->coords;
        // delete[] snaped_curves[i]->curve;
        delete snaped_curves[i];
    }
    return raw_data;
}

/*function to convert our point struct that correspond to a curve in the files to a class
Curve of Fred libray*/

Curve pointToCurve(point* my_point){
    // cout << my_point->id << endl;
    // for(int i=0; i<my_point->curve_d; i++){
    //     cout << my_point->curve[i] << " , ";
    // }

    vector<Point> curve_points; // vector to stoure all the points of the curve
    for(int p=0; p < my_point->curve_d; p++){
        Coordinates coords;
        coords.push_back(my_point->curve[p]); //save coordinate in coords-we have only one coordinate because we have an 1 dimensional curve
        Point* point = new Point(1); //we have 1 domensional curve points
        point->set(0,coords[0]); // set the only coordinate value in the point coordinate vector
        curve_points.push_back(*point); //store points on vector to copy later on the Points class
    }

    //create an instance of Points class, needed to create the Curve
    const curve_size_t size = my_point->curve_d;
    Points* points = new Points(1);
    for(int p=0; p<curve_points.size(); p++){
        // cout << "GETING COORDINATES FROM POINT " << curve_points[p]->get(0) << endl;
        points->add(curve_points[p]); //add points 
    }
    //create an instance Curve 
    Curve* curve = new Curve(*points, my_point->id);
    // cout << my_point->id <<endl;
    // for(int i=0; i<my_point->curve_d; i++){
    //     cout << curve->get(i).get(0) << ", ";
    // }
    // cout << endl;
    return *curve;
}



float continousFrechet(point a,point b)
{
    Curve query_c = pointToCurve(&a);
    Curve data_c = pointToCurve(&b);
    return  Frechet::Continuous::distance(query_c,data_c).value;
}
