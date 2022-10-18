#pragma once
#include <vector>
#include <string>
#include "../../Common/Common.hpp"
#include "../../Fred/include/frechet.hpp"


using namespace std;

//filters a vector of points(their curves) and returns
//remember each point is associated with the original curve through point->curve pointer
vector<point*> filter_data(vector<point*> raw_data,float e);
//same as before, only for snaping and padding
vector<point*> snap_pad(vector<point*> filtered_curves, vector<point*>raw_data,float delta, float t,float paddValue);
Curve pointToCurve(point* point);

