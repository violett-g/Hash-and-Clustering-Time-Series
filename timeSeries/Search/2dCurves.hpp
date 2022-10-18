#pragma once
#include "../../Common/Common.hpp"
#include <vector>
//initialise the vector of ts of the grid with uniform distr
grid2D createGrid(float delta,int d);
//snapping,duplicate removal and padding
float * TS_to_Vector(discrete2DCurve curve,grid2D G,int vectorSize,float paddValue);
//discrete frechet between 2dimensional curves through dynamic programmming
float getFrechet(discrete2DCurve A,discrete2DCurve B);
// frechet distance of associated curves
float getFrechet(point a,point b);
//usefull for input handling converts a point to a curve , 
//points coords becoem the curves y coords, the x coords
//of the curve are 1,2,3, .... ,curve size
discrete2DCurve pointTo2dCurve(point *a);


