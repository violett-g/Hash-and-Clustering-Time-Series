#pragma once
#include "Lloyd.hpp"
#include "../../NN/LSH/LSH.hpp"
#include "../../NN/Proj/Proj.hpp"
#include "../../Common/Common.hpp"
#include "../Search/2dCurves.hpp"
#include <string>
#include <math.h>
#include <iostream>
#include <cstring>
#include <vector>
#include <sstream>
#include <fstream>
#include <ctime>
//used for calculateing mean curve of n curves
struct treeNode{
    discrete2DCurve *curve;
    bool isLeaf;
    treeNode* left;
    treeNode* right;
    treeNode* parent;
};
/*function to initialize an empty tree structure
tree structure is an array where each node is marked as leaf or inner node and has pointer to its left and right child
leafs have NULL pointers for childs
then randomly scatter the curves to the leaves of the tree  */
treeNode *initialize_tree(std::vector<discrete2DCurve *> curves);
discrete2DCurve *meanOfNcurves(treeNode root);
/* function that takes two curves as arguments and return the mean curve */
discrete2DCurve *meanOfTwoCurves(discrete2DCurve *curve1, discrete2DCurve *curve2);
void deallocateTree(treeNode* tree,int num_of_curves);

