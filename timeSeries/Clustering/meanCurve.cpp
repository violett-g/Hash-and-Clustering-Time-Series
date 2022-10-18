#include "Lloyd.hpp"
#include "../../NN/LSH/LSH.hpp"
#include "../../NN/Proj/Proj.hpp"
#include "../../Common/Common.hpp"
#include <string>
#include <math.h>
#include <iostream>
#include <cstring>
#include <vector>
#include <sstream>
#include <fstream>
#include <ctime>
#include "meanCurve.hpp"
/*function to initialize an empty tree structure
tree structure is an array where each node is marked as leaf or inner node and has pointer to its left and right child
leafs have NULL pointers for childs
then randomly scatter the curves to the leaves of the tree  */
treeNode *initialize_tree(std::vector<discrete2DCurve *> curves)
{
    int num_of_curves = curves.size();
    int h = ceil(log2(num_of_curves));   //calculate tree height
    int size = pow(2,h+1) - 1;        //calculate tree structure size
    treeNode *tree = new treeNode[size]; //initialize tree structure = array of nodes
    for (int i = 0; i < size; i++)
    {
        if (2 * i + 1 < size && 2 * i + 2 < size)
        {
            tree[i].left = &tree[2 * i + 1];  //set right child f node
            tree[i].right = &tree[2 * i + 2]; //set left child of node
            tree[i].isLeaf = 0;               //set node an inner node
            tree[i].curve = NULL;
        }
        else
        {
            tree[i].right = NULL;
            tree[i].left = NULL;
            tree[i].curve = NULL;
            tree[i].isLeaf = 1; //set node as a leaf node
        }
    }
    //randomly scatter the curves to the leaves of the tree
    int j = 0;
    for (int i = 0; i < size; i++)
    { //traverse tree structure
        if (tree[i].isLeaf == 1){ //if a leaf node is found
            
            if (j < curves.size()){

                tree[i].curve=curves[j];           
            
                j++;
            }
        }
    }
    /*for(int i=0; i<size; i++){
        if (tree[i].isLeaf == 1){
        if(tree[i].curve == NULL) std::cout  << i << tree[i].curve << std::endl;
        if(tree[i].curve != NULL) {
            std::cout << i << std::endl;
        std::cout << tree[i].curve->vertices[0].coords[0] << std::endl;
        }
        }
    }*/
    return tree;
}

/*a function that calculates the mean curve of n given curves
this function corresponds to the PostOrderTraversal function shown in the corse 
****should be called with root(tree[0]) as parameter */
discrete2DCurve *meanOfNcurves(treeNode root)
{
    if (root.isLeaf == 1)
    {
        if(root.curve == NULL) 
        {
            return NULL;
        } 
        discrete2DCurve *res=new discrete2DCurve();
        res->d=root.curve->d;
        res->vertices=new point[res->d];
        for(int i=0;i<res->d;i++)
        {
            res->vertices[i].coords=new float[2];
            res->vertices[i].coords[0]=root.curve->vertices[i].coords[0];
            res->vertices[i].coords[1]=root.curve->vertices[i].coords[1];

        }  
        return res;
    }
    else
    {
        discrete2DCurve *right;
        discrete2DCurve *left =  meanOfNcurves(*root.left);
        
        if (root.right != NULL)
        {
            right = meanOfNcurves(*root.right);
        }
        else
        {
            right = NULL;
        }
        discrete2DCurve *res=meanOfTwoCurves(left, right);
        if(left!=NULL)
        {
            for(int i=0;i<left->d;i++)
                delete[] left->vertices[i].coords;
            delete[] left->vertices;
            delete left;
        }
        if(right!=NULL)
        {
           for(int i=0;i<right->d;i++)
                delete[] right->vertices[i].coords;
            delete[] right->vertices;
            delete right; 
        }
        return res;
    }
}

/* function that takes two curves as arguments and return the mean curve */
discrete2DCurve *meanOfTwoCurves(discrete2DCurve *curve1, discrete2DCurve *curve2)
{
    if(curve1 == NULL && curve2==NULL){
        return NULL;
    }
   
    if (curve1 == NULL ){
        discrete2DCurve *mean_curve=new discrete2DCurve();
        mean_curve->d=curve2->d;
        mean_curve->vertices = new point[curve2->d];
        for (int i = 0; i < curve2->d; i++)
        { //for each point of curve
            mean_curve->vertices[i].coords = new float[2];
            mean_curve->vertices[i].d=2;
            mean_curve->vertices[i].coords[0]=curve2->vertices[i].coords[0];
            mean_curve->vertices[i].coords[1]=curve2->vertices[i].coords[1];

        }
        return mean_curve;
    }
    if (curve2 == NULL){
        discrete2DCurve *mean_curve=new discrete2DCurve();
        mean_curve->d=curve1->d;
        mean_curve->vertices = new point[curve1->d];
        for (int i = 0; i < curve1->d; i++)
        { //for each point of curve
            mean_curve->vertices[i].coords = new float[2];
            mean_curve->vertices[i].d=2;
            mean_curve->vertices[i].coords[0]=curve1->vertices[i].coords[0];
            mean_curve->vertices[i].coords[1]=curve1->vertices[i].coords[1];

        }
        return mean_curve;
    }
   
    discrete2DCurve *mean_curve = new discrete2DCurve();
    mean_curve -> d = curve1->d;
    mean_curve->vertices = new point[curve1->d];
    for (int i = 0; i < curve1->d; i++)
    { //for each point of curve
        mean_curve->vertices[i].coords = new float[2];
        mean_curve->vertices[i].d=2;
        mean_curve->vertices[i].coords[0] = (curve1->vertices[i].coords[0] + curve2->vertices[i].coords[0]) / 2;
        mean_curve->vertices[i].coords[1] = (curve1->vertices[i].coords[1] + curve2->vertices[i].coords[1]) / 2;
    }
    /*for (int i = 0; i < curve1->d; i++)
    {
        std::cout << "x:" << mean_curve->vertices[i].coords[0] << " "
             << "y:" << mean_curve->vertices[i].coords[1];
    }*/
    ///std::cout << std::endl;
    return mean_curve;
}
/* h zvise to vector kai perase san orisma to int num_of_curves = curves.size() an exeis diagrapsei pio prin ta curves*/
void deallocateTree(treeNode* tree,int num_of_curves){
    int h = ceil(log2(num_of_curves));   //calculate tree height
    int size = pow(2,h+1) - 1; 

    for(int i=0; i< size; i++){
        if(tree[i].isLeaf)
        {
            continue; 
        }
        if(tree[i].curve !=NULL)
        {
            for(int d=0; d<tree[i].curve->d; d++ ){
                delete[] tree[i].curve->vertices[d].coords;
            }
            delete[] tree[i].curve->vertices; //vertices
            delete tree[i].curve; //discrete2dcurve
        }
        //delete &tree[i]; //treeNode
    }
    delete[] tree;
}