#ifndef GRAPH_H
#define GRAPH_H

#include <Eigen/Dense>
#include <vector>
#include "boost/graph/adjacency_matrix.hpp"
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/simple_point.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/circle_layout.hpp>
#include <boost/graph/fruchterman_reingold.hpp>
#include <boost/graph/kamada_kawai_spring_layout.hpp>
#include "kmean.h"
using namespace Eigen;
using namespace std;

struct vect{
    int index;
    complex<double> data;
};

class Graph
{
public:
    Graph(MatrixXd mat);
    MatrixXd mcl(int e,int r);
    int* spectral(int k);
    MatrixXd getAdj();

    MatrixXd inflate(MatrixXd grap,int r);
    MatrixXd expand(MatrixXd inf, int e);
    double error(MatrixXd ex, MatrixXd inf);
    vector<vector<double> > classify(MatrixXd mat);


private:
    MatrixXd adj;

};

#endif // GRAPH_H
