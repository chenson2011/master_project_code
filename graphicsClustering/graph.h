#ifndef GRAPH_H
#define GRAPH_H

#include <Eigen/Dense>
#include <vector>
using namespace Eigen;

class Graph
{
public:
    Graph(MatrixXd mat);
    MatrixXd mcl(int e,int r);
    MatrixXd spectral(MatrixXd mat);
    MatrixXd getAdj();


private:
    MatrixXd adj;
    MatrixXd inflate(MatrixXd grap,int r);
    MatrixXd expand(MatrixXd inf, int e);
    double error(MatrixXd ex, MatrixXd inf);
    vector classify(MatrixXd mat);
};

#endif // GRAPH_H
