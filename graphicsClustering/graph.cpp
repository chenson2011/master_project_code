#include "graph.h"
#include <Eigen/Dense>

using namespace Eigen;

Graph::Graph(MatrixXd mat)
{
    adj = mat;
}

MatrixXd Graph::getAdj()
{
    return adj;
}

MatrixXd Graph::mcl(int e,int r)
{
        int size = adj.rows();
        MatrixXd inf(size,size);
        MatrixXd ex(size,size);
        adj = adj + MatrixXd::Identity(size,size);
        inf = inflate(adj,1);
        double err = 1;

        while(err>0.0001)
        {
            ex = expand(inf,e);m
            inf = inflate(ex,r);
            err = error(ex,inf);
          //  cout<<err<<endl;

        }

        return inf;
}

// inflate
MatrixXd Graph::inflate(MatrixXd grap,int r)
{
    int size = grap.rows();
    MatrixXd sumM(1,size);
    MatrixXd xPower(size,size);
    MatrixXd result(size,size);

    for(int col=0; col<size; col++)
        for(int row=0; row<size; row++)
        {
            xPower(row,col) = pow(grap(row,col),r);
        }

    for(int col=0; col<size; col++)
    {
        double sum = 0;
        for(int row=0; row<size; row++)
        {
            sum = sum + xPower(row,col);
        }
        sumM(0,col) = sum;
    }

    for(int col=0; col<size; col++)
        for(int row=0; row<size; row++)
        {
            result(row,col) = xPower(row,col)/sumM(0,col);
        }

    return result;

}

MatrixXd Graph::expand(MatrixXd inf, int e)
{
    MatrixXd ex = inf;
    while(e > 1)
    {
        ex = ex*inf;
        e--;
    }

    return ex;
}

double Graph::error(MatrixXd ex, MatrixXd inf)
{
    int size = ex.rows();
    double err;
    for(int row=0; row<size; row++)
        for(int col=0; col<size; col++)
        {
            err = err + pow(ex(row,col) - inf(row,col),2);
        }

    return err;
}

vector classify(MatrixXd::mat)
{

}
