#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <list>
#include "kmean.h"
#include <float.h>

using namespace Eigen;
using namespace std;

VectorXd init(int groupNum, int pointNum)
{
    VectorXd randomV(groupNum);
    for(int i=0; i<groupNum; i++)
    {
        int value = rand()%pointNum;
        for(int j=0; j<=i; j++)
        {
            while(randomV(j)== value)
            {
                value = rand()%pointNum;
            }
        }
        randomV(i) = value;
    }

    return randomV;
}

/*VectorXd calCenter(vector<VectorXd> group)
{
    int pNum = group.size();
    int dim = group[0].size();
    VectorXd center(dim);

    double sum = 0.0;

    for(int i = 0; i<dim; i++)
    {
        sum = 0.0;
        for(int j=0; j<pNum; j++)
        {
            VectorXd temp = group[j];
            sum = sum + temp(i);
        }
        center(i) = sum/pNum;
    }

    return center;
}
*/

double calDistance(VectorXd v, VectorXd center)
{
    int size = v.rows();

    double sum = 0;
   for(int i=0; i<size; i++)
    {
        sum = sum + std::abs(v[i]-center[i]);
    }


   sum = pow(sum,1.0/2);

    return sum;
}

int* kmean(MatrixXd mat, int groupNum)
{
    int count = 100;
    int pointNum = mat.rows();
    int dim = mat.cols();
    //vector<VectorXd> *test[] = new vector<VectorXd>[groupNum];
    int *label = new int[pointNum];
    VectorXd initV = init(groupNum,pointNum);
    vector<VectorXd> center;

    for(int i=0; i<groupNum; i++)
    {
        center.push_back(mat.row(initV[i]).transpose());
    }

 /*   for(int i=0; i<groupNum; i++)
    {
        cout<<"init = "<<initV[i]<<endl;
        cout<<"group "<<i<<endl<<center[i]<<endl;
    }
    */


    while(count>0)
    {
        for(int i=0; i<pointNum; i++)
        {
            double min = DBL_MAX;
            VectorXd point = mat.row(i).transpose();
            int index;
            for(int j=0; j<groupNum; j++)
            {
                VectorXd C = center[j];
                double dist = calDistance(point,C);
               // cout<<"min = " <<min<<endl;
              //  cout<<"dist = "<<dist<<endl;
                if(dist<min)
                {
                    index = j;
                    min = dist;
                }
            }

            label[i] = index;
        }

      /*  for(int i=0; i<pointNum; i++)
        {
            cout<<label[i]<<endl;
        }
        */
        center.clear();
        for(int i=0; i<groupNum; i++)
        {
            VectorXd sum = VectorXd::Zero(dim);

            int flag = 0;
            for(int j=0; j<pointNum; j++)
            {
                if(label[j] == i)
                {
                  //  cout<<mat.row(j).transpose()<<endl;
                    sum = sum + mat.row(j).transpose();
                //    cout<<endl<<"sum = "<<sum<<endl;
                    flag++;
                }

            }

            sum = sum/flag;
           // cout<<sum<<endl;

            center.push_back(sum);
        }

        count--;
    }

    cout<<"intdex = "<<endl;
    for(int i=0; i<pointNum; i++)
    {
        cout<<label[i]<<endl;
    }
    return label;
}





