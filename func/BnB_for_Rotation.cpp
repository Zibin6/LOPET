#include "mex.h"
#include <cmath> 
#include <math.h>
#include <vector>
//#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#define M_PI 3.1415926

void remove_columns_with_NaN(const MatrixXd& inputMatrix, MatrixXd& outputMatrix) {
    int numCols = inputMatrix.cols();
    std::vector<int> validColumns;

    for (int i = 0; i < numCols; ++i) {
        if (!(inputMatrix.col(i).array().isNaN().any())) {
            validColumns.push_back(i);
        }
    }

    int numValidCols = validColumns.size();
    outputMatrix.conservativeResize(inputMatrix.rows(), numValidCols);

    for (int i = 0; i < numValidCols; ++i) {
        outputMatrix.col(i) = inputMatrix.col(validColumns[i]);
    }

    return;
}

void convert_input_eigen(MatrixXd& X, MatrixXd& Y, const mxArray *prhs0, const mxArray *prhs1)
{
    MatrixXd X_ori, Y_ori;
    int nx = mxGetN(prhs0);
    int ny = mxGetN(prhs1);
    X_ori = MatrixXd::Zero(3, nx);
    Y_ori = MatrixXd::Zero(3, ny);
    
    double *X_pt_dbl, *Y_pt_dbl;
    float *X_pt_flt, *Y_pt_flt;
    const mxArray *inputArray0 = prhs0;
    if (mxIsSingle(inputArray0))
    {
        X_pt_flt = (float *)mxGetData(prhs0);
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < nx; j++)
            {
                X_ori(i, j) = X_pt_flt[j*3+i];
            }
        }
    }
    else
    {
        X_pt_dbl = (double *)mxGetData(prhs0);
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < nx; j++)
            {
                X_ori(i, j) = X_pt_dbl[j*3+i];
            }
        }
    }
    const mxArray *inputArray1 = prhs1;
    if (mxIsSingle(inputArray1))
    {
        Y_pt_flt = (float *)mxGetData(prhs1);
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                Y_ori(i, j) = Y_pt_flt[j*3+i];
            }
        }
    }
    else
    {
        Y_pt_dbl = (double *)mxGetData(prhs1);
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                Y_ori(i, j) = Y_pt_dbl[j*3+i];
            }
        }
    }
    remove_columns_with_NaN(X_ori, X);
    remove_columns_with_NaN(Y_ori, Y);

    if (X.cols() == 0 || Y.cols() == 0)
        mexErrMsgTxt("all columns of the input(s) have NaN!");
    return;
}

int bitget(int num, int bit) {
    return (num >> (bit - 1)) & 1;
}

void obtain_new_branch(Eigen::Matrix3d& M, Eigen::Matrix<double,6,8>& new_branch)
{
    for(int i = 0; i < 8; i++)
    {
        new_branch(0, i) = M(0, bitget(i+1, 1));
        new_branch(1, i) = M(1, bitget(i+1, 2));
        new_branch(2, i) = M(2, bitget(i+1, 3));
        new_branch(3, i) = M(0, bitget(i+1, 1)+1);
        new_branch(4, i) = M(1, bitget(i+1, 2)+1);
        new_branch(5, i) = M(2, bitget(i+1, 3)+1);
    }
    return;
}

void rotationVectorToMatrix(const Eigen::Vector3d& rotationVector, Eigen::Matrix3d& rotationMatrix)
{
    double angle = rotationVector.norm();
    Eigen::Vector3d axis = rotationVector.normalized();
    rotationMatrix = Eigen::AngleAxisd(angle, axis).toRotationMatrix();
    rotationMatrix.transposeInPlace();
    return;
}

void combine_branch_matrix(
    Eigen::MatrixXd& branch,
    Eigen::Matrix<double,6,8>& new_branch)
{
    int n_col = branch.cols();
    branch.conservativeResize(branch.rows(), branch.cols()+8);
    branch.block<6, 8>(0, n_col) = new_branch;
    return;
}

void combine_branch_ul_matrix(
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>& branch_ul,
    Eigen::Matrix<int,1,8>& new_upper, 
    Eigen::Matrix<int,1,8>& new_lower)
{
    int n_col = branch_ul.cols();
    branch_ul.conservativeResize(branch_ul.rows(), branch_ul.cols()+8);
    branch_ul.block<1, 8>(0, n_col) = new_upper;
    branch_ul.block<1, 8>(1, n_col) = new_lower;
    return;
}

void rearrange_branch_data(
    Eigen::MatrixXd& branch,
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>& branch_ul,
    int ind_upper, int best_lower)
{
    int n = branch.cols();
    Eigen::VectorXi boolVector(n);
    boolVector.setOnes();
    boolVector(ind_upper) = 0;
    for(int i = 0; i < n; i++)
    {
        if(i == ind_upper)
            continue;
        if(branch_ul(0,i) < best_lower)
        {
            boolVector(i) = 0;
        }
    }

    int idx = 0;
    for(int i = 0; i < n; i++)
    {
        if(boolVector(i))
        {
            branch.block<6, 1>(0, idx) = branch.block<6, 1>(0, i);
            branch_ul.block<2, 1>(0, idx) = branch_ul.block<2, 1>(0, i);
            idx++;
        }
    }
    branch.conservativeResize(branch.rows(), idx);
    branch_ul.conservativeResize(branch_ul.rows(), idx);
    
    return;
}
Eigen::Matrix3d core_comp(MatrixXd& X, MatrixXd& Y, double epsilon, const Eigen::Matrix<double,6,1>& initial_best_branch)
{
    Eigen::Matrix3d R_opt;
    R_opt.setZero();
    Eigen::Matrix<double,3,1> R_best;
    R_best.setZero();
    Eigen::Matrix<double,6,1> best_branch = initial_best_branch; 
    Eigen::MatrixXd branch = Eigen::MatrixXd::Zero(6,8);
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> branch_ul;
    branch_ul.resize(2, 8);
    branch_ul.setZero();

    Eigen::Matrix<int,1,8> new_upper, new_lower;
    new_upper.setZero();
    new_lower.setZero();
    int best_lower = 0;

    Eigen::Matrix3d M;
    Eigen::Matrix<double,6,8> new_branch;
    Eigen::Matrix<double,6,1> B;
    Eigen::Matrix3d R_c;
    double delta_B;

    while(true)
    {
        M.col(0) = best_branch.head(3);
        M.col(1) = (best_branch.head(3) + best_branch.segment<3>(3))*0.5;
        M.col(2) = best_branch.segment<3>(3);
        obtain_new_branch(M, new_branch);

    for(int s = 0; s < 8; s++)
    {
        B = new_branch.col(s);
        Eigen::Vector3d rv = (B.head(3) + B.segment<3>(3)) * 0.5;
        rotationVectorToMatrix(rv, R_c);

        Eigen::Vector3d tmp = B.segment<3>(3) - B.head(3);
        delta_B = tmp.norm() * 0.5;

        
        Eigen::MatrixXd ematrix = (M_PI / 2 - Eigen::acos(((R_c * X).transpose() * Y).array())).array().abs();


        Eigen::Matrix<bool, Eigen::Dynamic, 1> temp_min_check = (ematrix.colwise().minCoeff().array() <= epsilon);

        Eigen::Matrix<bool, Eigen::Dynamic, 1> temp_upper_check = ((ematrix.colwise().minCoeff().array() - delta_B) <= epsilon);

        new_upper(s) = temp_upper_check.cast<int>().sum();
        new_lower(s) =  temp_min_check.cast<int>().sum();
    }
        
        combine_branch_matrix(branch, new_branch);
        combine_branch_ul_matrix(branch_ul, new_upper, new_lower);

        int best_upper = -1;
        int ind_upper = -1;
        int new_best_lower = -1;
        int ind_lower = -1;
        for(int i = 0; i < branch_ul.cols(); i++)
        {
            if(branch_ul(0, i) > best_upper)
            {
                best_upper = branch_ul(0, i);
                ind_upper = i;
            }
            if(branch_ul(1, i) > new_best_lower)
            {
                new_best_lower = branch_ul(1, i);
                ind_lower = i;
            }
        }

        if(best_lower < new_best_lower)
        {
            best_lower=new_best_lower;
            R_best = (branch.block<3, 1>(0, ind_lower) + branch.block<3, 1>(3, ind_lower))*0.5;
        }
        best_branch = branch.block<6, 1>(0, ind_upper);

        rearrange_branch_data(branch, branch_ul, ind_upper, best_lower);

        if(best_upper <= best_lower)
        {
            break;
        }

    }

    rotationVectorToMatrix(R_best, R_opt);
    return R_opt;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
    if (nrhs < 3)
        mexErrMsgTxt("at least 3 inputs are required!");
    if (nlhs < 1)
        mexErrMsgTxt("at least 1 output is required!");

    if (mxIsEmpty(prhs[0]) || mxIsEmpty(prhs[1]) || mxIsEmpty(prhs[2]))
        mexErrMsgTxt("input parameter should not be an empty array!");

    if (mxGetM(prhs[0])!=3)
        mexErrMsgTxt("the 1st input parameter should be 3 by n!");

    if (mxGetM(prhs[1])!=3)
        mexErrMsgTxt("the 2nd input parameter should be 3 by n!");

    double epsilon = mxGetScalar(prhs[2]);

    if (nrhs < 4) {  
        mexErrMsgTxt("at least 4 inputs are required for the initial best_branch!");  
    }  

    Eigen::Matrix<double,6,1> best_branch_eigen;  
    if (mxIsDouble(prhs[3]) && mxGetM(prhs[3]) == 6 && mxGetN(prhs[3]) == 1) {  
        double* data = mxGetPr(prhs[3]);  
        for (int i = 0; i < 6; ++i) {  
            best_branch_eigen(i) = data[i];  
        }  
    } else {  
        mexErrMsgTxt("The 4th input parameter should be a 6x1 double array for the initial best_branch!");  
    } 


    int n_col = mxGetN(prhs[0]);
    plhs[0] = mxCreateDoubleMatrix(3, 3, mxREAL);

    MatrixXd X, Y;

    convert_input_eigen(X, Y, prhs[0], prhs[1]);
    Eigen::Matrix3d R_opt = core_comp(X, Y, epsilon, best_branch_eigen);

    Eigen::MatrixXd ematrix = (M_PI / 2 - Eigen::acos(((R_opt * X).transpose() * Y).array())).array().abs();
    Eigen::VectorXd min_per_col = ematrix.colwise().minCoeff();  
    double r_index0 = min_per_col.mean();  



    double *rot_m = mxGetPr(plhs[0]);
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            rot_m[j*3+i] = R_opt(i, j);
        }
    }
    // 为 ematrix 分配输出矩阵（假设 ematrix 的大小与 X 或 Y 的列数相同）  
    plhs[1] = mxCreateDoubleMatrix(ematrix.rows(), ematrix.cols(), mxREAL);  
    double *temp_data = mxGetPr(plhs[1]);  
    for (int i = 0; i < ematrix.rows(); i++)  
    {  
        for (int j = 0; j < ematrix.cols(); j++)  
        {  
            temp_data[j*ematrix.rows() + i] = ematrix(i, j); // 注意行列索引的顺序  
        }  
    }  
    plhs[2] = mxCreateDoubleScalar(r_index0);  

    return;
}