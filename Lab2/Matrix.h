#pragma once
#include <iostream>
#include <vector>

using namespace std;

class Matrix 
{
protected:
    vector<vector<double>> data;
    int rows;
    int cols;

public:
    Matrix(int m, int n) : rows(m), cols(n) 
    {
        data.resize(m, vector<double>(n, 0.0));
    }

    Matrix operator+(const Matrix& other) 
    {
        if (rows != other.rows || cols != other.cols) 
        {
            throw std::runtime_error("Dimensions of matrices do not match");
        }

        Matrix result(rows, cols);
        for (int i = 0; i < rows; ++i) 
        {
            for (int j = 0; j < cols; ++j) 
            {
                result.data[i][j] = data[i][j] + other.data[i][j];
            }
        }

        return result;
    }

    Matrix operator-(const Matrix& other) 
    {
        if (rows != other.rows || cols != other.cols) 
        {
            throw std::runtime_error("Dimensions of matrices do not match");
        }

        Matrix result(rows, cols);
        for (int i = 0; i < rows; ++i) 
        {
            for (int j = 0; j < cols; ++j) 
            {
                result.data[i][j] = data[i][j] - other.data[i][j];
            }
        }

        return result;
    }

    Matrix operator*(double scalar) 
    {
        Matrix result(rows, cols);
        for (int i = 0; i < rows; ++i) 
        {
            for (int j = 0; j < cols; ++j) 
            {
                result.data[i][j] = data[i][j] * scalar;
            }
        }

        return result;
    }

    std::vector<double> solveGauss(std::vector<double> b) 
    {
        
        return std::vector<double>(rows, 0.0);
    }
};

class SquareMatrix : public Matrix 
{
public:
    SquareMatrix(int n) : Matrix(n, n) 
    {
    
    }

    SquareMatrix inverse() 
    {
        if (rows != cols) 
        {
            throw std::runtime_error("The inverse matrix can only be defined for square matrices");
        }

        return SquareMatrix(rows);
    }

    double determinant() 
    {
        if (rows != cols) 
        {
            throw std::runtime_error("Determinant can only be defined for square matrices");
        }

        return 0.0;
    }
};

int main() 
{
    

    return 0;
}
