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

    Matrix(const vector<vector<double>>& data) 
    {
        rows = data.size();
        
        if (rows > 0) 
        {
            cols = data[0].size();
        }
        else 
        {
            cols = 0;
        }

        this->data = data;

        for (int i = 1; i < rows; i++) 
        {
            if (data[i].size() != cols) 
            { 
                throw invalid_argument("Vector of vectors has different number of columns in different rows");
            }
        }
    }

    Matrix operator+(const Matrix& other) 
    {
        if (rows != other.rows || cols != other.cols) 
        {
            throw std::runtime_error("Dimensions of matrices do not match");
        }

        Matrix result(rows, cols);
        for (int i = 0; i < rows; i++) 
        {
            for (int j = 0; j < cols; j++) 
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
        for (int i = 0; i < rows; i++) 
        {
            for (int j = 0; j < cols; j++) 
            {
                result.data[i][j] = data[i][j] - other.data[i][j];
            }
        }

        return result;
    }

    Matrix operator*(double scalar) 
    {
        Matrix result(rows, cols);
        for (int i = 0; i < rows; i++) 
        {
            for (int j = 0; j < cols; j++) 
            {
                result.data[i][j] = data[i][j] * scalar;
            }
        }

        return result;
    }

    friend ostream& operator<<(ostream& os, const Matrix& matrix) 
    {
        for (int i = 0; i < matrix.rows; i++) 
        {
            for (int j = 0; j < matrix.cols; j++) 
            {
                os << matrix.data[i][j] << " ";
            }
            os << endl;
        }

        return os;
    }

    friend istream& operator>>(istream& is, Matrix& matrix) 
    {
        for (int i = 0; i < matrix.rows; i++) 
        {
            for (int j = 0; j < matrix.cols; j++) 
            {
                is >> matrix.data[i][j];
            }
        }
        return is;
    }

    double getElement(int row, int col) const 
    {
        if (row >= 0 && row < rows && col >= 0 && col < cols) 
        {
            return data[row][col];
        }
        else 
        {
            throw std::out_of_range("Invalid item index");
        }
    }

    void setElement(int row, int col, double value) 
    {
        if (row >= 0 && row < rows && col >= 0 && col < cols) 
        {
            data[row][col] = value;
        }
        else 
        {
            throw std::out_of_range("Invalid item index");
        }
    }

    int getRows() const 
    {
        return rows;
    }

    int getCols() const 
    {
        return cols;

    }

    void setColumn(int colIndex, const vector<double>& column) 
    {
        int n = getRows();
        
        if (colIndex < 0 || colIndex >= n || column.size() != n) 
        {
            throw out_of_range("Invalid indexes or column size");
        }

        for (int i = 0; i < n; ++i) 
        {
            setElement(i, colIndex, column[i]);
        }
    }
};

class SquareMatrix : public Matrix
{
public:
    SquareMatrix(int n) : Matrix(n, n)
    {

    }

    SquareMatrix(const vector<vector<double>>& data) : Matrix(data) 
    {
        if (getRows() != getCols()) 
        {
            throw invalid_argument("A square matrix must have the same number of rows and columns.");
        }
    }
    ///*
    void gaussianElimination()
    {
        for (int i = 0; i < rows; i++)
        {
            int maxRow = i;
            double maxElement = std::abs(data[i][i]);

            for (int j = i + 1; j < rows; j++)
            {
                double absVal = std::abs(data[j][i]);
                if (absVal > maxElement)
                {
                    maxRow = j;
                    maxElement = absVal;
                }
            }

            if (maxElement == 0.0)
            {
                throw std::runtime_error("Matrix has no unique solution");
            }

            if (maxRow != i)
            {
                std::swap(data[i], data[maxRow]);
            }

            for (int j = i + 1; j < rows; j++)
            {
                double factor = data[j][i] / data[i][i];

                for (int k = i; k < cols; k++)
                {
                    data[j][k] -= factor * data[i][k];
                }
            }
        }
    }
    //*/
    /*
    void gaussianElimination() {
        for (int i = 0; i < rows; ++i) {
            if (data[i][i] == 0.0) {
                throw std::runtime_error("ћатриц€ не маЇ ун≥кального розв'€зку");
            }

            for (int j = i + 1; j < rows; ++j) {
                double factor = data[j][i] / data[i][i];
                for (int k = i; k < cols; ++k) {
                    data[j][k] -= factor * data[i][k];
                }
            }
        }
    }
    */
    double determinant()
    {
        double det = 1.0;
        gaussianElimination();

        for (int i = 0; i < rows; i++)
        {
            det *= data[i][i];
        }
        return det;
    }
    
    vector<double> solveLinearSystemGauss(vector<double> b)
    {
        if (rows != b.size())
        {
            throw runtime_error("The system of equations has no unique solution");
        }

        gaussianElimination();

        vector<double> x(rows, 0.0);

        for (int i = rows - 1; i >= 0; i--)
        {
            double sum = 0.0;
            for (int j = i + 1; j < rows; j++)
            {
                sum += data[i][j] * x[j];
            }

            x[i] = (b[i] - sum) / data[i][i];
        }

        return x;
    }

    SquareMatrix findInverseMatrix() 
    {
        SquareMatrix A_inv(getRows());
        SquareMatrix A_copy(*this); 

        try 
        {
            double det = determinant();

            if (det == 0.0) 
            {
                throw runtime_error("Matrix has no inverse matrix (determinant = 0)");
            }

            A_copy.gaussianElimination();

            for (int i = 0; i < getRows(); i++) 
            {
                vector<double> b(getRows(), 0.0);
                b[i] = 1.0; 

                vector<double> x = A_copy.solveLinearSystemGauss(b);
                A_inv.setColumn(i, x);
            }
        }
        catch (const runtime_error& e) 
        {
            cerr << "Error: " << e.what() << endl;
        }

        return A_inv;
    }

};


