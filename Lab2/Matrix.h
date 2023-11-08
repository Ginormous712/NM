#pragma once
#include <iostream>
#include <vector>
#include <iomanip>

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
                os << setprecision(4) << fixed << setw(10) << matrix.data[i][j] << " ";
            }
            os << endl;
        }
        cout.unsetf(ios::fixed);
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

    vector<double> getColumn(int columnIndex) const 
    {
        if (columnIndex < 0 || columnIndex >= cols) 
        {
            throw out_of_range("Invalid column index");
        }

        vector<double> column;
        column.reserve(rows);

        for (int i = 0; i < rows; i++) 
        {
            column.push_back(data[i][columnIndex]);
        }

        return column;
    }

    void gaussianSimpleElimination()
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

    void gaussianElimination()
    {
        for (int i = 0; i < getRows(); i++)
        {
            double maxVal = abs(data[i][i]);
            int maxRow = i;
            for (int k = i + 1; k < getRows(); k++)
            {
                if (abs(data[k][i]) > maxVal)
                {
                    maxVal = abs(data[k][i]);
                    maxRow = k;
                }
            }

            if (maxRow != i)
            {
                for (int k = i; k < getCols(); k++)
                {
                    swap(data[i][k], data[maxRow][k]);
                }
            }

            double pivot = data[i][i];

            for (int k = i; k < getCols(); k++)
            {
                data[i][k] /= pivot;
            }

            for (int j = 0; j < getRows(); j++)
            {
                if (j != i)
                {
                    double factor = data[j][i];
                    for (int k = i; k < getCols(); k++)
                    {
                        data[j][k] -= factor * data[i][k];
                    }
                }
            }
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
    
    double determinant()
    {
        SquareMatrix copy(*this);
        double det = 1.0;
        copy.gaussianSimpleElimination();

        for (int i = 0; i < copy.rows; i++)
        {
            det *= copy.data[i][i];
        }
        return det;
    }
    
    vector<double> solveLinearSystemGauss(vector<double> b)
    {
        if (rows != b.size())
        {
            throw runtime_error("The system of equations has no unique solution");
        }

        Matrix system(rows, cols + 1);
        for (size_t i = 0; i < system.getRows(); i++)
        {
            for (size_t j = 0; j < system.getCols(); j++)
            {
                if (j == system.getCols() - 1)
                {
                    system.setElement(i, j, b[i]);
                }
                else
                {
                    system.setElement(i, j, getElement(i, j));
                }
            }
        }
        
        //cout << system;
        system.gaussianElimination();
        
        return system.getColumn(system.getCols() - 1);
    }
    
    SquareMatrix findInverseMatrixGauss() 
    {
        int n = getRows();
        vector<vector<double>> augmented(n, vector<double>(2 * n, 0.0));

        for (int i = 0; i < n; i++) 
        {
            for (int j = 0; j < n; j++) 
            {
                augmented[i][j] = getElement(i, j);
                if (i == j) 
                {
                    augmented[i][j + n] = 1.0;
                }
            }
        }

        Matrix augmetedMatrix(augmented);
        augmetedMatrix.gaussianElimination();

        SquareMatrix inverseMatrix(n);

        for (int i = 0; i < n; i++) 
        {
            for (int j = 0; j < n; j++) 
            {
                inverseMatrix.setElement(i, j, augmetedMatrix.getElement(i, j + n));
            }
        }

        return inverseMatrix;
    }
};


