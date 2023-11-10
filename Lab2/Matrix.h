#pragma once
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

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

    vector<vector<double>> getData() const
    {
        return data;
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

            if (maxVal == 0.0)
            {
                throw std::runtime_error("Matrix has no unique solution");
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
    
    bool isTridiagonal() const
    {
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                if (std::abs(i - j) > 1 && data[i][j] != 0)
                {
                    return false;
                }
            }
        }
        return true;
    }
    
    void solveAndPrint(const vector<double>& b, const vector<double>& solution) const
    {
        cout << "Vector b:\n";
        printVector(b, false);

        cout << "Solution using the specified method:\n";
        printVector(solution, true);
    }

    void printVector(const vector<double>& vector, bool flag) const 
    {
        char c;
        if (flag)
        {
            c = 'x';
        }
        else
        {
            c = 'b';
        }
        for (size_t i = 0; i < vector.size(); i++)
        {
            cout << c << i + 1 << " = " << vector[i] << ";\t";
        }
        cout << endl;
    }

    Matrix transpose() const 
    {
        int rows = data.size();
        int cols = data[0].size();

        Matrix result(vector<vector<double>>(cols, vector<double>(rows, 0.0)));

        for (int i = 0; i < rows; ++i) 
        {
            for (int j = 0; j < cols; ++j) 
            {
                result.data[j][i] = data[i][j];
            }
        }

        return result;
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
        cout << "The square matrix was created successfully." << endl;
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
    
    bool isSingular()
    {
        double d = determinant();
        double epsilon = 1e-10;

        return abs(d) < epsilon;
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

    static double calculateEuclideanNorm(const vector<double>& vector) 
    {
        double norm = 0.0;
        for (double component : vector) {
            norm += component * component;
        }
        return sqrt(norm);
    }

    static double calculateAbsoluteError(const vector<double>& x1, const vector<double>& x2)
    {
        int size = x1.size();
        double error = 0.0;

        for (int i = 0; i < size; i++) 
        {
            error = max(error, abs(x1[i] - x2[i]));
        }

        return error;
    }

    bool isSymmetric() const
    {
        Matrix transposed = transpose();
        return data == transposed.getData();
    }

    bool isConvergent() const
    {
        int size = data.size();
        for (int i = 0; i < size; i++)
        {
            double diagonalElement = abs(data[i][i]);
            double sum = 0.0;
            for (int j = 0; j < size; ++j)
            {
                if (j != i)
                {
                    sum += abs(data[i][j]);
                }
            }
            if (diagonalElement <= sum)
            {
                return false;
            }
        }
        return true;
    }



    vector<double> gaussSeidel(const vector<double>& b, int maxIterations, double epsilon) const 
    {
        int n = rows;
        vector<double> x(n, 0.0);

        cout << "Vector b:\n";
        printVector(b, false);

        double e = epsilon;
        int pow = 0;
        while (e != 1)
        {
            pow++;
            e *= 10;
        }

        cout << setprecision(pow) << fixed;
        double error1 = 0;
        for (int iter = 0; iter < maxIterations; iter++) 
        {
            vector<double> x_old = x;

            for (int i = 0; i < n; i++) 
            {
                double sigma = 0.0;

                for (int j = 0; j < n; j++) 
                {
                    if (j != i) 
                    {
                        sigma += data[i][j] * x[j];
                    }
                }

                x[i] = (b[i] - sigma) / data[i][i];
            }
            cout << iter + 1 << "| ";
            printVector(x, true);
            double error2 = calculateEuclideanNorm(x);
            cout << "Error: " << error2 - error1 << endl;
            //double error = calculateAbsoluteError(x, x_old);
            if (abs(error2 - error1) < epsilon) 
            {
                cout << "*****************     STOP     ******************\n";
                cout << "Result in " << iter + 1 << " iterations.\n";
                cout.unsetf(ios::fixed);
                return x;
            }
            error1 = error2;
        }

        cout << "We reached the maximum number of iterations.\n";
        cout.unsetf(ios::fixed);
        return x;
    }
};


class TridiagonalMatrix : public SquareMatrix
{
public:
    TridiagonalMatrix(const vector<vector<double>>& matrixData) : SquareMatrix(matrixData)
    {
        if (!isTridiagonal())
        {
            throw invalid_argument("The matrix is not tridiagonal");
        }
    }
    
    static bool checkStability(vector<vector<double>> data) 
    {
        vector<double> a, b, c;
        a.reserve(data[0].size());
        b.reserve(data[0].size());
        c.reserve(data[0].size());

        for (int i = 0; i < data[0].size(); i++)
        {
            if (i > 0)
            {
                a.push_back(data[i][i - 1]);
            }
            else
            {
                a.push_back(0.0);
            }

            b.push_back(data[i][i]);

            if (i < data[0].size() - 1)
            {
                c.push_back(data[i][i + 1]);
            }
            else
            {
                c.push_back(0.0);
            }
        }

        bool flag = false;
        for (int i = 0; i < a.size(); i++)
        {
            if (abs(b[i]) < abs(a[i]) + abs(c[i]))
            {
                return false;
            }
            if (abs(b[i]) > abs(a[i]) + abs(c[i]))
            {
                flag = true;
            }

        }
        return flag;
    }
    /*
    vector<double> tridiagonalSolve(const vector<double>& d) const 
    {
        

        if (rows != d.size()) {
            throw invalid_argument("The dimension of the right-hand side vector does not correspond to the dimension of the matrix");
        }

        vector<double> a, b, c;
        a.reserve(rows);
        b.reserve(rows);
        c.reserve(rows);

        for (int i = 0; i < rows; i++)
        {
            if (i > 0) 
            {
                a.push_back(data[i][i - 1]);
            }
            else 
            {
                a.push_back(0.0);
            }
            
            b.push_back(data[i][i]);
            
            if (i < rows - 1) 
            {
                c.push_back(data[i][i + 1]);
            }
            else 
            {
                c.push_back(0.0);
            }
        }
        if (!checkStability(a, b, c))
        {
            throw std::runtime_error("The tridiagonal matrix is not diagonally dominant, and the stability condition is not satisfied.");
        }
        vector<double> alpha(rows);
        vector<double> beta(rows);
        vector<double> x(rows);

        alpha[0] = b[0];
        beta[0] = d[0] / alpha[0];

        for (int i = 1; i < rows; i++) 
        {
            alpha[i] = b[i] - a[i] * c[i - 1] / alpha[i - 1];
            beta[i] = (d[i] - a[i] * beta[i - 1]) / alpha[i];
        }

        x[rows - 1] = beta[rows - 1];

        for (int i = rows - 2; i >= 0; i--) 
        {
            x[i] = beta[i] - c[i] * x[i + 1] / alpha[i];
        }

        return x;
    }
    //*/
    
    vector<double> tridiagonalSolve(const vector<double>& d) const 
    {
        if (rows != d.size() || rows != data[0].size() || rows != data.size()) 
        { 
            throw invalid_argument("Matrix and vector dimensions do not match");
        }
        /*
        if (!checkStability(getData()))
        {
            throw std::runtime_error("The tridiagonal matrix is not diagonally dominant, and the stability condition is not satisfied.");
        }
        */
        vector<double> alpha(rows);
        vector<double> beta(rows);
        vector<double> z(rows);

        alpha[0] = -data[0][1] / data[0][0];
        beta[0] = d[0] / data[0][0];

        for (int i = 1; i < rows - 1; i++) 
        {
            z[i] = -data[i][i] - alpha[i - 1] * data[i][i - 1];
            alpha[i] = data[i][i + 1] / z[i];
            beta[i] = (-d[i] + data[i - 1][i] * beta[i - 1]) / z[i];
            
            
        }
        z[rows - 1] = -data[rows - 1][rows - 1] - alpha[rows - 2] * data[rows - 1][rows - 2];

        vector<double> x(rows);
        x[rows - 1] = (-d[rows-1] + data[rows - 1][rows - 2]*beta[rows-2])/z[rows - 1];

        for (int i = rows - 2; i >= 0; i--) 
        {
            x[i] = alpha[i] * x[i + 1] + beta[i];
        }

        return x;
    }

};


