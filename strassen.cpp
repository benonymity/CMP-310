// 
// By Benjamin Bassett (benonymity) on 11.7.24
// 
// For CMP-310 with Dr. Serang at Hillsdale College
// https://github.com/benonymity/CMP-310
// 

#include <vector>
#include <iostream>
#include <stdexcept>

using namespace std;

// I stole most of the matrix stuff from the internet——C++ not supporting matrices was the hardest part of doing this
template<typename T>
class Matrix {
private:
    vector<vector<T>> data;
    size_t rows, cols;

public:
    Matrix(size_t r, size_t c) : rows(r), cols(c), data(r, vector<T>(c)) {}
    
    Matrix(const vector<vector<T>>& input) : data(input) {
        rows = input.size();
        cols = rows > 0 ? input[0].size() : 0;
    }

    T& operator()(size_t i, size_t j) {
        return data[i][j];
    }

    const T& operator()(size_t i, size_t j) const {
        return data[i][j];
    }

    Matrix operator+(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw invalid_argument("Matrix dimensions must match for addition");
        }
        Matrix result(rows, cols);
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                result(i,j) = data[i][j] + other.data[i][j];
            }
        }
        return result;
    }

    Matrix operator-(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw invalid_argument("Matrix dimensions must match for subtraction");
        }
        Matrix result(rows, cols);
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                result(i,j) = data[i][j] - other.data[i][j];
            }
        }
        return result;
    }

    Matrix operator*(const Matrix& other) const {
        if (cols != other.rows) {
            throw invalid_argument("Matrix dimensions must match for multiplication");
        }
        Matrix result(rows, other.cols);
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < other.cols; j++) {
                result(i,j) = T();
                for (size_t k = 0; k < cols; k++) {
                    result(i,j) += data[i][k] * other.data[k][j];
                }
            }
        }
        return result;
    }

    size_t getRows() const { return rows; }
    size_t getCols() const { return cols; }
    const vector<vector<T>>& getData() const { return data; }

    Matrix getQuadrant(int quad) const {
        // Check if matrix can be divided into quadrants
        if (rows % 2 != 0 || cols % 2 != 0) {
            throw invalid_argument("Matrix dimensions must be even for Strassen algorithm");
        }
        
        size_t half_rows = rows/2;
        size_t half_cols = cols/2;
        Matrix result(half_rows, half_cols);
        
        size_t row_start = (quad == 2 || quad == 4) ? half_rows : 0;
        size_t col_start = (quad == 3 || quad == 4) ? half_cols : 0;
        
        for(size_t i = 0; i < half_rows; i++) {
            for(size_t j = 0; j < half_cols; j++) {
                result(i,j) = data[i + row_start][j + col_start];
            }
        }
        return result;
    }

    static Matrix join(const Matrix& c11, const Matrix& c12, 
                      const Matrix& c21, const Matrix& c22) {
        size_t n = c11.getRows();
        Matrix result(2*n, 2*n);
        
        for(size_t i = 0; i < n; i++) {
            for(size_t j = 0; j < n; j++) {
                result(i,j) = c11(i,j);
                result(i,j + n) = c12(i,j);
                result(i + n,j) = c21(i,j);
                result(i + n,j + n) = c22(i,j);
            }
        }
        return result;
    }

    void print(bool first_ten = false) const {
        if (first_ten) {
            int count = 0;
            for (const auto& row : data) {
                cout << " [";
                for (const auto& val : row) {
                    cout << val << " ";
                    if (++count >= 10) {
                        cout << "...]";
                        return;
                    }
                }
            }
        } else {
            cout << "[" << endl;
            for (const auto& row : data) {
                cout << " [";
                for (const auto& val : row) {
                    cout << val << " ";
                }
                cout << "]" << endl;
            }
        }
    }
};

// Not sure if the templating stuff actually works, but figured I'd try to make it flexible
template<typename T>
Matrix<T> strassen(const Matrix<T>& A, const Matrix<T>& B) {
    size_t n = A.getRows();
    
    // Check if matrices are square and dimensions match
    if (A.getRows() != A.getCols() || B.getRows() != B.getCols() || A.getCols() != B.getRows()) {
        throw invalid_argument("Matrices must be square and dimensions must match");
    }
    
    // Check if dimensions are powers of 2
    if ((n & (n-1)) != 0) {
        throw invalid_argument("Matrix dimensions must be powers of 2 for Strassen algorithm");
    }
    
    if(n <= 64) { // Leaf size
        return A * B;
    }
    
    // Get quadrants
    Matrix<T> A11 = A.getQuadrant(1);
    Matrix<T> A12 = A.getQuadrant(2);
    Matrix<T> A21 = A.getQuadrant(3);
    Matrix<T> A22 = A.getQuadrant(4);
    
    Matrix<T> B11 = B.getQuadrant(1);
    Matrix<T> B12 = B.getQuadrant(2);
    Matrix<T> B21 = B.getQuadrant(3);
    Matrix<T> B22 = B.getQuadrant(4);
    
    // Calculate the 7 products recursively
    Matrix<T> P1 = strassen(A11 + A12, B22);
    Matrix<T> P2 = strassen(B12 - B22, A11);
    Matrix<T> P3 = strassen(A21 + A22, B11);
    Matrix<T> P4 = strassen(B21 - B11, A22);
    Matrix<T> P5 = strassen(A11 + A22, B11 + B22);
    Matrix<T> P6 = strassen(A12 - A22, B21 + B22);
    Matrix<T> P7 = strassen(A21 - A11, B11 + B12);
    
    // Calculate the quadrants of the result
    Matrix<T> C11 = P5 + P6 + P4 - P1;
    Matrix<T> C12 = P1 + P2;
    Matrix<T> C21 = P3 + P4;
    Matrix<T> C22 = P5 + P7 + P2 - P3;
    
    // Join the quadrants into the final result
    return Matrix<T>::join(C11, C12, C21, C22);
}

int main() {
    // Test different matrix sizes
    std::vector<size_t> sizes = {2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};
    
    for (size_t n : sizes) {
        // Create random matrices
        Matrix<int> A(n, n);
        Matrix<int> B(n, n);
        
        // Fill with random integers between 0 and 256
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                A(i,j) = rand() % 256;
                B(i,j) = rand() % 256;
            }
        }
        
        // Time Strassen's algorithm
        auto start = std::chrono::high_resolution_clock::now();
        Matrix<int> C_strassen = strassen(A, B);
        auto end = std::chrono::high_resolution_clock::now();
        auto strassen_duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        
        // Time naive multiplication
        start = std::chrono::high_resolution_clock::now();
        Matrix<int> C_naive = A * B;
        end = std::chrono::high_resolution_clock::now();
        auto naive_duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

        // Verify that both methods produce the same result
        bool matrices_equal = true;
        for (size_t i = 0; i < n && matrices_equal; i++) {
            for (size_t j = 0; j < n && matrices_equal; j++) {
                if (C_strassen(i,j) != C_naive(i,j)) {
                    matrices_equal = false;
                }
            }
        }
        
        if (!matrices_equal) {
            std::cout << "ERROR: Results do not match for size " << n << "x" << n << std::endl;
            C_strassen.print(true);
            C_naive.print(true);
            exit(1);
        }
        
        std::cout << "Matrix size: " << n << "x" << n << std::endl;
        std::cout << "Strassen: " << strassen_duration.count() << " microseconds" << std::endl;
        std::cout << "Naive: " << naive_duration.count() << " microseconds" << std::endl;
        C_strassen.print(true);
        C_naive.print(true);
        std::cout << "\n-------------------------" << std::endl;
    }
    
    return 0;
}
