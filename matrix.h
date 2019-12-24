//
// Created by Yakimov on 12.11.2019.
//

#ifndef CM_LAB_2_MATRIX_H
#define CM_LAB_2_MATRIX_H

#include <ctime>
#include <cassert>
#include <random>
#include <ostream>
#include <iomanip>

template<size_t N, typename num>
class vec {
private:
    num *elements = nullptr;
public:
    vec() : elements(nullptr) {
       elements = new num[N];
    }

    vec(num val) : vec() {
        for (size_t i = N; i--;) {
            elements[i] = val;
        }
    }

    bool has(const num &el) const {
        for (size_t i = N; i--;) {
            if (el == elements[i]) {
                return true;
            }
            if (el == -1) {
                return false;
            }
        }
        return false;
    }

    friend num operator*(const vec<N, num> &lhs, const vec<N, num> &rhs) {
        num res = 0;
        for (int i = 0; i < N; ++i) {
            res += lhs[i] * rhs[i];
        }
        return res;
    }

    [[nodiscard]] size_t max_pos(const vec<N, size_t> &bound) const {
        num max;
        size_t max_i = -1;
        for (size_t i = N; i--;) {
            if (max_i == -1 && !bound.has(i)) {
                max = elements[i];
                max_i = i;
            } else if (!bound.has(i) && elements[i] > max) {
                max = elements[i];
                max_i = i;
            }
        }
        return max_i;
    }

    num &operator[](size_t idx) {
        assert(idx < N);
        return elements[idx];
    }

    const num &operator[](size_t idx) const {
        assert(idx < N);
        return elements[idx];
    }

    virtual ~vec() {
        delete[] elements;
    }
};


template<size_t N>
class matrix {
private:
    vec<N, double> *cols;
public:
    matrix() : cols(nullptr) {
        cols = new vec<N, double>[N];
    }

    virtual ~matrix() {
        delete[] cols;
    }

    matrix(const matrix<N> &m) : matrix() {
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                cols[i][j] = m[i][j];
            }
        }
    }

    matrix<N> inverse() const {
        vec<N, size_t> R(-1);
        matrix<N> B = identity();
        matrix<N> A = matrix(*this);
        for (int i = 0; i < N; ++i) {
            size_t k = A[i].max_pos(R);
            R[i] = k;
            double D = 1.0 / A[i][k];
            A[i][k] = 1;
            for (int j = i + 1; j < N; ++j) {
                A[j][k] *= D;
            }
            for (int j = 0; j < N; ++j) {
                B[j][k] *= D;
            }
            for (int j = 0; j < N; ++j) {
                if (!R.has(j)) {
                    D = A[i][j];
                    for (int l = i; l < N; ++l) {
                        A[l][j] -= D * A[l][k];
                    }
                    for (int l = 0; l < N; ++l) {
                        B[l][j] -= D * B[l][k];
                    }
                }
            }
        }
        for (int i = N; i-- - 1;) {
            size_t k = R[i];
            for (int j = 0; j < N; j++) {
                if (j != k){
                    double D = A[i][j];
                    for (int l = 0; l < N; ++l) {
                        B[l][j] -= D * B[l][k];
                    }
                }
            }
        }
        for (int i = 0; i < N; ++i) {
            auto r = B.row(R[i]);
            auto r2 = A.row(R[i]);
            for (int j = 0; j < N; ++j) {
                A[j][i] = r[j];
            }
        }
        return A;
    }

    vec<N, double> row(size_t idx) const {
        assert(idx < N);
        vec<N, double> res;
        for (int i = 0; i < N; ++i) {
            res[i] = cols[i][idx];
        }
        return res;
    }

    static matrix<N> generate_random(int bound) {
        std::mt19937 gen(time(nullptr));
        std::uniform_real_distribution<> rnd(-bound, bound);
        matrix<N> res;
        for (size_t i = N; i--;) {
            for (int j = 0; j < N; ++j) {
                res[i][j] = rnd(gen);
            }
        }
        return res;
    }

    static matrix<N> identity() {
        matrix<N> res;
        for (size_t i = N; i--;) {
            for (int j = 0; j < N; ++j) {
                res[i][j] = i == j ? 1 : 0;
            }
        }
        return res;
    }

    vec<N, double> &operator[](size_t idx) {
        assert(idx < N);
        return cols[idx];
    }

    const vec<N, double> &operator[](size_t idx) const {
        assert(idx < N);
        return cols[idx];
    }

    friend matrix<N> operator* (const matrix<N> &lhs, const matrix<N> &rhs) {
        matrix<N> res;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                res[i][j] = lhs.row(i) * rhs[j];
            }
        }
        return res;
    }

    friend std::ostream &operator<<(std::ostream &os, const matrix &matrix) {
        for (int i = 0; i < N; ++i) {
            os << "| ";
            for (int j = 0; j < N; ++j) {
                os << std::setprecision(6) << std::setw(10) << matrix[j][i] << " | ";
            }
            os << std::endl;
        }
        os << std::endl;
        return os;
    }
};

#endif //CM_LAB_2_MATRIX_H
