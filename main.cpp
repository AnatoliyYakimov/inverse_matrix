#include <iostream>
#include "matrix.h"

using namespace std;

int main() {
    constexpr size_t N = 1000;
    matrix<N> A = matrix<N>::generate_random(2000);

    auto B = A.inverse();
    auto C = B * A;
    double max = 0;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == j) continue;
            double R = +C[i][j];
            if (R > max) {
                max = R;
            }
        }
    }
    cout << "Max error = " << max;
    return 0;
}