#include "ProbabilityMatrix.h"
#include<bits/stdc++.h>




void ProbabilityMatrix::setProbability(int state1, int state2, float value) {

    matrix[state1][state2] = value;
}

float ProbabilityMatrix::getProbability(int state1, int state2) {
    return matrix[state1][state2];
}