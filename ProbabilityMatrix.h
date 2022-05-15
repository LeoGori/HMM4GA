#ifndef HMM4GA_CLION_PROBABILITYMATRIX_H
#define HMM4GA_CLION_PROBABILITYMATRIX_H

#include<vector>
#include<map>
#include<iostream>

/*
 * Generic probability matrix containing floating points values.
 * It is supposed to be inherited by DynamicMatrix and StateTransitionMatrix, which share the
 * following properties and methods
 */

using namespace std;

class ProbabilityMatrix {
protected:

    vector<vector<float>> matrix;

public:

    // set the value of probability in position (state1 - state2) inside the matrix
    virtual void setProbability(int state1, int state2, float value);

    // returns the value of probability in position (state1 - state2) inside the matrix
    virtual float getProbability(int state1, int state2);
};

#endif //HMM4GA_CLION_PROBABILITYMATRIX_H