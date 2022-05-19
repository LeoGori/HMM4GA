#ifndef HMM4GA_CLION_STATETRANSITIONMATRIX_H
#define HMM4GA_CLION_STATETRANSITIONMATRIX_H

#include <vector>
#include <map>
#include <string>
#include "ProbabilityMatrix.h"

/*
    The state transition matrix holds the probability according to which a state switch to another.
    As described in the reference document, the states that are considered when analyzing
    a base of the read and comparing it with a base of the haplotype can be resumed in:
    - Matching state (M)
    - Insertion state (I)
    - Deletion state (D)
 */

using namespace std;


class StateTransitionMatrix : public ProbabilityMatrix{

private:

    // char to int mapping that is useful when iterating over the states ('M':0, 'I':1, 'D':2)
    map<char, int> states;

public:
    StateTransitionMatrix();

    // returns the transition probability for switching from state char1 to state char2
    virtual float getProbability(char char1, char char2);

    // returns the whole transition matrix
    vector<vector<float>> getTransitionMatrix();

    // fills the matrix with random probability values, making sure that transition probabilities
    // from the same states sum up to 1
    // possible improvement: assign defined transition probabilities, as the method describe in "references/2017_fpga_pairhmm.pdf" (?)
    void populateMatrix();

    // output stream operator overloading, for easily reading the transition matrix
    friend ostream& operator<<(ostream& out, const StateTransitionMatrix& stm);

//    void print();

};

#endif //HMM4GA_CLION_STATETRANSITIONMATRIX_H