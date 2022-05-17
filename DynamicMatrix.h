#ifndef HMM4GA_CLION_DYNAMICMATRIX_H
#define HMM4GA_CLION_DYNAMICMATRIX_H

#include "ProbabilityMatrix.h"

/*
 * the dynamic matrix that is supposed to be dynamically built up by the algorithm
 * and that is supposed to grow up while new nucleotides are introduced in the sequences (?)
 */

class DynamicMatrix : public ProbabilityMatrix {

private:

    // the value that is set for initializing the imaginary positioned cells of the matrix
    // it "corresponds to starting at an imaginary position one base before the read start in a deletion state"
    float initializationValue;
public:

    explicit DynamicMatrix(float initValue);

    // returns the initializationValue
    float getInitializationValue() const;

    // sets the initializationValue
    void setInitializationValue(float initializationValue);

    // adds a 0-filled row to the matrix
    void addRow();

    // adds a 0-filled column to the matrix
    void addColumn();

    // set a value in position (pos1 - pos2) of the matrix
    void setValue(int pos1, int pos2, float value);

    // returns the value in position (pos1 - pos2) of the matrix
    float getValue(int pos1, int pos2);

    // output stream operator overloading, for easily reading the dynamic matrix
    friend ostream& operator<<(ostream& out, const DynamicMatrix& stm);

    void setDimension(int nRows, int nColumns);
};

#endif //HMM4GA_CLION_DYNAMICMATRIX_H