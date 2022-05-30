#include "DynamicMatrix.h"
#include<bits/stdc++.h>

/*// adds a row filled by 0s
 * this method can be useful in the case a read sequence is generated at runtime, inside the execution of the
 * hmm pair forward algorithm, for instance in the case a new characters of the read needs to be generated inside the
 * for loop of the algorithm.
 * In such a case, a new row needs to be added to the M, I and D matrices.
 * In the current implementation we consider the haplotype and read sequences already defined before the execution of the algorithm
 * * For this reason this method has been commented
void DynamicMatrix::addRow() {

    vector<float> values;
    values.reserve(matrix.size());

    if(!matrix[0].empty()) {
        for (int i = 0; i < matrix[0].size(); i++)
            values.push_back(0);
    }
    else
        values.push_back(0);

    matrix.push_back(values);
}

// add a column filled by 0s
  * this method can be useful in the case a haplotype sequence is generated at runtime, inside the execution of the
 * hmm pair forward algorithm, for instance in the case a new character of the haplotype needs to be generated inside the
 * for loop of the algorithm.
 * In such a case, a new column needs to be added to the M, I and D matrices.
 * In the current implementation we consider the haplotype and read sequences already defined before the execution of the algorithm
 * For this reason this method has been commented
void DynamicMatrix::addColumn() {
    if(!matrix.empty()) {
        for (int i = 0; i < matrix.size(); i++)
            matrix[i].push_back(0);
    }
    else {
        vector<float> values;
        values.push_back(0);
        matrix.push_back(values);
    }
}*/


// returns the i+1, j+1 value of the matrix, hiding the user the fact that elements:
// - [0,0];
// - [0,j] for 1 <= j < length of the haplotype sequence;
// - [i,0] for 1 <= i < length of the read sequence
// are all reserved for initialization purposes, and do not logically belong to the matrix
float DynamicMatrix::getValue(int pos1, int pos2) {

    return this->getProbability(pos1 + 1, pos2 + 1);
}


ostream& operator<<(ostream& out, const DynamicMatrix& dm) {

    out << "      |";
    for (int i = 0; i < dm.matrix[0].size(); i++) {
        if(i>0)
            out<< " ";
        out << "     " << i-1 << "       |";
    }
    out << endl;

    for(int i=0; i < dm.matrix.size(); i++){
        out << "|";
        if(i>0)
            out<< " ";
        out << " " << i-1 << "  |";
        for(float value : dm.matrix[i]){

            out << " " << fixed << setprecision(10) <<  value << " |";
        }
        out << "\n";
    }
    return out;
}


// sets the i+1, j+1 value of the matrix, hiding the user the fact that elements:
// - [0,0];
// - [0,j] for 1 <= j < length of the haplotype sequence;
// - [i,0] for 1 <= i < length of the read sequence
// are all reserved for initialization purposes, and do not logically belong to the matrix
void DynamicMatrix::setValue(int pos1, int pos2, float value) {
    matrix[pos1 + 1][pos2 + 1] = value;
}

/*
 * Sets the dimension of the matrix such that:
 * - number of rows = length of the read + 1;
 * - number of columns = length of the haplotype + 1
 * the additional column and row are reserved for initialization purposes
 */
void DynamicMatrix::setDimension(int readLength, int haplotypeLength) {

    vector<float> values;
    for(int i = 0; i<readLength + 1; i++) {
        for (int j = 0; j < haplotypeLength + 1; j++)
            values.push_back(-1);
        matrix.push_back(values);
        values.clear();
    }

}
