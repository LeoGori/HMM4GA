#include "DynamicMatrix.h"
#include<bits/stdc++.h>

DynamicMatrix::DynamicMatrix(float initValue) {
    initializationValue = initValue;

    vector<float> cell;
    cell.push_back(0);
    matrix.push_back(cell);
}

// adds a row filled by 0s
void DynamicMatrix::addRow() {

    vector<float> values;
    values.reserve(matrix.size());

    for(int i=0; i<matrix[0].size(); i++)
        values.push_back(0);

    matrix.push_back(values);
}

// add a column filled by 0s
void DynamicMatrix::addColumn() {
    for(int i=0; i<matrix.size(); i++)
        matrix[i].push_back(0);
}

float DynamicMatrix::getInitializationValue() const {
    return initializationValue;
}

void DynamicMatrix::setInitializationValue(float initValue) {
    initializationValue = initValue;
}

float DynamicMatrix::getValue(int pos1, int pos2) {

    if(pos1 > 0 && pos2 > 0)
        return this->getProbability(pos1, pos2);
    else
        return initializationValue;
}


ostream& operator<<(ostream& out, const DynamicMatrix& dm) {

    out << "    |";
    for (int i = 0; i < dm.matrix[0].size(); i++) {
        out << "   " << i << "   |";
    }
    out << endl;

    for(int i=0; i < dm.matrix.size(); i++){
        out << "|";
        out << " " << i << " |";
        for(float value : dm.matrix[i]){

            out << " " << fixed << setprecision(3) <<  value << " |";
        }
        out << "\n";
    }
    return out;
}

void DynamicMatrix::setValue(int pos1, int pos2, float value) {
    matrix[pos1][pos2] = value;
}
