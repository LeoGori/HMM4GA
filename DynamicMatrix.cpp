#include "DynamicMatrix.h"
#include<bits/stdc++.h>

/*// adds a row filled by 0s
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

void DynamicMatrix::setValue(int pos1, int pos2, float value) {
    matrix[pos1 + 1][pos2 + 1] = value;
}

void DynamicMatrix::setDimension(int readLength, int haplotypeLength) {

    vector<float> values;
    for(int i = 0; i<readLength + 1; i++) {
        for (int j = 0; j < haplotypeLength + 1; j++)
            values.push_back(0);
        matrix.push_back(values);
        values.clear();
    }

}
