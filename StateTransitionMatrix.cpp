#include "StateTransitionMatrix.h"
#include<iostream>
#include<random>
#include<bits/stdc++.h>

using namespace std;


StateTransitionMatrix::StateTransitionMatrix() {

    nucleotides = {
            { 'A', 0 },
            { 'C', 1 },
            { 'G', 2 },
            { 'T', 3 },
    };

}

float StateTransitionMatrix::getTransitionProbability(char char1, char char2) {

    float probability = transitionMatrix[nucleotides[char1]][nucleotides[char2]];

    return probability;
}

vector<vector<float>> StateTransitionMatrix::getTransitionMatrix() {
    return transitionMatrix;
}

void StateTransitionMatrix::populateMatrix() {

    random_device rd;   // non-deterministic generator
    mt19937 gen(rd());  // to seed mersenne twister.
    uniform_int_distribution<int> dist(1,1000); // distribute results between 1 and 1000 inclusive.


    vector<float> values;

    for(int i=0; i<nucleotides.size(); i++){

        int sum = 0;
        for (int j=0; j< nucleotides.size(); j++) {

            int value = dist(gen);

            sum += value;
            values.push_back(static_cast<float>(value));

        }
        for (int j = 0; j < nucleotides.size(); j++) {
            values[j] /= (float) sum;
        }
        transitionMatrix.push_back(values);
        values.clear();
    }
}

/*void StateTransitionMatrix::print(){
    for(auto & iterator : transitionMatrix)
        for (int j=0; j < transitionMatrix.size(); j++)
            cout << iterator[j] << endl;
}*/

ostream& operator<<(ostream& out, const StateTransitionMatrix& stm) {

    char nuc[4] = {'A', 'C', 'G', 'T'};
    out << "    |";
    for (int i = 0; i < stm.nucleotides.size(); i++) {
        out << "   " << nuc[i] << "   |";
    }
    out << endl;

    for(int i=0; i < stm.nucleotides.size(); i++){
        out << "|";
        out << " " << nuc[i] << " |";
        for(int j=0; j < stm.nucleotides.size(); j++){

            float value = stm.transitionMatrix[i][j];
            out << " " << fixed << setprecision(3) <<  value << " |";
        }
        out << "\n";
    }
    return out;
}
