#include "StateTransitionMatrix.h"
#include<iostream>
#include<random>
#include<bits/stdc++.h>




StateTransitionMatrix::StateTransitionMatrix() {

//    states = {
//            { 'A', 0 },
//            { 'C', 1 },
//            { 'G', 2 },
//            { 'T', 3 },
//    };

    states = {
            { 'M', 0 },
            { 'I', 1 },
            { 'D', 2 },
    };

}

float StateTransitionMatrix::getTransitionProbability(char char1, char char2) {

    float probability = transitionMatrix[states[char1]][states[char2]];

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

    for(int i=0; i < states.size(); i++){

        int sum = 0;
        for (int j=0; j < states.size(); j++) {

            int value = dist(gen);

            sum += value;
            values.push_back(static_cast<float>(value));

        }
        for (int j = 0; j < states.size(); j++) {
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

    char state[3] = {'M', 'I', 'D'};
    out << "    |";
    for (int i = 0; i < stm.states.size(); i++) {
        out << "   " << state[i] << "   |";
    }
    out << endl;

    for(int i=0; i < stm.states.size(); i++){
        out << "|";
        out << " " << state[i] << " |";
        for(int j=0; j < stm.states.size(); j++){

            float value = stm.transitionMatrix[i][j];
            out << " " << fixed << setprecision(3) <<  value << " |";
        }
        out << "\n";
    }
    return out;
}
