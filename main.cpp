#include "Sequence.h"
#include <iostream>
#include "omp.h"
#include "SequenceGenerator.h"
#include "StateTransitionMatrix.h"

using namespace std;


int main() {

    StateTransitionMatrix T;
    SequenceGenerator seqGen;
    Sequence seq(seqGen.getSymbol());

    seq.printSequence();

    T.populateMatrix();

    char val[4] = {'A', 'C', 'G', 'T'};

    /*for(char i : val) {
        for (char j : val)
            cout << T.getTransitionProbability(i, j) << " ";
        cout << endl;

    }*/

    cout << T << endl;

    omp_set_num_threads(3);
    for(int i = 0; i < 5; i++) {

        seq.addResidue(seqGen.getSymbol());
        seq.printSequence();

        #pragma omp parallel for default(none) shared(seq, i, cout)
            for (int j = 0; j < seq.getLength(); j++) {
                #pragma omp critical
                {
                    int ID = omp_get_thread_num();
                    cout << "Lettera " << seq.getChar(j) << ", Thread ID: " << ID << " relativa al ciclo numero " << i << endl;
                }
            }

    }
}