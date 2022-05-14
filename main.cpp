#include "Sequence.h"
#include <iostream>
#include "omp.h"
#include "SequenceGenerator.h"
#include "StateTransitionMatrix.h"

using namespace std;


int main() {

    StateTransitionMatrix T;
    SequenceGenerator haplotypeSeqGen;
    Sequence haplotypeSeq;

    T.populateMatrix();

    cout << T << endl;

    omp_set_num_threads(3);
    for(int i = 0; i < 5; i++) {

        haplotypeSeq.addResidue(haplotypeSeqGen.getSymbol());
        haplotypeSeq.printSequence();

        #pragma omp parallel for default(none) shared(haplotypeSeq, i, cout)
            for (int j = 0; j < haplotypeSeq.getLength(); j++) {
                #pragma omp critical
                {
                    int ID = omp_get_thread_num();
                    cout << "Lettera " << haplotypeSeq.getChar(j) << ", Thread ID: " << ID << " relativa al ciclo numero " << i << endl;
                }
            }

    }
}