#include "PairHMM.h"
#include "omp.h"
#include<iostream>
#include <cmath>

PairHMM::PairHMM(): matchMatrix(0), insertionMatrix(0), deletionMatrix(INFINITY) {

    omp_set_num_threads(3);

    transitionMatrix.populateMatrix();
}

void PairHMM::run() {

    for(int i = 0; i < 5; i++) {

        haplotype.getNext();

        haplotype.printSequence();

        matchMatrix.addRow();
        matchMatrix.addColumn();
        matchMatrix.addRow();

        matchMatrix.setValue(0,0, 1.2);

        insertionMatrix.addColumn();
        insertionMatrix.addRow();

        insertionMatrix.setValue(1,0, 1.2);

        deletionMatrix.addColumn();
        deletionMatrix.addRow();

        deletionMatrix.setValue(0,1, 1.2);

        cout << matchMatrix << endl;

        cout << insertionMatrix << endl;

        cout << deletionMatrix << endl;

#pragma omp parallel for default(none) shared(haplotype, i, cout)
        for (int j = 0; j < haplotype.getLength(); j++) {
#pragma omp critical
            {
                int ID = omp_get_thread_num();
                cout << "Lettera " << haplotype.getChar(j) << ", Thread ID: " << ID << " relativa al ciclo numero " << i << endl;
            }
        }

    }
}


