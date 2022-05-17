#include "PairHMM.h"
#include "omp.h"
#include<iostream>
#include <cmath>

PairHMM::PairHMM(): matchMatrix(1), insertionMatrix(0), deletionMatrix(0) {

    omp_set_num_threads(3);

    transitionMatrix.populateMatrix();

    for(int i=0;i<5; i++)
        haplotype.getNext();

    for(int i=0;i<3; i++)
        read.getNext();

    matchMatrix.setDimension(read.getLength(), haplotype.getLength());

    insertionMatrix.setDimension(read.getLength(), haplotype.getLength());

    deletionMatrix.setDimension(read.getLength(), haplotype.getLength());

}

void PairHMM::run() {

    float prior;
    float value;

    read.printSequence();
    haplotype.printSequence();

    for(int i = 0; i < read.getLength(); i++) {

#pragma omp parallel for default(none) shared(haplotype, i, cout) private(prior, value)
        for (int j = 0; j < haplotype.getLength(); j++) {
#pragma omp critical
            {
                prior = getPerBaseEmission(read.getChar(i), haplotype.getChar(j));

                value = prior *
                              (matchMatrix.getValue(i-1, j-1) * transitionMatrix.getProbability('M', 'M') +
                               insertionMatrix.getValue(i-1, j-1) * transitionMatrix.getProbability('I', 'M') +
                               deletionMatrix.getValue(i-1, j-1) * transitionMatrix.getProbability('D', 'M'));

                matchMatrix.setValue(i,j, value);

                value = matchMatrix.getValue(i-1, j) * transitionMatrix.getProbability('M', 'I') +
                        insertionMatrix.getValue(i-1, j) * transitionMatrix.getProbability('I', 'I');

                insertionMatrix.setValue(i, j, value);

                value = matchMatrix.getValue(i, j-1) * transitionMatrix.getProbability('M', 'D') +
                        deletionMatrix.getValue(i, j-1) * transitionMatrix.getProbability('D', 'D');

                deletionMatrix.setValue(i, j, value);

                int ID = omp_get_thread_num();
                cout << "Lettera " << read.getChar(i) << " e "<< haplotype.getChar(j) << ", Thread ID: " << ID << " relativa al ciclo numero " << i << endl;

                cout << "Matrice M\n\n" << matchMatrix << "\n, Thread ID: " << ID << " relativa al ciclo numero " << i << endl;

                cout << "Matrice I\n\n" << insertionMatrix << "\n, Thread ID: " << ID << " relativa al ciclo numero " << i << endl;

                cout << "Matrice D\n\n" << deletionMatrix << "\n, Thread ID: " << ID << " relativa al ciclo numero " << i << endl;


            }
        }

    }


}

float PairHMM::getPerBaseEmission(char read, char haplotype) {
    if(read == haplotype)
        return 1 - getErrorRate();
    return getErrorRate()/3;
}

float PairHMM::getErrorRate() {
    float errorRate = 10 ^ (-baseQuality/10);

    return errorRate;
}




