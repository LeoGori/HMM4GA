#include "PairHMM.h"
#include "omp.h"
#include <iostream>
#include <cmath>

using namespace std;


PairHMM::PairHMM() {

    omp_set_num_threads(8);

    baseQuality = 10;

    transitionMatrix.populateMatrix();

    for(int i=0;i<5; i++)
        haplotype.getNext();

    for(int i=0;i<3; i++)
        read.getNext();

    // set the dimension of M, I and D matrices to (n x m), where n is the length of the read, and m is the length of the haplotype
    // all the matrices are completely filled by -1 values in this step
    matchMatrix.setDimension(read.getLength(), haplotype.getLength());

    insertionMatrix.setDimension(read.getLength(), haplotype.getLength());

    deletionMatrix.setDimension(read.getLength(), haplotype.getLength());


    // following convention described in "references/PHMMFA_Ren_Bertels_Al-Ars.pdf"

    for(int i=0; i< read.getLength(); i++) {
        deletionMatrix.setValue(i, -1, 0);
        matchMatrix.setValue(i, -1, 0);
        insertionMatrix.setValue(i, -1, 0);
    }

    float value = (float) 1 / (float) haplotype.getLength();
    for(int j=-1; j< haplotype.getLength(); j++) {
        deletionMatrix.setValue(-1, j, value);
        matchMatrix.setValue(-1, j, 0);
        insertionMatrix.setValue(-1, j, 0);
    }

    cout << " The state transition matrix is:\n" << transitionMatrix << endl;
    cout << "###########################################################################################################################################################" << endl;

}

void PairHMM::run() {

    cout << "Matrice M\n\n" << matchMatrix << endl;

    cout << "Matrice I\n\n" << insertionMatrix << endl;

    cout << "Matrice D\n\n" << deletionMatrix << endl;

    cout
            << "###########################################################################################################################################################"
            << endl;

    cout << "The haplotype is: " << haplotype.getSequence() << endl;
    cout << "The read is: " + read.getSequence() << endl;

    cout
            << "###########################################################################################################################################################"
            << endl;



/* access to the upper left triangle of the matrix
 * A = accessed
 * N = Not accessed
 *
 * Example with matrix M(3x5)
 *
 * | A | A | A | N | N |
 * | A | A | N | N | N |
 * | A | N | N | N | N |
 *
 */
    for (int k = 0; k < read.getLength(); k++) {

#pragma omp parallel for default(none) shared(k)
        for (int i = k; i >= 0; i--) {
            execMatricesComputation(i, k - i);
        }

    }

/* access to the inner part of the matrix (just in case the matrix is not squared)
 * A = accessed
 * N = Not accessed
 *
 * Example with matrix M(3x5)
 *
 * | N | N | N | A | A |
 * | N | N | A | A | N |
 * | N | A | A | N | N |
 *
 */
    for (int k = 1; k <= haplotype.getLength() - read.getLength(); k++) {
        cout << "inizio ciclo" << endl;
#pragma omp parallel for default(none) shared(k, cout, read)
        for (int i = read.getLength() - 1; i >= 0; i--) {
            execMatricesComputation(i, k + read.getLength() - i - 1);
        }
        cout << "fine ciclo" << endl;
    }

/* access to the lower triangle of the matrix
 * A = accessed
 * N = Not accessed
 *
 * Example with matrix M(3x5)
 *
 * | N | N | N | N | N |
 * | N | N | N | N | A |
 * | N | N | N | A | A |
 *
 */
    for (int k = haplotype.getLength() - read.getLength() + 1; k<haplotype.getLength(); k++) {

#pragma omp parallel for default(none) shared(haplotype, read, cout, k)
        for (int j = k; j < haplotype.getLength(); j++) {
            execMatricesComputation( read.getLength() - (j - k) -1 , j);
        }

    }


    // result of the algorithm according to the book:
    // “Biological Sequence Analysis: Probabilistic Models of Proteins and Nucleic Acids”, page 88, chapter 4.2
    /*float result = matchMatrix.getValue(read.getLength()-1, haplotype.getLength()-1) +
                   insertionMatrix.getValue(read.getLength()-1, haplotype.getLength()-1);*/


    // result took according to both the reference documents
    float result = 0;

    for (int j = 0; j < haplotype.getLength(); j++)
        result += matchMatrix.getValue(read.getLength() - 1, j) +
                  insertionMatrix.getValue(read.getLength() - 1, j);


    cout << "the probability P(R|H) of read R (" << read.getSequence() << ") to be sequenced from haplotype H ("
         << haplotype.getSequence() << ") is: " << result << endl;

}


float PairHMM::getPerBaseEmission(char readChar, char haplotypeChar) const {
    float errorRate = getErrorRate();

    if(readChar == haplotypeChar)
        return 1 - errorRate;
    return errorRate/3;
}


float PairHMM::getErrorRate() const {
    float ratio = (float) baseQuality / 10;

    auto errorRate = (float) pow(0.1, ratio);

    return errorRate;
}

void PairHMM::execMatricesComputation(int i, int j) {
    float prior;
    float value;
    // define the prior probability, based on the value of the base quality score
    prior = getPerBaseEmission(read.getChar(i), haplotype.getChar(j));

    // compute the value of the element (i, j) of the matrix M
    value = prior *
            (matchMatrix.getValue(i - 1, j - 1) * transitionMatrix.getProbability('M', 'M') +
             insertionMatrix.getValue(i - 1, j - 1) * transitionMatrix.getProbability('I', 'M') +
             deletionMatrix.getValue(i - 1, j - 1) * transitionMatrix.getProbability('D', 'M'));

    // assign the value to the element (i, j) of the matrix M
    matchMatrix.setValue(i, j, value);

    // compute the value of the element (i, j) of the matrix I
    value = matchMatrix.getValue(i - 1, j) * transitionMatrix.getProbability('M', 'I') +
            insertionMatrix.getValue(i - 1, j) * transitionMatrix.getProbability('I', 'I');

    // assign the value to the element (i, j) of the matrix I
    insertionMatrix.setValue(i, j, value);

    // compute the value of the element (i, j) of the matrix D
    value = matchMatrix.getValue(i, j - 1) * transitionMatrix.getProbability('M', 'D') +
            deletionMatrix.getValue(i, j - 1) * transitionMatrix.getProbability('D', 'D');

    // assign the value to the element (i, j) of the matrix D
    deletionMatrix.setValue(i, j, value);

    int ID = omp_get_thread_num();

#pragma omp critical
    {

        cout << "Analyzing chars " << read.getChar(i) << " and " << haplotype.getChar(j)
             << ", executed by ThreadID: " << ID << ", during cycle i=" << i << ", j=" << j << endl;

        cout << "Matrice M\n\n" << matchMatrix << endl;

        cout << "Matrice I\n\n" << insertionMatrix << endl;

        cout << "Matrice D\n\n" << deletionMatrix << endl;

        cout
                << "###########################################################################################################################################################"
                << endl;
    }
}





