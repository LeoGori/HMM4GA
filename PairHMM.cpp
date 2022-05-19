#include "PairHMM.h"
#include "omp.h"
#include<iostream>
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
    // all the matrices are completely 0-filled in this step
    matchMatrix.setDimension(read.getLength(), haplotype.getLength());

    insertionMatrix.setDimension(read.getLength(), haplotype.getLength());

    deletionMatrix.setDimension(read.getLength(), haplotype.getLength());


    // following convention described in "references/PHMMFA_Ren_Bertels_Al-Ars.pdf"
    // the other inizialization values in the matrix
    // are already initialized, since the matrices M, I and D are completely 0-filled by inizialization
    float value = (float) 1 / (float) haplotype.getLength();


    for(int i=0; i< haplotype.getLength(); i++) {
        deletionMatrix.setValue(-1, i, value);
    }

    // following convention described in "Biological Sequence Analysis: Probabilistic Models of Proteins and Nucleic Acids", page 88 (chapter 4.2),
    /*// set the value of the element in position (-1,-1) of matrix M to 1, while the other inizialization values in the matrix
    // are already initialized, since the matrices M, I and D are completely 0-filled by inizialization
    for(int i=0; i< haplotype.getLength(); i++) {

        matchMatrix.setValue(-1, -1, 1);
    }*/

    cout << "The state transition matrix is:\n" << transitionMatrix << endl;

}

void PairHMM::run() {

    float prior;
    float value;

    cout << "The haplotype is: " << haplotype.getSequence() << endl;
    cout << "The read is: " + read.getSequence() << endl;

    for(int i = 0; i < read.getLength(); i++) {

#pragma omp parallel for default(none) shared(haplotype, i, cout) private(prior, value)
        for (int j = 0; j < haplotype.getLength(); j++) {
#pragma omp critical
            {

                // define the prior probability, based on the value of the base quality score
                prior = getPerBaseEmission(read.getChar(i), haplotype.getChar(j));

                // compute the value of the element (i, j) of the matrix M
                value = prior *
                              (matchMatrix.getValue(i-1, j-1) * transitionMatrix.getProbability('M', 'M') +
                               insertionMatrix.getValue(i-1, j-1) * transitionMatrix.getProbability('I', 'M') +
                               deletionMatrix.getValue(i-1, j-1) * transitionMatrix.getProbability('D', 'M'));

                // assign the value to the element (i, j) of the matrix M
                matchMatrix.setValue(i,j, value);

                // compute the value of the element (i, j) of the matrix I
                value = matchMatrix.getValue(i-1, j) * transitionMatrix.getProbability('M', 'I') +
                        insertionMatrix.getValue(i-1, j) * transitionMatrix.getProbability('I', 'I');

                // assign the value to the element (i, j) of the matrix I
                insertionMatrix.setValue(i, j, value);

                // compute the value of the element (i, j) of the matrix D
                value = matchMatrix.getValue(i, j-1) * transitionMatrix.getProbability('M', 'D') +
                        deletionMatrix.getValue(i, j-1) * transitionMatrix.getProbability('D', 'D');

                // assign the value to the element (i, j) of the matrix D
                deletionMatrix.setValue(i, j, value);

                int ID = omp_get_thread_num();
                cout << "Analyzing chars " << read.getChar(i) << " and "<< haplotype.getChar(j) << ", executed by ThreadID: " << ID << ", during cycle i=" << i << ", j=" << j << endl;

                cout << "Matrice M\n\n" << matchMatrix << endl;

                cout << "Matrice I\n\n" << insertionMatrix << endl;

                cout << "Matrice D\n\n" << deletionMatrix << endl;


            }
        }

    }

    // result of the algorithm according to the book:
    // “Biological Sequence Analysis: Probabilistic Models of Proteins and Nucleic Acids”, page 88, chapter 4.2
    /*float result = matchMatrix.getValue(read.getLength()-1, haplotype.getLength()-1) +
                   insertionMatrix.getValue(read.getLength()-1, haplotype.getLength()-1);*/


    // result took according to both the reference documents
    float result = 0;

    for(int j=0;j<haplotype.getLength(); j++)
        result += matchMatrix.getValue(read.getLength()-1, j) +
            insertionMatrix.getValue(read.getLength()-1, j);


    cout << "the probability P(R|H) of read R (" << read.getSequence() <<") to be sequenced from haplotype H (" << haplotype.getSequence() << ") is: " << result << endl;

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




