#include "Sequence.h"
#include <iostream>
#include "omp.h"
#include "SequenceGenerator.h"

using namespace std;


int main() {

    SequenceGenerator seqGen;
    Sequence seq(seqGen.getSymbol());

    int N = seq.getLength();

    cout << N << endl;

    seq.printSequence();

    omp_set_num_threads(3);

    for(int i = 0; i < 5; i++) {

        seq.addResidue(seqGen.getSymbol());
        seq.printSequence();

        #pragma omp parallel for
            for (int j = 0; j < seq.getLength(); j++) {
                #pragma omp critical
                {
                    int ID = omp_get_thread_num();
                    cout << "Lettera " << seq.getChar(j) << ", Thread ID: " << ID << " relativa al ciclo numero " << i
                         << endl;
                }
            }

    }
}