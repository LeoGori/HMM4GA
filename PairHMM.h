#include "StateTransitionMatrix.h"
#include "DynamicMatrix.h"
#include "Sequence.h"

class PairHMM {
private:

    // haplotype's sequence
    Sequence haplotype;

    // read's sequence
    Sequence read;

    // transition matrix T
    StateTransitionMatrix transitionMatrix;

    // matching matrix M
    DynamicMatrix matchMatrix;

    // insetion matrix I
    DynamicMatrix insertionMatrix;

    // deletion matrix M
    DynamicMatrix deletionMatrix;


public:
    PairHMM();

    // TODO
    // runs the algorithm described in the reference document
    void run();
};

