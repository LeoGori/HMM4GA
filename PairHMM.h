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

    // the base quality score
    // possible improvement: the quality may be an attribute of the class Sequence (in particular of the read)
    // and should be different for each base and its position in the sequence (?)
    int baseQuality;



public:
    PairHMM();

    // runs the algorithm described in the reference document
    void run();

    // returns the emission probability, based on the base quality score and if the characters match or not
    float getPerBaseEmission(char readChar, char haplotypeChar) const;

    // returns the errorRate function, according to the value of the base quality score
    float getErrorRate() const;
};

