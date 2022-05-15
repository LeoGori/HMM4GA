#include "SequenceGenerator.h"
#include<string>

using namespace std;

class Sequence {
private:
    string data;
    SequenceGenerator sg;
public:

    explicit Sequence(char symbol);
    explicit Sequence(string sequence="");

    // adds an user defined nucleotide to the bottom of the sequence
    void addResidue(char residue);

    // prints the sequence
    void printSequence();

    // returns the char of the sequence that is in position pos
    char getChar(int pos);

    // returns the length of the sequence
    int getLength();

    // adds a random nucleotide to the sequence
    // the emission probability is defined randomly by SequenceGenerator
    void getNext();

    // returns the emission probability of the input nucleotide
    float getEmissionProbability(char nucleotide);
};