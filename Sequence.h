#include "SequenceGenerator.h"
#include<string>

/*
 * Sequence of nucleotides, defined casually through a random emission distribution defined by SequenceGenerator
 * that act as a private attribute in this class, and that is initialized in construction of the sequence entity
 */

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

    // returns the whole sequence
    string getSequence();

    // returns the emission probability of the input nucleotide
    //float getEmissionProbability(char nucleotide);
};