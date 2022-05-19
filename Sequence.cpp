#include "Sequence.h"
#include <iostream>
#include <utility>

using namespace std;

Sequence::Sequence(char symbol) {
    data = symbol;
}

Sequence::Sequence(string sequence) {
    data = std::move(sequence);
}

void Sequence::addResidue(char residue) {
    data = data + residue;
}

void Sequence::printSequence() {
    cout << data << endl;
}

int Sequence::getLength() {
    return data.length();
}

char Sequence::getChar(int pos) {
    if (pos >= 0 && pos < data.length())
        return data[pos];
    return 0;
}

void Sequence::getNext() {
    data += sg.getSymbol();
}

string Sequence::getSequence() {
    return data;
}



// this method returns the emission probability of the input nucleotide
// Since the emission probability is computed exclusively through the phred quality score,
// as defined in the reference document, this method has been commented. The emission probability is
// computed inside the PairHMM class, inside getPerBaseEmission(...) method

/*float Sequence::getEmissionProbability(char nucleotide) {
    return sg.getEmissionProbability(nucleotide);
}*/

