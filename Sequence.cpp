#include "Sequence.h"
#include <iostream>
#include <utility>

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
    if (pos > 0 && pos < data.length())
        return data[pos];
    return 0;
}

void Sequence::getNext() {
    data += sg.getSymbol();
}

float Sequence::getEmissionProbability(char nucleotide) {
    return sg.getEmissionProbability(nucleotide);
}

