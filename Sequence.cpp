#include "Sequence.h"
#include <iostream>

Sequence::Sequence(char symbol) {
    data = symbol;
}

Sequence::Sequence(string sequence) {
    data = sequence;
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
}