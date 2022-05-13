#include "SequenceGenerator.h"
#include "time.h"

SequenceGenerator::SequenceGenerator() {
    srand(time(NULL));
}

char SequenceGenerator::getSymbol() {

    int index = rand() % sizeof(nucleotide);
    return nucleotide[index];
}