#include "SequenceGenerator.h"
#include <random>
#include <algorithm>
#include <iostream>


SequenceGenerator::SequenceGenerator() {

    nucleotide.push_back('A');
    nucleotide.push_back('C');
    nucleotide.push_back('G');
    nucleotide.push_back('T');

    for(int i=0;i<nucleotide.size(); i++)
        nucleotides.insert({nucleotide[i], i});

    totalOdds = 0;
    random_device rd;   // non-deterministic generator
    mt19937 gen(rd());  // to seed mersenne twister.
    uniform_int_distribution<int> dist(1,100); // distribute results between 1 and 1000 inclusive.

    for(char & i : nucleotide) {
        int value = dist(gen);

        totalOdds += value;
        odds.insert({i, value});

        for(int j=0; j<value; j++)
            fakeSequence.push_back(i);
    }

    char nuc;
    float num;
    auto denum = static_cast<float>(totalOdds);

    for(char i : nucleotide){

        nuc = i;

        num = static_cast<float>(odds.find(nuc)->second);

        emissionProbability.push_back(num/denum);

    }

    auto rng = default_random_engine {rd()};
    std::shuffle(std::begin(fakeSequence), std::end(fakeSequence), rng);

}

char SequenceGenerator::getSymbol() {

    random_device rd;   // non-deterministic generator
    mt19937 gen(rd());  // to seed mersenne twister.
    uniform_int_distribution<int> dist(0,totalOdds - 1); // distribute results between 1 and 1000 inclusive.

    int index = dist(gen);

    return fakeSequence[index];
}

float SequenceGenerator::getEmissionProbability(char base) {

    return emissionProbability[nucleotides.find(base)->second];
}

