#include<string>

//enum Nucleotide {A = 0, C = 1, G = 2, T = 3};

class SequenceGenerator {
private:
    char nucleotide[4] = {'A', 'C', 'G', 'T'};
public:
    SequenceGenerator();
    char getSymbol();
};