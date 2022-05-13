#include<string>


class SequenceGenerator {
private:
    char nucleotide[4] = {'A', 'C', 'G', 'T'};
public:
    SequenceGenerator();
    char getSymbol();
};