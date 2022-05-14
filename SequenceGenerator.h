#include<string>
#include<map>
#include <vector>

//enum Nucleotide {A = 0, C = 1, G = 2, T = 3};

using namespace std;

class SequenceGenerator {
private:
    map<char, int> odds;
    vector<char> nucleotide;
    map<char, int> nucleotides;
    int totalOdds;
    vector<char> priorDistribution;
    vector<float> emissionProbability;


public:
    SequenceGenerator();
    char getSymbol();
    float getEmissionProbability(char symbol);
};