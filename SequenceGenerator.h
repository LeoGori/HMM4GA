#include<string>
#include<map>
#include <vector>

using namespace std;

/* The class that is responsible for defining a random distribution of emission probabilities of
 * a sequence of nucleotides
 * */

class SequenceGenerator {
private:

    // in order to simulate a prior distribution, it saves the number of occurrences 'previously' seen for each nucleotide
    map<char, int> odds;

    // int to char relationship that is useful when iterating over the nucleotides
    vector<char> nucleotide;

    // char to int relationship that is useful when iterating over the nucleotides
    map<char, int> nucleotides;

    // total number of occurrences of the nucleotide in the simulated sequence - i.e. the sum of all the occurrences of the nucleotides
    int totalOdds;

    // the actual simulated sequence, it is used in order to casually retrieve the nucleotide for the casual generation
    // it is used for the implementation of the method getNext() of class Sequence
    string fakeSequence;

    // the actual emission distribution
    vector<float> emissionProbability;

    // to be defined (?); the phred quality score introduced in the reference paper
    vector<int> baseQuality;

public:

    SequenceGenerator();

    // returns a casual nucleotide, according to the casual emission distribution
    char getSymbol();

    // returns the emission probability of the input nucleotide
    //float getEmissionProbability(char symbol);

    //void setEmissionProbability(map<char, float> emissionDistribution);
};