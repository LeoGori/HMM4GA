#include <vector>
#include <map>
#include <string>


using namespace std;


class StateTransitionMatrix {

private:
    vector<vector<float>> transitionMatrix;
    map<char, int> nucleotides;

public:
    StateTransitionMatrix();

    float getTransitionProbability(char char1, char char2);

    vector<vector<float>> getTransitionMatrix();

    void populateMatrix();

//    void print();

    friend ostream& operator<< (std::ostream& stream, const StateTransitionMatrix& matrix);

};