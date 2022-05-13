using namespace std;
#include<string>

class Sequence {
private:
    string data;
public:
    Sequence(char symbol);
    Sequence(string sequence="");
    void addResidue(char residue);
    void printSequence();
    /* string get_list();*/
    char getChar(int pos);
    int getLength();
};