# Hidden Markov Models for Genome Analysis
The project's goal is the development of a basic implementation of the **pair hidden Markov Model (HMM) forward algorithm** for genomic sequence analysis (described in the
[**reference document**](https://github.com/LeoGori/HMM4GA/blob/main/references/PHMMFA_Benjamin.pdf)), with the introduction of parallel
computation through the use of [**OpenMP**](https://www.openmp.org/) APIs.

Further details are provided in these articles ([**[1]**](https://github.com/LeoGori/HMM4GA/blob/main/references/PHMMFA_Ren_Bertels_Al-Ars.pdf), [**[2]**](https://github.com/LeoGori/HMM4GA/blob/main/references/2017_fpga_pairhmm.pdf)) and on the book **Biological Sequence Analysis: Probabilistic Models of Proteins and Nucleic Acids** (page 88, §4.2).

# Code Description
The main files that build up the project are:

- [**_Sequence.h_**](https://github.com/LeoGori/HMM4GA/blob/main/Sequence.h): class that represents the sequence of nucleotides, it contains the string of characters that compose the sequence and the class SequenceGenerator
- [**_SequenceGenerator.h_**](https://github.com/LeoGori/HMM4GA/blob/main/SequenceGenerator.h): class that defines a random emission probability distribution of a sequence of nucleotides. Currently, an instance of the class Sequence is randomly generated according to its SequenceGenerator.
- [**_ProbabilityMatrix.h_**](https://github.com/LeoGori/HMM4GA/blob/main/ProbabilityMatrix): class that represents a generic matrix of floating point values, from which the classes **DynamicMatrix** and **StateTransitionMatrix** inherit common attributes and methods. DynamicMatrix adds the possibility of adding rows and columns dynamically, while StateTransitionMatrix provides a series of states, and a mapping between them and the indexes of the matrix
- [**_PairHMM.h_**](https://github.com/LeoGori/HMM4GA/blob/main/PairHMM.h): the class that implements the pair HMM forwarding algorithm, it encloses 2 instances of the class **Sequence** (one for defining the read sequence, and one for defining the haplotype sequence), 1 instance of the class **StateTransitionMatrix** (for defining matrix T), and 3 instances of the class **DynamicMatrix** (for the definition of matrices M, I and D)
- [**_main.cpp_**](https://github.com/LeoGori/HMM4GA/blob/main/main.cpp): the entry point of the program, contains an instance of the class **PairHMM** and the call of its method for the execution of the PairHMM forwarding algorithm

# Language and APIs
The code is entirely written in C++ programming language, with the use of the following libraries and APIs (omitting the standard ones):

- [**_random_**](https://www.cplusplus.com/reference/random/): used for the random generation of sequences and the random definition of state transition probabilities
- [**_algorithm_**](https://www.cplusplus.com/reference/algorithm/): used for the shuflling of sequences, used for randomization purposes
- [**_OpenMP_**](https://www.openmp.org/): used in [**_PairHMM.cpp_**](https://github.com/LeoGori/HMM4GA/blob/main/PairHMM.cpp) for introducing thread level computation in the algorithm


# How to run the code (Windows)

1. Install [MinGw64](https://winlibs.com/) version > 9.2 (otherwise the random generated sequence will be the same at each execution, as reported [here](https://stackoverflow.com/questions/18880654/why-do-i-get-the-same-sequence-for-every-run-with-stdrandom-device-with-mingw) and [here](https://en.cppreference.com/w/cpp/numeric/random/random_device))
2. Install [CMake](https://cmake.org/download/)
3. Create folder for building project
```
  mkdir build
  cd build
```
4.Generate the makefiles
```
  cmake -G “MinGW Makefiles” ..
```
5. build the project
```
  cmake --build .
```
6. run the program
```
  ./HMM4GA.exe
```
