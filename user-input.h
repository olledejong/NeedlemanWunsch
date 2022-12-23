#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <filesystem>
#include <algorithm>
#include <chrono>
#include <vector>
using namespace std;
using namespace chrono;

// ====== Globals ====== //
string filename;          // relative path to the FASTA file
bool doPrintMatrix;      // user setting; print computed matrix or not
int gapPenalty;          // user setting; height of gap penalty
char globalOrLocal;      // user setting; determines whether doing global or local alignment
bool proteinSeqs;        // true; working with protein seqs, else; working with dna
// ===================== //

/**
 * Collects user settings through the command line and stores these in globals.
 */
void getUserSettings() {
    printf("NOTE: If the sequences are substantial, the alignment might take a while.\n");
    printf("Please enter the ALBOSULTE path to the FASTA file containing both sequences:\n");
    getline(cin, filename);
    printf("Would you like perform a global, or local alignment? [ g/l ]: "); cin >> globalOrLocal;
    printf("Do you want the matrix to be printed? [ 0/1 ]: "); cin >> doPrintMatrix;
    printf("What should be the (linear) gap penalty?: "); cin >> gapPenalty;
    if (gapPenalty < 1) gapPenalty = 1;
}

/**
 * Reads the first two sequences from a by the user given FASTA file.
 * If the file does not exist, or if the file does not contain at least two sequences, then
 * an error is thrown and the program is terminated.
 *
 * @return a pair containing sequence A and B
 */
pair<string, string> readFastaFile() {
    struct stat buffer{};
    if (stat(filename.c_str(), &buffer) != 0) throw logic_error("That file does not exist, aborting..");

    // check if file is in correct format and contains two sequences
    ifstream fastaFile(filename);
    string line, SeqA, SeqB;
    int count = 0;
    while ( getline(fastaFile, line) ) {
        if (line[0] != '>') {
            switch(count) {
                case 1:
                    SeqA += line;
                    break;
                case 2:
                    SeqB += line;
                    break;
                default:
                    continue;
            }
        } else {
            ++count;
        }
    }
    // close the file
    fastaFile.close();

    // if not both strings are filled, throw an error
    if (SeqA.empty() || SeqB.empty()) throw logic_error("The file does not contain two sequences, aborting..");
    proteinSeqs = (SeqA.find('M') != std::string::npos && SeqB.find('M') != std::string::npos);

    return {SeqA, SeqB};
}