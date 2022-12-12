/*
 * Initial, simplistic version of the Needleman Wunsch algorithm for alignment.
 * No fancy stuff, just a scoring of 1 for a match, -1 for a mismatch and a linearly
 * increasing gap penalty.
 */
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <filesystem>
#include <algorithm>
#include <chrono>
#include <vector>
#include <map>

using namespace std;
using namespace chrono;

// Globals
int gapPenalty;

// Map for base-wise match scores. Purine -> Purine or Pyrimidine -> Pyrimidine swaps produce a smaller penalty
map<pair<char,char>, int> baseWiseScore {
        {{'A','A'}, 10},
        {{'A','G'}, -1},
        {{'A','C'}, -3},
        {{'A','T'}, -4},
        {{'G','A'}, -1},
        {{'G','G'}, 7},
        {{'G','C'}, -5},
        {{'G','T'}, -3},
        {{'C','A'}, -3},
        {{'C','G'}, -5},
        {{'C','C'}, 9},
        {{'C','T'}, 0},
        {{'T','A'}, -4},
        {{'T','G'}, -3},
        {{'T','C'}, 0},
        {{'T','T'}, 8},
};

/**
 * Compute the Needleman-Wunsch matrix
 * For example: CACATA and CAGCTA, with scores 1, -1 and -1, for match, mismatch and gap respectively:
 *         C   A   C   A   T   A
 *     0  -1  -2  -3  -4  -5  -6
 * C  -1   1   0  -1  -2  -3  -4
 * A  -2   0   2   1   0  -1  -2
 * G  -3  -1   1   1   2   1   0
 * C  -4  -2   0   0   1   1   2
 * T  -5  -3  -1  -1   0   2   1
 * A  -6  -4  -2  -2  -1   1   3
 *
 * @return the maximum score of the matrix
 */
vector<vector<int>> needlemanWunsch(const string &SeqA, const string &SeqB) {
    // create the 'empty' matrix
    const int lengthA = SeqA.size(), lengthB = SeqB.size();
    vector<int> row((lengthB + 1), 0);
    vector<vector<int>> nwMatrix((lengthA + 1), row);
    // Set the values of the first row and the first column to a multiple of the gap penalty
    for (int i = 0; i <= lengthA; i++) nwMatrix[i][0] = -i * gapPenalty;
    for (int i = 0; i <= lengthB; i++) nwMatrix[0][i] = -i * gapPenalty;
    // Loop through every cell except for the already filled first row and column. This is done column by column.
    for (int i = 1; i <= lengthA; i++) {
        for (int j = 1; j <= lengthB; j++) {
            // grab the score for the two compared chars from the map
            int S = baseWiseScore[{SeqA[i - 1], SeqB[j - 1]}];
            // Check which movement (down, right or diagonal) results in the highest score for the iterated cell.
            // and set the iterated cell to that value. There are three options:
            // 1. A match/mismatch results in the highest score (previous diagonal + matchScore or - mismatchScore)
            //    Depending on whether it is a match or a mismatch.
            // 2. A vertical gap results in the highest score (the above neighbouring cell - gapPenalty)
            // 3. A horizontal gap results in the highest score (the left neighbouring cell - gapPenalty)
            nwMatrix[i][j] = max({
                nwMatrix[i - 1][j - 1] + S,
                nwMatrix[i - 1][j] - gapPenalty,
                nwMatrix[i][j - 1] - gapPenalty
            });
        }
    }
    // only return the maximum score of the alignmnent
    return nwMatrix;
}

/**
 * Simple function that generates visual representation of the Needleman-Wunsch matrix and
 * prints it to the console.
 */
void printMatrix(const vector<vector<int>> &nwMatrix) {
    for (int i = 0; i < nwMatrix.size(); ++i) {
        for (int j = 0; j < nwMatrix[0].size(); j++) {
            printf("%d\t", nwMatrix[i][j]);
        }
        printf("\n");
    }
}


/**
 * Once it has been noticed that there are two optimal alignments to be retrieved from the nwMatrix,
 * this function gets is executed. This function builds the alternate optimal alignment (which has the same
 * needleman-wunsch score) by taking the already aligned part, but where the 'original' optimal alignment adds a
 * gap to Sequence A (first sequence), this function adds the gap to Sequence B (second sequence).
 *
 * @param nwMatrix the computed Needleman-Wunsch matrix
 * @param iA row index of needleman-wunsch matrix where two optimal paths split
 * @param iB column index of needleman-wunsch matrix where two optimal paths split
 * @param SeqA the original sequence A
 * @param SeqB the original sequence B
 * @param alignPartA the first part of the already aligned strings (alt alignment is picked up from there)
 * @param alignPartB the second part of the already aligned strings (alt alignment is picked up from there)
 * @return a pair containing the alternate optimal alignment
 */
pair<string,string> getAltAlignment(const vector<vector<int>> &nwMatrix,
                                    int iA,
                                    int iB,
                                    const string &SeqA,
                                    const string &SeqB,
                                    string altOutA,
                                    string altOutB) {

    // instead of adding gap in seq A, add a gap in seq B
    altOutA += SeqA[iA - 1];
    altOutB += '-';
    iA--;

    // finish the alignment
    while (iA != 0 || iB != 0) {
        if (iA == 0) {
            altOutA += '-';
            altOutB += SeqB[iB - 1];
            iB--;
            // if the end of string B has been hit, then that strand can only produce a gap.
        } else if (iB == 0) {
            altOutA += SeqA[iA - 1];
            altOutB += '-';
            iA--;
            // compare two characters
        } else {
            // if the next characters of strings are equal, set S = matchScore, else S = mismatchScore
            // grab the score for the two compared chars from the map
            int S = baseWiseScore[{SeqA[iA - 1], SeqB[iB - 1]}];
            // if the previous diagonal cell + S is equal to the iterated cell, then this must be the
            if (nwMatrix[iA][iB] == (nwMatrix[iA - 1][iB - 1] + S )) {
                altOutA += SeqA[iA - 1];
                altOutB += SeqB[iB - 1];
                iA--; iB--;
                // if the next characters of the strings are not equal, then check if the value of the above cell has a
                // higher value than the value of the cell to the left. If so, add a gap to sequence B.
            } else if (nwMatrix[iA - 1][iB] > nwMatrix[iA][iB - 1]) {
                altOutA += SeqA[iA - 1];
                altOutB += '-';
                // only subtract from the iterator you add an actual character to
                iA--;
                // if the value of cell to the left of the iterated cell is higher than the value of
                // the cell above the iterated cell, then add a gap in sequence A.
            } else {
                altOutA += '-';
                altOutB += SeqB[iB - 1];
                iB--;
            }
        }
    }
    reverse(altOutA.begin(), altOutA.end());
    reverse(altOutB.begin(), altOutB.end());
    return {altOutA, altOutB};
}

/**
 * This function loops until there are no characters left to process. During the loops we step through
 * the already generated Needleman-Wunsch matrix called nwMatrix. Every loop we determine what direction
 * has the highest value, we take that route, and logically add the characters to the return strings.
 *
 * @param nwMatrix the computed Needleman-Wunsch matrix
 * @param SeqA the original sequence A
 * @param SeqB the original sequence B
 * @return the optimal alignment based on the scores defined in main()
 */
vector<pair<string, string>> getAlignment(const vector<vector<int>> &nwMatrix, const string &SeqA, const string &SeqB) {
    vector<pair<string,string>> possibleAlignments;
    string outA, outB;
    int iA = SeqA.size(), iB = SeqB.size();
    // loop as long as both strings have not been fully processed
    while (iA != 0 || iB != 0) {
        // if the end of string A has been hit, then that strand can only produce a gap.
        if (iA == 0) {
            outA += '-';
            outB += SeqB[iB - 1];
            iB--;
            // if the end of string B has been hit, then that strand can only produce a gap.
        } else if (iB == 0) {
            outA += SeqA[iA - 1];
            outB += '-';
            iA--;
            // compare two characters
        } else {
            // grab the score for the two compared chars from the map
            int S = baseWiseScore[{SeqA[iA - 1], SeqB[iB - 1]}];
            // if the previous diagonal cell + S is equal to the iterated cell, then this must be the optimal way
            if (nwMatrix[iA][iB] == (nwMatrix[iA - 1][iB - 1] + S )) {
                outA += SeqA[iA - 1];
                outB += SeqB[iB - 1];
                iA--; iB--;
                // if the next characters of the strings are not equal, then check if the value of the above cell has a
                // higher value than the value of the cell to the left. If so, add a gap to sequence B.
            } else if (nwMatrix[iA - 1][iB] > nwMatrix[iA][iB - 1]) {
                outA += SeqA[iA - 1];
                outB += '-';
                // only subtract from the iterator you add an actual character to
                iA--;
                // if the value of cell to the left of the iterated cell is higher than the value of
                // the cell above the iterated cell, then add a gap in sequence A.
            } else if (nwMatrix[iA - 1][iB] < nwMatrix[iA][iB - 1]) {
                outA += '-';
                outB += SeqB[iB - 1];
                iB--;
            } else if (nwMatrix[iA - 1][iB] == nwMatrix[iA][iB - 1]) {
                // every time there are two equal options in gap creation, then get the alternate alignment as well
                possibleAlignments.push_back(getAltAlignment(nwMatrix, iA, iB, SeqA, SeqB, outA, outB));
                // Prefer gap creation in sequence A.
                outA += '-';
                outB += SeqB[iB - 1];
                iB--;
            }
        }
    }
    // since we handled the sequences backwards, we have to reverse them for correct representation
    reverse(outA.begin(), outA.end());
    reverse(outB.begin(), outB.end());
    possibleAlignments.emplace_back(outA, outB);

    // return the optimal alignment pair
    return possibleAlignments;
}

/**
 * Reads the first two sequences from a by the user given FASTA file.
 * If the file does not exist, or if the file does not contain at least two sequences, then
 * an error is thrown and the program is terminated.
 *
 * @return a pair containing sequence A and B
 */
pair<string, string> readFastaFile() {
    // check if the input file exists at all
    printf("Your current working directory is: %s\n", std::filesystem::current_path().c_str());
    printf("NOTE: If the sequences are substantial, the alignment might take a while.\n"
           "Please enter the relative path to the FASTA file: ");
    string filename; cin >> filename;
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

    return {SeqA, SeqB};
}

// program starts here
int main() {
    try {
        pair<string, string> inputSeqs = readFastaFile();
        string SeqA = inputSeqs.first, SeqB = inputSeqs.second;
        int lengthA = SeqA.size(), lengthB = SeqB.size();

        // set gap score
        gapPenalty = 5;

        // get time at beginning of alignment
        auto t1 = high_resolution_clock::now();

        // compute the needleman-wunsch matrix for the given sequences
        auto nwMatrix = needlemanWunsch(SeqA, SeqB);

        // get the optimal alignment
        vector<pair<string, string>> alignments = getAlignment(nwMatrix, SeqA, SeqB);
        printf("The given optimal alignment(s) with a Needleman-Wunsch score of %d is/are:\n", nwMatrix[lengthA][lengthB]);
        for (const auto &item: alignments) {
            printf("\n\t%s\n\t%s\n", item.first.c_str(), item.second.c_str());
        }

        // report duration of alignment
        auto t2 = high_resolution_clock::now();
        duration<double, std::milli> ms_double = t2 - t1;
        printf("\nThis alignment took %f ms.", ms_double.count());

    } catch (const logic_error &error) {
        cerr << error.what() << endl;
        exit(EXIT_FAILURE);
    }
    exit(EXIT_SUCCESS);
}