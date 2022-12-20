/*
 * Implementation of the Needleman and Wunsch algorithm for alignment of two strings.
 * The user can either choose to perform a global, or local alignment. The user can also set the gap penalty.
 *
 * Contents
 *    [ charchar-maps.h ][ contains two scoring matrices (in map form), one for DNA, and one for protein ]
 *    [ user-input.h ][ contains logic related to user-input (filepath, gap-penalty, FASTA file reading) ]
 *    [ main.cpp ][ contains the main logic responsible for generating the (global/local) alignment ]
 *
 * By Olle de Jong on 20/12/2022
 */

// ====== Includes ======= //
#include <iostream>        //
#include <algorithm>       //
#include <chrono>          //
#include <vector>          //
#include "charchar-maps.h" // scoring matrixes for global alignment
#include "user-input.h"    // contains functions related to user-input
using namespace std;       //
using namespace chrono;    //
// ======================= //

// ====== Globals ====== //
vector<pair<string,string>> possibleAlignments; // vector that holds pairs of all possible optimal alignments

/**
 * Simple function that generates visual representation of the Needleman-Wunsch matrix and
 * prints it to the console.
 *
 * @param nwMatrix the Needleman-Wunsch matrix
 */
void printMatrix(const vector<vector<int>> &nwMatrix) {
    printf("\n");
    for (int i = 0; i < nwMatrix.size(); ++i) {
        for (int j = 0; j < nwMatrix[0].size(); j++) {
            printf("%d\t", nwMatrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

/**
 * Computes the Needleman-Wunsch matrix for the purpose of global alignment.
 * Scores are taken from the bbScores or aaScores map, which are located in the
 * 'charchar-maps.h' header file.
 *
 * @param SeqA the first string of the alignment
 * @param SeqB the second string of the alignment
 * @return the maximum score of the matrix
 */
vector<vector<int>> needlemanWunsch(const string &SeqA, const string &SeqB) {
    // create the 'empty' matrix
    const int lengthA = SeqA.size(), lengthB = SeqB.size();
    vector<int> row((lengthB + 1), 0);
    vector<vector<int>> nwMatrix((lengthA + 1), row);

    // set values of first and second column to multiple of gap penalty
    for (int i = 0; i <= lengthA; i++) nwMatrix[i][0] = -i * gapPenalty;
    for (int i = 0; i <= lengthB; i++) nwMatrix[0][i] = -i * gapPenalty;

    // loop through all other cells and fill those
    for (int i = 1; i <= lengthA; i++) {
        for (int j = 1; j <= lengthB; j++) {
            // grab the score for the pair of chars
            pair<char, char> charPair = {SeqA[i - 1], SeqB[j - 1]};
            int S = proteinSeqs ? aaScores[charPair] : bbScores[charPair];
            // check which movement (down, right or diagonal) results in the
            // highest score and set iterated cell to it.
            nwMatrix[i][j] = max({
                nwMatrix[i - 1][j - 1] + S,
                nwMatrix[i - 1][j] - gapPenalty,
                nwMatrix[i][j - 1] - gapPenalty
            });
        }
    }
    return nwMatrix;
}

/**
 * This function loops until both sequences have been fully processed. During the loops we step through
 * the already generated Needleman-Wunsch matrix called nwMatrix. Every loop we determine what direction
 * has the highest value, we take that route, and logically add the characters to the return strings.
 *
 * @param nwMatrix the computed Needleman-Wunsch matrix
 * @param SeqA the original sequence A
 * @param SeqB the original sequence B
 * @return the optimal alignment based on the scores defined in main()
 */
void getAlignmentNW(const vector<vector<int>> &nwMatrix,
                    const string &SeqA,
                    const string &SeqB,
                    int iA,
                    int iB,
                    string outA = "",
                    string outB = "",
                    bool altOptAlignment = false) {
    // if flag is true, instead of adding gap in seq A, add a gap in seq B and continue with the alignment
    if (altOptAlignment) {
        outA += SeqA[iA - 1];
        outB += '-';
        iA--;
    }
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
            pair<char, char> charPair = {SeqA[iA - 1], SeqB[iB - 1]};
            int S = proteinSeqs ? aaScores[charPair] : bbScores[charPair];
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
                iA--;
            // if the value of cell to the left of the iterated cell is higher than the value of
            // the cell above the iterated cell, then add a gap in sequence A.
            } else if (nwMatrix[iA - 1][iB] < nwMatrix[iA][iB - 1]) {
                outA += '-';
                outB += SeqB[iB - 1];
                iB--;
            } else if (nwMatrix[iA - 1][iB] == nwMatrix[iA][iB - 1]) {
                // every time there are two equal options in gap creation, then get the alternate alignment as well
                getAlignmentNW(nwMatrix, SeqA, SeqB, iA, iB, outA, outB, true);
                // prefer gap creation in sequence A.
                outA += '-';
                outB += SeqB[iB - 1];
                iB--;
            }
        }
    }
    // since we handled the sequences backwards, we have to reverse them for correct representation
    reverse(outA.begin(), outA.end());
    reverse(outB.begin(), outB.end());

    // add the optimal alignment pair
    possibleAlignments.emplace_back(outA, outB);
}

/**
 * Computes the Smith-Waterman matrix for the purpose of local alignment.
 * In contrary to Needleman Wunsch, this algorithm usualy does not work with scoring
 * matrices, instead, the scores are +1 and -1 for a match and mismatch, respectively.
 * Also, all negative outcome values are set to 0. So are the first row and column.
 * NOTE: The gap penalty can still be set by the user.
 *
 * @param SeqA
 * @param SeqB
 * @return
 */
vector<vector<int>> smithWaterman(const string &SeqA, const string &SeqB) {
    // create the 'empty' matrix
    const int lengthA = SeqA.size(), lengthB = SeqB.size();
    vector<int> row((lengthB + 1), 0);
    vector<vector<int>> swMatrix((lengthA + 1), row);

    // set values of first and second column to multiple of gap penalty
    for (int i = 0; i <= lengthA; i++) swMatrix[i][0] = 0;
    for (int i = 0; i <= lengthB; i++) swMatrix[0][i] = 0;

    // loop through all other cells and fill those
    for (int i = 1; i <= lengthA; i++) {
        for (int j = 1; j <= lengthB; j++) {
            // grab the score for the pair of chars
            int S = (SeqA[i - 1] == SeqB[j - 1]) ? 1 : -1;
            // check which movement (down, right or diagonal) results in the
            // highest score and set iterated cell to it.
            int maxVal = max({
                swMatrix[i - 1][j - 1] + S,
                swMatrix[i - 1][j] - gapPenalty,
                swMatrix[i][j - 1] - gapPenalty
            });
            swMatrix[i][j] = maxVal < 0 ? 0 : maxVal;
        }
    }
    return swMatrix;
}

/**
 * This function finds the highest value in the matrix. After the highest value is found,
 * the matrix is scanned for all positions that hold this value. These positions are stored
 * within a vector and returned.
 *
 * @param swMatrix the Smith-Waterman matrix
 * @return a vector of pairs that contains the locations of the max values in the matrix
 */
vector<pair<int, int>> findMaxValPos(const vector<vector<int>> &swMatrix) {
    // find the highest score(s) in the matrix
    vector<pair<int,int>> maxValPoses;
    pair<int,int> maxValPos;
    int maxVal = INT_MIN;

    // find the highest value in the matrix
    for (int i = 0; i < swMatrix.size(); ++i) {
        for (int j = 0; j < swMatrix[0].size(); ++j) {
            if (swMatrix[i][j] > maxVal) {
                maxVal = swMatrix[i][j];
            }
        }
    }
    // store all positions that have that highest value
    for (int i = 0; i < swMatrix.size(); ++i) {
        for (int j = 0; j < swMatrix[0].size(); ++j) {
            if (swMatrix[i][j] == maxVal) {
                maxValPoses.emplace_back(i, j);
            }
        }
    }
    return maxValPoses;
}

/**
 * This function generates all possible local alignments from the swMatrix. It does this from
 * every single one of the maximum value positions in the matrix.
 *
 * @param swMatrix the smith-waterman matrix for local alignments
 * @param SeqA the first to be aligned sequence
 * @param SeqB the second to be aligned sequence
 */
void getAlignmentSW(const vector<vector<int>> &swMatrix,
                    const string &SeqA,
                    const string &SeqB) {

    // get the all positions in the matrix that hold the maximum value
    vector<pair<int,int>> maxValPoses = findMaxValPos(swMatrix);

    // while not all possible alignments have been handled
    while (!maxValPoses.empty()) {
        string outA;
        string outB;
        // the startposition of the local alignment
        pair<int, int> startCell = maxValPoses[maxValPoses.size() - 1];
        // sequence iterators
        int iA = startCell.first - 1;
        int iB = startCell.second - 1;
        // add the first chars to the aligment
        outA += SeqA[iA];
        outB += SeqB[iB];
        // handle the rest of the sequences
        while (swMatrix[iA][iB] != 0) {
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
                int S = (SeqA[iA - 1] == SeqB[iB - 1]) ? 1 : -1;
                // if the previous diagonal cell + S is equal to the iterated cell, then this must be the optimal way
                if (swMatrix[iA][iB] == (swMatrix[iA - 1][iB - 1] + S)) {
                    outA += SeqA[iA - 1];
                    outB += SeqB[iB - 1];
                    iA--;
                    iB--;
                // if the next characters of the strings are not equal, then check if the value of the above cell has a
                // higher value than the value of the cell to the left. If so, add a gap to sequence B.
                } else if (swMatrix[iA - 1][iB] > swMatrix[iA][iB - 1]) {
                    outA += SeqA[iA - 1];
                    outB += '-';
                    iA--;
                // if the value of cell to the left of the iterated cell is higher than the value of
                // the cell above the iterated cell, then add a gap in sequence A.
                } else {
                    outA += '-';
                    outB += SeqB[iB - 1];
                    iB--;
                }
            }
        }
        // since we handled the sequences backwards, we have to reverse them for correct representation
        reverse(outA.begin(), outA.end());
        reverse(outB.begin(), outB.end());
        // add the alignment pair to the vector
        possibleAlignments.emplace_back(outA, outB);
        // remove the handled starting location of the alignment from the vector
        maxValPoses.pop_back();
    }
}

//=== program starts here ===/
int main() {
    try {
        getUserSettings();
        pair<string, string> inputSeqs = readFastaFile();
        string SeqA = inputSeqs.first, SeqB = inputSeqs.second;
        int lengthA = SeqA.size(), lengthB = SeqB.size();

        //==== time the duration of the actual logic; start timer ====//
        auto t1 = high_resolution_clock::now();

        // report whether we are working with nucleotide or protein sequences
        proteinSeqs ? printf("We're working with amino-acids.\n") : printf("We're working with nucleotides.\n");

        // compute the matrix for the given sequences and get the possible alignments
        if (globalOrLocal == 'g') {
            auto nwMatrix = needlemanWunsch(SeqA, SeqB);
            if (doPrintMatrix) printMatrix(nwMatrix);
            getAlignmentNW(nwMatrix, SeqA, SeqB, lengthA, lengthB);

            printf("The given optimal alignment(s) with a Needleman-Wunsch score of %d is/are:\n",
                   nwMatrix[lengthA][lengthB]);
            for (const auto &item: possibleAlignments) {
                printf("\n\t%s\n\t%s\n", item.first.c_str(), item.second.c_str());
            }
        } else if (globalOrLocal == 'l') {
            auto swMatrix = smithWaterman(SeqA, SeqB);
            if (doPrintMatrix) printMatrix(swMatrix);
            getAlignmentSW(swMatrix, SeqA, SeqB);

            printf("The given local alignment(s) is/are:\n");
            for (const auto &item: possibleAlignments) {
                printf("\n\t%s\n\t%s\n", item.first.c_str(), item.second.c_str());
            }
        } else {
            throw logic_error("Incorrect alignment type input. Aborting..");
        }

        //==== stop timer and report duration ====//
        auto t2 = high_resolution_clock::now();
        duration<double, std::milli> ms_double = t2 - t1;
        printf("\nThis alignment took %f ms.", ms_double.count());

    } catch (const logic_error &error) {
        cerr << error.what() << endl;
        exit(EXIT_FAILURE);
    }
    exit(EXIT_SUCCESS);
}