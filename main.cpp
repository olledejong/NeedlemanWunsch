/*
 * Initial, simplistic version of the Needleman Wunsch algorithm for alignment.
 * No fancy stuff, just a scoring of 1 for a match, -1 for a mismatch and a linearly
 * increasing gap penalty.
 */
#include <iostream>
#include <algorithm>
#include <chrono>
#include <vector>

using namespace std;
using namespace chrono;

// Globals
int matchScore, mismatchScore, gapScore;

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
    // TODO; can this be done more efficiently?
    for (int i = 0; i <= lengthA; i++) nwMatrix[i][0] = -i * gapScore;
    for (int i = 0; i <= lengthB; i++) nwMatrix[0][i] = -i * gapScore;
    // Loop through every cell except for the already filled first row and column. This is done column by column.
    for (int i = 1; i <= lengthA; i++) {
        for (int j = 1; j <= lengthB; j++) {
            // If the previous characters are equal, assign S to matchScore (1), else to mismatchScore (-1)
            bool charsMatch = SeqA[i - 1] == SeqB[j - 1];
            int S = charsMatch ? matchScore : -mismatchScore;
            // Check which movement (down, right or diagonal) results in the highest score for the iterated cell.
            // and set the iterated cell to that value. There are three options:
            // 1. A match/mismatch results in the highest score (previous diagonal + matchScore or - mismatchScore)
            //    Depending on whether it is a match or a mismatch.
            // 2. A vertical gap results in the highest score (the above neighbouring cell - gapScore)
            // 3. A horizontal gap results in the highest score (the left neighbouring cell - gapScore)
            nwMatrix[i][j] = max({
                nwMatrix[i - 1][j - 1] + S,
                nwMatrix[i - 1][j] - gapScore,
                nwMatrix[i][j - 1] - gapScore
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
 * This function loops until there are no characters left to process. During the loops we step through
 * the already generated Needleman-Wunsch matrix called nwMatrix. Every loop we determine what direction
 * has the highest value, we take that route, and logically add the characters to the return strings.
 *
 * @return the optimal alignment based on the scores defined in main()
 */
pair<string, string> getAlignment(const vector<vector<int>> &nwMatrix, const string &SeqA, const string &SeqB) {
    string alignmentPartA, alignmentPartB;
    int iA = SeqA.size(), iB = SeqB.size();
    // loop as long as both strings have not been fully processed
    while (iA != 0 || iB != 0) {
        // if the end of string A has been hit, then that strand can only produce a gap.
        if (iA == 0) {
            alignmentPartA += '-';
            alignmentPartB += SeqB[iB - 1];
            iB--;
            // if the end of string B has been hit, then that strand can only produce a gap.
        } else if (iB == 0) {
            alignmentPartA += SeqA[iA - 1];
            alignmentPartB += '-';
            iA--;
            // compare two characters
        } else {
            // if the next characters of strings are equal, set S = matchScore, else S = mismatchScore
            bool charsMatch = SeqA[iA - 1] == SeqB[iB - 1];
            int S = charsMatch ? matchScore : -mismatchScore;
            // if the previous diagonal cell + S is equal to the iterated cell, then this must be the
            if (nwMatrix[iA][iB] == (nwMatrix[iA - 1][iB - 1] + S )) {
                alignmentPartA += SeqA[iA - 1];
                alignmentPartB += SeqB[iB - 1];
                iA--; iB--;
                // if the next characters of the strings are not equal, then check if the value of the above cell has a
                // higher value than the value of the cell to the left. If so, add a gap to sequence B.
            } else if (nwMatrix[iA - 1][iB] > nwMatrix[iA][iB - 1]) {
                alignmentPartA += SeqA[iA - 1];
                alignmentPartB += '-';
                // only subtract from the iterator you add an actual character to
                iA--;
                // if the value of cell to the left of the iterated cell is higher than the value of
                // the cell above the iterated cell, then add a gap in sequence A.
            } else {
                alignmentPartA += '-';
                alignmentPartB += SeqB[iB - 1];
                iB--;
            }
        }
    }
    // since we handled the sequences backwards, we have to reverse them for correct representation
    reverse(alignmentPartA.begin(), alignmentPartA.end());
    reverse(alignmentPartB.begin(), alignmentPartB.end());

    // return the optimal alignment pair
    return {alignmentPartA, alignmentPartB};
}

// program starts here
int main() {
    matchScore = 1;
    mismatchScore = 1;
    gapScore = 1;
    string SeqA = "CACATA", SeqB = "CAGCTATA";
    int lengthA = SeqA.size(), lengthB = SeqB.size(); // set lengths of seqs globally

    // get time at beginning of alignment
    auto t1 = high_resolution_clock::now();

    // compute the needleman-wunsch matrix for the given sequences
    auto nwMatrix = needlemanWunsch(SeqA, SeqB);
    // uncomment the next line to see the generated matrix as well
//    printMatrix(nwMatrix);

    // get the optimal alignment
    pair<string, string> alignment = getAlignment(nwMatrix, SeqA, SeqB);
    printf("With rules; match: +%d, mismatch: -%d and gap: -%d, the given optimal alignment with a Needleman-Wunsch"
           " score of %d is:\n\n\t%s\n\t%s\n\n", matchScore, mismatchScore, gapScore, nwMatrix[lengthA][lengthB],
           alignment.first.c_str(), alignment.second.c_str());

    // get time at end of alignment
    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;
    printf("This alignment took %f ms.", ms_double.count());

    exit(EXIT_SUCCESS);
}