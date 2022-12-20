## Sequence Alignment
### A C++ implementation of both the Needleman-Wunsch and Smith-Waterman algorithms

The Needleman-Wunsch and Smith-Waterman algorithms are two algorithms that are very famous and are still widely used when performing sequence alignment.
The former is used for performing a global aligment, while the latter is used to perform local alignements. Below is a table that holds their specific differences.  


|                | Needleman-Wunsch algorithm                                                 | Smith-Waterman algorithm                                |
|----------------|----------------------------------------------------------------------------|---------------------------------------------------------|
| Initialization | First row and first column are subject to gap penalty                      | First row and first column are set to 0                 |
| Scoring        | Score can be negative                                                      | Negative score is set to 0                              |
| Traceback      | Begin with the cell at the lower right of the matrix, end at top left cell | Begin with the highest score, end when 0 is encountered |

