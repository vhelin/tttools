#ifndef HUFFMAN_H
#define HUFFMAN_H

/* Build canonical Huffman codes from symbol frequencies. lengths[i] is 0 when frequencies[i] is 0.
   Every used alphabet gets at least two non-zero lengths, each at most 15, so a single used symbol
   is stored as a one-bit code. Returns 0 on success and -1 if the arguments are unusable, memory
   runs out, or a valid 15-bit code cannot be built. */
int huffmanBuildCodes(const int *frequencies, int n, int *codeLengths, int *codes);

#endif
