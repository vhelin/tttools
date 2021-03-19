
/*
 * deflateTT, compresses data that can be uncompressed using inflateTT.
 * This is almost according to the RFC-1951 (DEFLATE Compressed Data Format
 * Specification version 1.3), only the header is different. I was too lazy
 * to implement them the standard way.
 *
 * Programmed by Ville Helin <vhelin#iki.fi> in 2007.
 *
 * This code is under GNU General Public Licence (GPL), version 2, June 1991.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "defines.h"
#include "main.h"


/* the huffman frequencies */
int freqLiterals[286];
int freqDistances[30];
int freqCombined[119];

/* code lengths */
int codeLengthLiterals[286];
int codeLengthDistances[30];
int codeLengthCombined[119];
int codeLengths[286+30];

/* codes */
int codeLiterals[286];
int codeDistances[30];
int codeCombined[119];

/* the huffman trees */
struct node *treeLiterals = NULL;
struct node *treeDistances = NULL;
struct node *treeCombined = NULL;

/* the priority node queue */
struct node *priorityQueue = NULL;

/* the number of extra bits in the compressed data */
const int extraBitsLengths[] = {
  0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
  1, 1, 2, 2, 2, 2, 3, 3, 3, 3,
  4, 4, 4, 4, 5, 5, 5, 5, 0
};
const int extraBitsDistances[] = {
  0, 0, 0, 0, 1, 1, 2, 2, 3, 3,
  4, 4, 5, 5, 6, 6, 7, 7, 8, 8,
  9, 9, 10, 10, 11, 11, 12, 12, 13, 13
};

/* limits for length encoding */
const int extraBitsLimitsLength[] = {
  11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258
};

/* limits for distance encoding */
const int extraBitsLimitsDistance[] = {
  5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577, 32769
};


static void _write_u32(FILE *f, int data) {

  fprintf(f, "%c%c%c%c", data & 0xFF, (data >> 8) & 0xFF, (data >> 16) & 0xFF, (data >> 24) & 0xFF);
}


static void _write_u8(FILE *f, int data) {

  fprintf(f, "%c", data & 0xFF);
}


static void _write_bits(FILE *f, int *outBitsN, int *outBits, int code, int codeLength) {

  int i;

  for (i = 0; i < codeLength; i++) {
    *outBits = (*outBits << 1) | ((code >> (codeLength - i - 1)) & 1);
    *outBitsN += 1;

    if (*outBitsN == 8) {
      _write_u8(f, *outBits);
      *outBits = 0;
      *outBitsN = 0;
    }
  }
}


/*
  Extra               Extra               Extra
  Code Bits Length(s) Code Bits Lengths   Code Bits Length(s)
  ---- ---- ------     ---- ---- -------   ---- ---- -------
  257   0     3       267   1   15,16     277   4   67-82
  258   0     4       268   1   17,18     278   4   83-98
  259   0     5       269   2   19-22     279   4   99-114
  260   0     6       270   2   23-26     280   4  115-130
  261   0     7       271   2   27-30     281   5  131-162
  262   0     8       272   2   31-34     282   5  163-194
  263   0     9       273   3   35-42     283   5  195-226
  264   0    10       274   3   43-50     284   5  227-257
  265   1  11,12      275   3   51-58     285   0    258
  266   1  13,14      276   3   59-66

  Extra           Extra               Extra
  Code Bits Dist  Code Bits   Dist     Code Bits Distance
  ---- ---- ----  ---- ----  ------    ---- ---- --------
  0   0    1     10   4     33-48    20    9   1025-1536
  1   0    2     11   4     49-64    21    9   1537-2048
  2   0    3     12   5     65-96    22   10   2049-3072
  3   0    4     13   5     97-128   23   10   3073-4096
  4   1   5,6    14   6    129-192   24   11   4097-6144
  5   1   7,8    15   6    193-256   25   11   6145-8192
  6   2   9-12   16   7    257-384   26   12  8193-12288
  7   2  13-16   17   7    385-512   27   12 12289-16384
  8   3  17-24   18   8    513-768   28   13 16385-24576
  9   3  25-32   19   8   769-1024   29   13 24577-32768
*/


void priorityQueuePush(struct node *node) {

  struct node *queue = priorityQueue;

  node->smaller = NULL;
  node->larger = NULL;

  if (queue == NULL) {
    priorityQueue = node;
    return;
  }

  while (1) {
    if (node->weight > queue->weight) {
      if (queue->larger != NULL)
        queue = queue->larger;
      else {
        queue->larger = node;
        node->smaller = queue;
        return;
      }
    }
    else {
      node->smaller = queue->smaller;
      node->larger = queue;
      if (queue->smaller != NULL)
        queue->smaller->larger = node;
      queue->smaller = node;

      if (queue == priorityQueue)
        priorityQueue = node;

      return;
    }
  }
}


struct node *priorityQueuePop(void) {

  struct node *node = priorityQueue;

  if (node == NULL)
    return NULL;

  priorityQueue = node->larger;
  if (priorityQueue != NULL)
    priorityQueue->smaller = NULL;

  return node;
}


void propagateCodes(struct node *node, int code, int codeLength) {

  if (node->leaf == YES) {
    node->code = code;
    node->codeLength = codeLength;
    return;
  }

  if (node->left != NULL)
    propagateCodes(node->left, (code << 1) | 1, codeLength + 1);
  if (node->right != NULL)
    propagateCodes(node->right, (code << 1) | 0, codeLength + 1);
}


void huffmanGrowTree(int *frequencies, int n, struct node **tree, int *codeLengths, int *codes) {

  struct node *nodes, *node1, *node2, *node;
  int i;

  /* init the nodes */
  nodes = malloc(sizeof(struct node) * n);
  if (nodes == NULL) {
    fprintf(stderr, "huffmanGrowTree(): Out of memory error [1].\n");
    return;
  }

  for (i = 0; i < n; i++) {
    nodes[i].left = NULL;
    nodes[i].right = NULL;
    nodes[i].literal = i;
    nodes[i].weight = frequencies[i];
    nodes[i].leaf = YES;
    nodes[i].code = 0;
    nodes[i].codeLength = 0;
  }

  /* push the nodes into a priority queue */
  priorityQueue = NULL;

  for (i = 0; i < n; i++) {
    if (nodes[i].weight > 0)
      priorityQueuePush(&nodes[i]);
  }

  /* DEBUG */
  /*
    {
    struct node *node = priorityQueue;

    while (node != NULL) {
    fprintf(stderr, "%d ", node->weight);
    node = node->larger;
    }
    }
  */

  /* create a tree */
  while (1) {
    /* pop the two rarest nodes */
    node1 = priorityQueuePop();
    node2 = priorityQueuePop();

    if (node2 == NULL) {
      /* the tree has been built */
      *tree = node1;
      break;
    }

    /* create a parent node */
    node = malloc(sizeof(struct node));
    if (node == NULL) {
      fprintf(stderr, "huffmanGrowTree(): Out of memory error [2].\n");
      return;
    }

    node->code = 0;
    node->codeLength = 0;
    node->leaf = NO;
    node->literal = 0xDEADBEEF;
    node->weight = node1->weight + node2->weight;
    node->left = node1;
    node->right = node2;

    /* push the parent back to the queue */
    priorityQueuePush(node);
  }

  /* calculate the codes and code lengths */
  propagateCodes(*tree, 0, 0);

  /* collect the codes and code lengths */
  for (i = 0; i < n; i++) {
    codes[i] = nodes[i].code;
    codeLengths[i] = nodes[i].codeLength;
  }

  /* DEBUG */
  /*
    for (i = 0; i < n; i++)
    fprintf(stderr, "%d/%d ", frequencies[i], codeLengths[i]);
  */
}


/* the number of bits we can have in a code */
#define HUFFMAN_CODE_MAX_BITS 32

/* tmp buffers for huffmanRecreateCodes() */
int tmpCount[HUFFMAN_CODE_MAX_BITS], nextCode[HUFFMAN_CODE_MAX_BITS+1];


void huffmanRecreateCodes(int n, int *lengths, int *codes) {

  int i, code, bits, length;

  /* for more information about this, see RFC-1951 section 3.2.2 */

  /* step 1: count the lengths */
  for (i = 0; i < HUFFMAN_CODE_MAX_BITS; i++)
    tmpCount[i] = 0;

  for (i = 0; i < n; i++)
    tmpCount[lengths[i]]++;

  /*
    for (i = 0; i < HUFFMAN_CODE_MAX_BITS; i++)
    fprintf(stderr, "%d ", tmpCount[i]);
    fprintf(stderr, "\n");
  */

  /* step 2: find the numerical value of the smallest code for each code length */
  code = 0;
  tmpCount[0] = 0;
  for (bits = 1; bits <= HUFFMAN_CODE_MAX_BITS; bits++) {
    code = (code + tmpCount[bits - 1]) << 1;
    nextCode[bits] = code;
  }

  /* step 3: assign numerical values to all codes, using consecutive
     values for all codes of the same length with the base
     values determined at step 2. */
  for (i = 0; i < n; i++) {
    length = lengths[i];
    if (length != 0) {
      /*
        fprintf(stderr, "oldCode = %d newCode = %d\n", codes[i], nextCode[length]);
      */
      codes[i] = nextCode[length];
      nextCode[length]++;
    }
  }
}


int main(int argc, char *argv[]) {

  int fileSize, *lz77, lz77Size, j, i, k, m, n, lz77Best, lz77Length, lz77Matches, lz77DuplicateBytes, codeLengthsN, codeLengthMax, codesN;
  unsigned char *data;
  FILE *f;

  if (argc != 3) {
    fprintf(stderr, "deflateTT v1.1 Written by Ville Helin 2007\n");
    fprintf(stderr, "USAGE: %s <IN RAW> <OUT DEF>\n", argv[0]);
    return 1;
  }

  /********************************************************************************/
  /* INPUT */
  /********************************************************************************/

  /* read the RAW file */
  f = fopen(argv[1], "rb");
  if (f == NULL) {
    fprintf(stderr, "main(): Could not open file \"%s\" for reading.\n", argv[1]);
    return 1;
  }

  /* get the file size */
  fseek(f, 0, SEEK_END);
  fileSize = ftell(f);
  fseek(f, 0, SEEK_SET);

  data = malloc(fileSize);
  if (data == NULL) {
    fprintf(stderr, "main(): Out of memory error [1].\n");
    fclose(f);
    return 1;
  }

  fread(data, 1, fileSize, f);
  fclose(f);

  /********************************************************************************/
  /* ~LZ77 */
  /********************************************************************************/

  /* allocate room for the lz77 compressed data */
  lz77 = malloc(sizeof(int) * fileSize * 3);
  if (lz77 == NULL) {
    fprintf(stderr, "main(): Out of memory error [2].\n");
    return 1;
  }

  /* LZ77 */
  i = 0;
  lz77Size = 0;
  lz77Matches = 0;
  lz77DuplicateBytes = 0;

  while (i < fileSize) {
    /* find the longest match */
    j = i - 0x7FFF + 1;
    if (j < 0)
      j = 0;

    lz77Best = -1;
    lz77Length = -1;

    while (j < i) {
      k = j;
      m = i;
      n = 0;
      while (m < fileSize && data[k] == data[m] && n < 258) {
        k++;
        m++;
        n++;
      }

      if (n >= lz77Length) {
        lz77Best = j;
        lz77Length = n;
      }

      j++;
    }

    if (lz77Length > 2) {
      /* we found a good match -> store */
      lz77[lz77Size++] = lz77Length + 254; /* 3 -> 257 */
      lz77[lz77Size++] = i - lz77Best;

      /* count statistics */
      lz77Matches++;
      lz77DuplicateBytes += lz77Length;

      /* move to pointer over the copied area */
      i += lz77Length - 1;
    }
    else {
      /* just output the data byte */
      lz77[lz77Size++] = data[i];
    }

    i++;
  }

  /* output the end marker */
  lz77[lz77Size++] = 256;

  /* print statistics */
  fprintf(stderr, "main(): LZ77: %d utilized matches | %d duplicate bytes.\n", lz77Matches, lz77DuplicateBytes);

  /********************************************************************************/
  /* PREPROCESS LZ77 -> HUFFMAN */
  /********************************************************************************/

  for (i = 0; i < lz77Size; i++) {
    /*
      fprintf(stderr, "main(): %d\n", lz77[i]);
    */

    /* skip plain data bytes (and the end marker) */
    if (lz77[i] < 257)
      continue;

    /* we have a length value -> transform */
    n = lz77[i] - 254;

    /*
      fprintf(stderr, "  length = %d\n", n);
    */

    if (n > 10) {
      for (j = 0; j < 20; j++) {
        if (n < extraBitsLimitsLength[j + 1]) {
          /* store the extra bits in the upper 16 bits */
          lz77[i] = (265 + j) | ((n - extraBitsLimitsLength[j]) << 16);
          break;
        }
      }
      if (j == 20) {
        if (n == 258)
          lz77[i] = 285;
        else
          fprintf(stderr, "main(): Unsupported length %d in LZ77 preprocess.\n", n);
      }
    }

    /* move to distance */
    i++;

    /* transform the distance */
    n = lz77[i];

    /*
      fprintf(stderr, "  distance = %d\n", n);
    */

    if (n == 1)
      lz77[i] = 0;
    else if (n == 2)
      lz77[i] = 1;
    else if (n == 3)
      lz77[i] = 2;
    else if (n == 4)
      lz77[i] = 3;
    else {
      for (j = 0; j < 26; j++) {
        if (n < extraBitsLimitsDistance[j + 1]) {
          /* store the extra bits in the upper 16 bits */
          lz77[i] = (4 + j) | ((n - extraBitsLimitsDistance[j]) << 16);
          break;
        }
      }
      if (j == 26)
        fprintf(stderr, "main(): Unsupported distance %d in LZ77 preprocess.\n", n);
    }
  }

  /********************************************************************************/
  /* HUFFMAN */
  /********************************************************************************/

  /* zero the frequencies */
  for (i = 0; i < 286; i++)
    freqLiterals[i] = 0;
  for (i = 0; i < 30; i++)
    freqDistances[i] = 0;

  /* calculate the frequencies */
  i = 0;
  while (i < lz77Size) {
    /* plain data bytes, and the end marker */
    n = lz77[i++] & 0xFFFF;

    if (n < 257) {
      freqLiterals[n]++;
    }
    else {
      /* length and distance */
      freqLiterals[n]++;
      n = lz77[i++] & 0xFFFF;
      freqDistances[n]++;
    }
  }

  /* DEBUG */
  /*
    for (i = 0; i < 286; i++)
    fprintf(stderr, "main(): Literal %.3d: %d occurrences.\n", i, freqLiterals[i]);
    for (i = 0; i < 30; i++)
    fprintf(stderr, "main(): Distance %.2d: %d occurrences.\n", i, freqDistances[i]);
  */

  /* create the trees */
  huffmanGrowTree(freqLiterals, 286, &treeLiterals, codeLengthLiterals, codeLiterals);
  huffmanGrowTree(freqDistances, 30, &treeDistances, codeLengthDistances, codeDistances);

  /* rewrite the codes, so that the decoder can create them as well using the same code */
  huffmanRecreateCodes(286, codeLengthLiterals, codeLiterals);
  huffmanRecreateCodes(30, codeLengthDistances, codeDistances);

  /********************************************************************************/
  /* COMPRESS CODE LENGTHS */
  /********************************************************************************/

  /* merge the code lengths */
  j = 0;
  for (i = 0; i < 286; i++)
    codeLengths[j++] = codeLengthLiterals[i];
  for (i = 0; i < 30; i++)
    codeLengths[j++] = codeLengthDistances[i];

  /* find the longest code */
  codeLengthMax = 0;
  for (i = 0; i < 286+30; i++) {
    if (codeLengthMax < codeLengths[i])
      codeLengthMax = codeLengths[i];
  }
  codesN = codeLengthMax + 4;

  /* RLE compress the code lengths */
  i = 0;
  j = 0;
  while (i < 286 + 30) {
    m = codeLengths[i];

    if (m > 115) {
      fprintf(stderr, "main(): Got a code length of %d bits (> 115)!\n", m);
      return 1;
    }

    if (m != 0) {
      /* output nonzeros as they are */
      codeLengths[j++] = m;

      /* ... but does the symbol repeat? */
      k = 0;
      while (i + k < 286+30) {
        if (codeLengths[i + k] != m)
          break;
        if (k == 7)
          break;
        k++;
      }

      i++;
      k--;

      if (k < 3)
        continue;

      /*
        fprintf(stderr, "%d repeats %d\n", m, k);
      */

      codeLengths[j++] = (codeLengthMax + 1) | ((k - 3) << 16);
      i += k;
    }
    else {
      /* RLE compress zeros */
      k = 1;
      i++;
      while (i < 286+30 && k < 138) {
        if (codeLengths[i] != 0)
          break;
        k++;
        i++;
      }

      /* got k zeros */
      if (k < 3) {
        while (k > 0) {
          codeLengths[j++] = 0;
          k--;
        }
      }
      else if (k < 11)
        codeLengths[j++] = (codeLengthMax + 2) | ((k -  3) << 16);
      else
        codeLengths[j++] = (codeLengthMax + 3) | ((k - 11) << 16);
    }
  }

  /* now we have k items in the code length array */
  codeLengthsN = j;

  /* calculate the frequencies */
  for (i = 0; i < codesN; i++)
    freqCombined[i] = 0;

  /* calculate the frequencies */
  i = 0;
  while (i < codeLengthsN) {
    /* plain data bytes, and the end marker */
    freqCombined[codeLengths[i++] & 0xFFFF]++;
  }

  /* create the trees */
  huffmanGrowTree(freqCombined, codesN, &treeCombined, codeLengthCombined, codeCombined);

  /* rewrite the codes, so that the decoder can create them as well using the same code */
  huffmanRecreateCodes(codesN, codeLengthCombined, codeCombined);

  /********************************************************************************/
  /* OUTPUT (DEF) */
  /********************************************************************************/

  /* write the DEF file */
  f = fopen(argv[2], "wb");
  if (f == NULL) {
    fprintf(stderr, "main(): Could not open file \"%s\" for writing.\n", argv[2]);
    return 1;
  }

  /* header */
  _write_u8(f, 'D');
  _write_u8(f, 'E');
  _write_u8(f, 'F');
  _write_u8(f, 'c');

  /* unpacked size */
  _write_u32(f, fileSize);

  j = 0;
  m = 0;

  /* the number of code lengths */
  _write_u8(f, codesN);

  /* find the longest code */
  n = 0;
  for (i = 0; i < codesN; i++) {
    if (n < codeLengthCombined[i])
      n = codeLengthCombined[i];
  }

  /* calculate the number of bits we need in order to write out the code length bits */
  if (n < 4)
    k = 2;
  else if (n < 8)
    k = 3;
  else if (n < 16)
    k = 4;
  else if (n < 32)
    k = 5;
  else if (n < 64)
    k = 6;
  else
    k = 7;

  /* the number of bits per code length */
  _write_bits(f, &j, &m, k, 3);

  /* combined code lengths */
  for (i = 0; i < codesN; i++) {
    if (codeLengthCombined[i] > 115)
      fprintf(stderr, "main(): Combined code length entry %d > 115!\n", codeLengthCombined[i]);
    _write_bits(f, &j, &m, codeLengthCombined[i], k);
    /*
      fprintf(stderr, "main(): codeLengthCombined[%d] = %d\n", i, codeLengthCombined[i]);
    */
  }

  /* compress combined code lengths */
  i = 0;
  while (i < codeLengthsN) {
    n = codeLengths[i] & 0xFFFF;

    _write_bits(f, &j, &m, codeCombined[n], codeLengthCombined[n]);

    /*
      fprintf(stderr, "i = %.3d: got %d\n", i, n);
    */

    /* extra bits? */
    if (n > codeLengthMax) {
      if (n == codeLengthMax + 1)
        _write_bits(f, &j, &m, codeLengths[i] >> 16, 2);
      else if (n == codeLengthMax + 2)
        _write_bits(f, &j, &m, codeLengths[i] >> 16, 3);
      else if (n == codeLengthMax + 3)
        _write_bits(f, &j, &m, codeLengths[i] >> 16, 7);
      else
        fprintf(stderr, "main(): Internal error, combined code %d is not supported!\n", n);

      /*
        fprintf(stderr, "main(): n = %d\n", codeLengths[i] >> 16);
      */
    }

    i++;
  }

  /* compress payload */
  i = 0;
  while (i < lz77Size) {
    n = lz77[i] & 0xFFFF;

    _write_bits(f, &j, &m, codeLiterals[n], codeLengthLiterals[n]);

    if (n > 256) {
      /* length extra bits */
      k = extraBitsLengths[n - 257];
      if (k > 0) {
        n = lz77[i] >> 16;
        _write_bits(f, &j, &m, n, k);
      }
      i++;

      /* distance */
      n = lz77[i] & 0xFFFF;
      _write_bits(f, &j, &m, codeDistances[n], codeLengthDistances[n]);

      /* distance extra bits */
      k = extraBitsDistances[n];
      if (k > 0) {
        n = lz77[i] >> 16;
        _write_bits(f, &j, &m, n, k);
      }
      i++;
    }
    else
      i++;
  }

  /* write out the last, remaining bits */
  if (j != 0)
    _write_u8(f, m << (8 - j));

  i = ftell(f);

  fclose(f);

  fprintf(stderr, "main(): Original size = %dB, deflated size = %dB -> Got rid of %.2f%%.\n", fileSize, i, 100 - (i*100.0f / fileSize));

  return 0;
}
