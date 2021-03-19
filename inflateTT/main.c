
/*
 * inflateTT, decompresses data that has been compressed using deflateTT.
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


/* code lengths */
int codeLengthLiterals[286];
int codeLengthDistances[30];
int codeLengthCombined[119];

/* codes */
int codeLiterals[286];
int codeDistances[30];
int codeCombined[119];

/* the huffman trees */
struct node *treeLiterals = NULL;
struct node *treeDistances = NULL;
struct node *treeCombined = NULL;

/* the free tree nodes */
int treeNodesFreeCurrent = 0;
struct node treeNodesFree[1024];

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

/* base values */
const int baseValueLengths[] = {
  3, 4, 5, 6, 7, 8, 9, 10, 11, 13,
  15, 17, 19, 23, 27, 31, 35, 43, 51, 59,
  67, 83, 99, 115, 131, 163, 195, 227, 258
};
const int baseValueDistances[] = {
  1, 2, 3, 4, 5, 7, 9, 13, 17, 25,
  33, 49, 65, 97, 129, 193, 257, 385, 513, 769,
  1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577
};


/*
  Extra               Extra               Extra
  Code Bits Length(s) Code Bits Lengths   Code Bits Length(s)
  ---- ---- ------     ---- ---- -------   ---- ---- -------
  257   0     3        267   1   15,16     277   4   67-82
  258   0     4        268   1   17,18     278   4   83-98
  259   0     5        269   2   19-22     279   4   99-114
  260   0     6        270   2   23-26     280   4  115-130
  261   0     7        271   2   27-30     281   5  131-162
  262   0     8        272   2   31-34     282   5  163-194
  263   0     9        273   3   35-42     283   5  195-226
  264   0    10        274   3   43-50     284   5  227-257
  265   1  11,12       275   3   51-58     285   0    258
  266   1  13,14       276   3   59-66

  Extra           Extra               Extra
  Code Bits Dist  Code Bits   Dist     Code Bits Distance
  ---- ---- ----  ---- ----  ------    ---- ---- --------
    0   0    1     10   4     33-48     20    9   1025-1536
    1   0    2     11   4     49-64     21    9   1537-2048
    2   0    3     12   5     65-96     22   10   2049-3072
    3   0    4     13   5     97-128    23   10   3073-4096
    4   1   5,6    14   6    129-192    24   11   4097-6144
    5   1   7,8    15   6    193-256    25   11   6145-8192
    6   2   9-12   16   7    257-384    26   12  8193-12288
    7   2  13-16   17   7    385-512    27   12 12289-16384
    8   3  17-24   18   8    513-768    28   13 16385-24576
    9   3  25-32   19   8   769-1024    29   13 24577-32768
*/ 


void huffmanConstructTree(struct node **root, int *codes, int *codeLengths, int n) {

  struct node *node;
  int i, j, code;

  /* init the root node */
  *root = &treeNodesFree[treeNodesFreeCurrent++];
  (*root)->left = NULL;
  (*root)->right = NULL;
  (*root)->literal = -1;

  /* add the leaves to the tree */
  for (i = 0; i < n; i++) {
    /* skip dead leaves */
    if (codeLengths[i] == 0)
      continue;

    /* walk to the leaf's position */
    node = *root;
    code = codes[i];
    for (j = codeLengths[i] - 1; j >= 0; j--) {
      if (((code >> j) & 1) == 0) {
        /* left */
        if (node->left != NULL)
          node = node->left;
        else {
          node->left = &treeNodesFree[treeNodesFreeCurrent++];
          node = node->left;
          node->left = NULL;
          node->right = NULL;
          node->literal = -1;
        }
      }
      else {
        /* right */
        if (node->right != NULL)
          node = node->right;
        else {
          node->right = &treeNodesFree[treeNodesFreeCurrent++];
          node = node->right;
          node->left = NULL;
          node->right = NULL;
          node->literal = -1;
        }
      }
    }

    /* turn this node into a leaf */
    node->literal = i;
  }

  /*
    fprintf(stderr, "HUFFMANCONSTRUCTTREE: Used %d nodes.\n", treeNodesFreeCurrent);
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

  int fileSize, i, j, k, m, n, length, b, e, distance, inflatedSize, codesN, bPrevious;
  unsigned char *data, *tmp;
  struct node *node;
  FILE *f;

  if (argc != 3) {
    fprintf(stderr, "inflateTT v1.1 Written by Ville Helin 2007\n");
    fprintf(stderr, "USAGE: %s <IN DEF> <OUT RAW>\n", argv[0]);
    return 1;
  }

  /********************************************************************************/
  /* INPUT */
  /********************************************************************************/

  /* read the DEF file */
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
  /* HUFFMAN */
  /********************************************************************************/

  /* check header */
  if (data[0] != 'D' || data[1] != 'E' || data[2] != 'F' || data[3] != 'c') {
    fprintf(stderr, "main(): File \"%s\" doesn't start with \"DEFc\".\n", argv[1]);
    return 1;
  }

  i = 4;

  /* parse inflated size */
  inflatedSize = data[i] | (data[i+1] << 8) | (data[i+2] << 16) | (data[i+3] << 24);
  i += 4;

  tmp = malloc(inflatedSize);
  if (tmp == NULL) {
    fprintf(stderr, "main(): Out of memory error [2].\n");
    return 1;
  }

  /*
    fprintf(stderr, "main(): Inflated size = %d\n", inflatedSize);
  */

  /* read the number of code lengths */
  codesN = data[i++];

  /* read bits per code length */
  k = 7;

  m = 0;
  for (n = 0; n < 3; n++) {
    /* parse one bit */
    b = (data[i] >> k) & 1;
    k--;

    m = (m << 1) | b;
  }

  /*
    fprintf(stderr, "main(): Number of items = %d. Bits per item = %d.\n", codesN, m);
  */

  /* read the combined code lengths */
  for (j = 0; j < codesN; j++) {
    codeLengthCombined[j] = 0;

    for (n = 0; n < m; n++) {
      /* parse one bit */
      if (k == -1) {
        k = 7;
        i++;
      }

      b = (data[i] >> k) & 1;
      k--;

      codeLengthCombined[j] = (codeLengthCombined[j] << 1) | b;
    }

    /*
      fprintf(stderr, "main(): codeLengthCombined[%d] = %d\n", j, codeLengthCombined[j]);
    */
  }

  /* create the codes from code lengths */
  huffmanRecreateCodes(codesN, codeLengthCombined, codeCombined);

  /* free all huffman tree nodes */
  treeNodesFreeCurrent = 0;

  /* build the huffman tree */
  huffmanConstructTree(&treeCombined, codeCombined, codeLengthCombined, codesN);

  /* inflate */
  j = 0;
  bPrevious = 0;
  while (j < 286 + 30) {
    node = treeCombined;

    while (node->literal < 0) {
      /* parse one bit */
      if (k == -1) {
        k = 7;
        i++;
      }

      if ((data[i] >> k) & 1)
        node = node->right;
      else
        node = node->left;

      k--;
    }

    b = node->literal;

    /*
      fprintf(stderr, "i = %.3d: got %d\n", j, b);
    */

    if (b <= codesN - 4)
      n = 1;
    else if (b == codesN - 4 + 1)
      n = 2;
    else if (b == codesN - 4 + 2)
      n = 3;
    else
      n = 7;

    /* pump more bits? */
    if (n > 1) {
      if (n == 2 || n == 3)
        m = 3;
      else
        m = 11;

      e = 0;

      while (n > 0) {
        /* parse one bit */
        if (k == -1) {
          k = 7;
          i++;
        }

        e = (e << 1) | ((data[i] >> k) & 1);

        k--;
        n--;
      }

      /*
        fprintf(stderr, "e = %d\n", e);
      */

      n = e + m;

      if (b == codesN - 4 + 1) {
        b = bPrevious;
        /*
          fprintf(stderr, "%d repeats %d\n", b, n);
        */
      }
      else
        b = 0;
    }

    while (n > 0) {
      if (j < 286)
        codeLengthLiterals[j] = b;
      else
        codeLengthDistances[j - 286] = b;
      j++;
      n--;
    }

    bPrevious = b;
  }

  /* create the codes from code lengths */
  huffmanRecreateCodes(286, codeLengthLiterals, codeLiterals);
  huffmanRecreateCodes(30, codeLengthDistances, codeDistances);

  /* free all huffman tree nodes */
  treeNodesFreeCurrent = 0;

  /* build the huffman trees */
  huffmanConstructTree(&treeLiterals, codeLiterals, codeLengthLiterals, 286);
  huffmanConstructTree(&treeDistances, codeDistances, codeLengthDistances, 30);

  /* inflate */
  j = 0;
  while (1) {
    node = treeLiterals;

    while (node->literal < 0) {
      /* parse one bit */
      if (k == -1) {
        k = 7;
        i++;
      }

      if ((data[i] >> k) & 1)
        node = node->right;
      else
        node = node->left;

      k--;
    }

    b = node->literal;

    /*
      fprintf(stderr, "main(): %d\n", b);
    */

    /* a loose literal? */
    if (b < 256) {
      tmp[j++] = b;
      continue;
    }

    /* end of archive? */
    if (b == 256)
      break;

    /* ... so it is a [length, distance] tuple... */
    b -= 257;

    /* get length */
    length = baseValueLengths[b];

    /* parse the extra bits */
    b = extraBitsLengths[b];
    if (b > 0) {
      e = 0;
      while (b > 0) {
        if (k == -1) {
          k = 7;
          i++;
        }

        e = (e << 1) | ((data[i] >> k) & 1);
        k--;
        b--;
      }

      /* add the extra bits */
      length += e;
    }

    /*
      fprintf(stderr, "  length = %d\n", length);
    */

    /* parse distance */
    node = treeDistances;

    while (node->literal < 0) {
      /* parse one bit */
      if (k == -1) {
        k = 7;
        i++;
      }

      if ((data[i] >> k) & 1)
        node = node->right;
      else
        node = node->left;

      k--;
    }

    b = node->literal;

    /* get length */
    distance = baseValueDistances[b];

    /* parse the extra bits */
    b = extraBitsDistances[b];
    if (b > 0) {
      e = 0;
      while (b > 0) {
        if (k == -1) {
          k = 7;
          i++;
        }

        e = (e << 1) | ((data[i] >> k) & 1);
        k--;
        b--;
      }

      /* add the extra bits */
      distance += e;
    }

    /*
      fprintf(stderr, "  distance = %d\n", distance);
    */

    /* de-lz77 */
    n = j - distance;
    for (b = 0; b < length; b++)
      tmp[j++] = tmp[n++];
  }

  fprintf(stderr, "main(): Orginal size = %d, uncompressed size = %d.\n", inflatedSize, j);

  /********************************************************************************/
  /* OUTPUT (RAW) */
  /********************************************************************************/

  /* write the RAW file */
  f = fopen(argv[2], "wb");
  if (f == NULL) {
    fprintf(stderr, "main(): Could not open file \"%s\" for writing.\n", argv[2]);
    return 1;
  }

  fwrite(tmp, 1, j, f);

  fclose(f);

  return 0;
}
