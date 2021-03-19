
/*
 * inflateTT-MP, decompresses data that has been compressed using deflateTT.
 * This is almost according to the RFC-1951 (DEFLATE Compressed Data Format
 * Specification version 1.3), only the header is different. I was too lazy
 * to implement them the standard way.
 *
 * Programmed by Ville Helin <vhelin#iki.fi> in 2007.
 *
 * This code is under GNU Lesser General Public Licence (LGPL), version 2.1,
 * February 1999.
 */

#include <stdlib.h>
#include <stdio.h>

#include "inflate.h"


/* the number of extra bits in the compressed data */
static const int extraBitsLengths[] = {
  0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
  1, 1, 2, 2, 2, 2, 3, 3, 3, 3,
  4, 4, 4, 4, 5, 5, 5, 5, 0
};
static const int extraBitsDistances[] = {
  0, 0, 0, 0, 1, 1, 2, 2, 3, 3,
  4, 4, 5, 5, 6, 6, 7, 7, 8, 8,
  9, 9, 10, 10, 11, 11, 12, 12, 13, 13
};

/* base values */
static const int baseValueLengths[] = {
  3, 4, 5, 6, 7, 8, 9, 10, 11, 13,
  15, 17, 19, 23, 27, 31, 35, 43, 51, 59,
  67, 83, 99, 115, 131, 163, 195, 227, 258
};
static const int baseValueDistances[] = {
  1, 2, 3, 4, 5, 7, 9, 13, 17, 25,
  33, 49, 65, 97, 129, 193, 257, 385, 513, 769,
  1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577
};


static void huffmanConstructTree(struct InflateNode **root, int *codes, int *codeLengths, int n, struct InflateContext *context) {

  struct InflateNode *node;
  int i, j, code;

  /* init the root node */
  *root = &context->treeNodesFree[context->treeNodesFreeCurrent++];
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
          node->left = &context->treeNodesFree[context->treeNodesFreeCurrent++];
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
          node->right = &context->treeNodesFree[context->treeNodesFreeCurrent++];
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
    fprintf(stderr, "HUFFMANCONSTRUCTTREE: Used %d nodes.\n", context->treeNodesFreeCurrent);
  */
}


static void huffmanRecreateCodes(int n, int *lengths, int *codes, struct InflateContext *context) {

  int i, code, bits, length;

  /* for more information about this, see RFC-1951 section 3.2.2 */

  /* step 1: count the lengths */
  for (i = 0; i < HUFFMAN_CODE_MAX_BITS; i++)
    context->tmpCount[i] = 0;

  for (i = 0; i < n; i++)
    context->tmpCount[lengths[i]]++;

  /*
    for (i = 0; i < HUFFMAN_CODE_MAX_BITS; i++)
    fprintf(stderr, "%d ", context->tmpCount[i]);
    fprintf(stderr, "\n");
  */

  /* step 2: find the numerical value of the smallest code for each code length */
  code = 0;
  context->tmpCount[0] = 0;
  for (bits = 1; bits <= HUFFMAN_CODE_MAX_BITS; bits++) {
    code = (code + context->tmpCount[bits - 1]) << 1;
    context->nextCode[bits] = code;
  }

  /* step 3: assign numerical values to all codes, using consecutive
     values for all codes of the same length with the base
     values determined at step 2. */
  for (i = 0; i < n; i++) {
    length = lengths[i];
    if (length != 0) {
      /*
        fprintf(stderr, "oldCode = %d newCode = %d\n", codes[i], context->nextCode[length]);
      */
      codes[i] = context->nextCode[length];
      context->nextCode[length]++;
    }
  }
}


int inflate(unsigned char *data, unsigned char *output, struct InflateContext *context) {

  int i, j, k, m, n, length, b, e, distance, inflatedSize, codesN, bPrevious;
  struct InflateNode *node;

  /* reset the context */
  context->treeLiterals = NULL;
  context->treeDistances = NULL;
  context->treeCombined = NULL;

  /********************************************************************************/
  /* HUFFMAN */
  /********************************************************************************/

  /* check header */
  if (data[0] != 'D' || data[1] != 'E' || data[2] != 'F' || data[3] != 'c')
    return INFLATE_WRONG_HEADER;

  i = 4;

  /* parse inflated size */
  inflatedSize = data[i] | (data[i+1] << 8) | (data[i+2] << 16) | (data[i+3] << 24);
  i += 4;

  /*
    fprintf(stderr, "inflate(): Inflated size = %d\n", inflatedSize);
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
    fprintf(stderr, "inflate(): Number of items = %d. Bits per item = %d.\n", codesN, m);
  */

  /* read the combined code lengths */
  for (j = 0; j < codesN; j++) {
    context->codeLengthCombined[j] = 0;

    for (n = 0; n < m; n++) {
      /* parse one bit */
      if (k == -1) {
        k = 7;
        i++;
      }

      b = (data[i] >> k) & 1;
      k--;

      context->codeLengthCombined[j] = (context->codeLengthCombined[j] << 1) | b;
    }

    /*
      fprintf(stderr, "inflate(): context->codeLengthCombined[%d] = %d\n", j, context->codeLengthCombined[j]);
    */
  }

  /* create the codes from code lengths */
  huffmanRecreateCodes(codesN, context->codeLengthCombined, context->codeCombined, context);

  /* free all huffman tree nodes */
  context->treeNodesFreeCurrent = 0;

  /* build the huffman tree */
  huffmanConstructTree(&context->treeCombined, context->codeCombined, context->codeLengthCombined, codesN, context);

  /* inflate */
  j = 0;
  bPrevious = 0;
  while (j < 286 + 30) {
    node = context->treeCombined;

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
        context->codeLengthLiterals[j] = b;
      else
        context->codeLengthDistances[j - 286] = b;
      j++;
      n--;
    }

    bPrevious = b;
  }

  /* create the codes from code lengths */
  huffmanRecreateCodes(286, context->codeLengthLiterals, context->codeLiterals, context);
  huffmanRecreateCodes(30, context->codeLengthDistances, context->codeDistances, context);

  /* free all huffman tree nodes */
  context->treeNodesFreeCurrent = 0;

  /* build the huffman trees */
  huffmanConstructTree(&context->treeLiterals, context->codeLiterals, context->codeLengthLiterals, 286, context);
  huffmanConstructTree(&context->treeDistances, context->codeDistances, context->codeLengthDistances, 30, context);

  /* inflate */
  j = 0;
  while (1) {
    node = context->treeLiterals;

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
      fprintf(stderr, "inflate(): %d\n", b);
    */

    /* a loose literal? */
    if (b < 256) {
      output[j++] = b;
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
    node = context->treeDistances;

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
      output[j++] = output[n++];
  }

  /*
    fprintf(stderr, "inflate(): Orginal size = %d, uncompressed size = %d.\n", inflatedSize, j);
  */

  return INFLATE_OK;
}
