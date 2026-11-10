/*
 * inflateTT-MP, decompresses data that has been compressed using deflateTT.
 * Literals, lengths and distances follow RFC 1951. The container is DEFc:
 * MSB-first bits, one block, custom code-length header.
 *
 * Programmed by Ville Helin <vhelin#iki.fi> in 2007.
 * Bounds checks and the 16-bit writer entry point were tightened in 2026.
 *
 * This code is under GNU Lesser General Public Licence (LGPL), version 2.1,
 * February 1999.
 */

#include <limits.h>

#include "inflate.h"


static const int extraBitsLengths[29] = {
  0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
  1, 1, 2, 2, 2, 2, 3, 3, 3, 3,
  4, 4, 4, 4, 5, 5, 5, 5, 0
};
static const int extraBitsDistances[30] = {
  0, 0, 0, 0, 1, 1, 2, 2, 3, 3,
  4, 4, 5, 5, 6, 6, 7, 7, 8, 8,
  9, 9, 10, 10, 11, 11, 12, 12, 13, 13
};
static const int baseValueLengths[29] = {
  3, 4, 5, 6, 7, 8, 9, 10, 11, 13,
  15, 17, 19, 23, 27, 31, 35, 43, 51, 59,
  67, 83, 99, 115, 131, 163, 195, 227, 258
};
static const int baseValueDistances[30] = {
  1, 2, 3, 4, 5, 7, 9, 13, 17, 25,
  33, 49, 65, 97, 129, 193, 257, 385, 513, 769,
  1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577
};


struct bitreader {
  const unsigned char *data;
  int size;
  int index;
  int bit;
};


struct sink {
  unsigned char *bytes;
  volatile unsigned short *vram;
  int n;
  int cap;
};


static int bit_get(struct bitreader *reader, int *out) {

  if (reader->bit < 0) {
    if (reader->size < 1 || reader->index >= reader->size - 1)
      return -1;
    reader->index++;
    reader->bit = 7;
  }
  if (reader->index >= reader->size)
    return -1;
  *out = (reader->data[reader->index] >> reader->bit) & 1;
  reader->bit--;
  return 0;
}


static int bit_getn(struct bitreader *reader, int count, int *out) {

  int value, i, bit;

  value = 0;
  for (i = 0; i < count; i++) {
    if (bit_get(reader, &bit) != 0)
      return -1;
    value = (value << 1) | bit;
  }
  *out = value;
  return 0;
}


static struct InflateNode *alloc_node(struct InflateContext *context) {

  struct InflateNode *node;

  if (context->treeNodesFreeCurrent < 0 ||
      context->treeNodesFreeCurrent >= INFLATE_MAX_TREE_NODES)
    return NULL;

  node = &context->treeNodesFree[context->treeNodesFreeCurrent];
  context->treeNodesFreeCurrent++;
  node->left = NULL;
  node->right = NULL;
  node->literal = -1;
  return node;
}


static int assign_canonical(struct InflateContext *context, int n,
                            const int *lengths, int *codes) {

  int i, bits, length, maxBits;
  unsigned int left, ucode, sum;

  if (n < 0)
    return -1;

  for (i = 0; i <= INFLATE_MAX_BITS; i++)
    context->tmpCount[i] = 0;

  maxBits = 0;
  for (i = 0; i < n; i++) {
    length = lengths[i];
    if (length < 0 || length > INFLATE_MAX_BITS)
      return -1;
    context->tmpCount[length]++;
    if (length > maxBits)
      maxBits = length;
  }

  left = 1;
  for (bits = 1; bits <= maxBits; bits++) {
    left <<= 1;
    if (left < (unsigned int)context->tmpCount[bits])
      return -1;
    left -= (unsigned int)context->tmpCount[bits];
  }

  ucode = 0;
  context->tmpCount[0] = 0;
  for (bits = 1; bits <= maxBits; bits++) {
    sum = ucode + (unsigned int)context->tmpCount[bits - 1];
    if (sum > 0x7FFFFFFFU)
      return -1;
    ucode = sum << 1;
    if (ucode > 0x7FFFFFFFU)
      return -1;
    context->nextCode[bits] = (int)ucode;
  }

  for (i = 0; i < n; i++) {
    length = lengths[i];
    if (length != 0) {
      codes[i] = context->nextCode[length];
      context->nextCode[length]++;
    }
    else
      codes[i] = 0;
  }

  return 0;
}


static int construct_tree(struct InflateContext *context, struct InflateNode **root,
                          const int *codes, const int *codeLengths, int n) {

  struct InflateNode *node;
  int i, j, code;

  *root = alloc_node(context);
  if (*root == NULL)
    return -1;

  for (i = 0; i < n; i++) {
    if (codeLengths[i] == 0)
      continue;

    node = *root;
    code = codes[i];
    for (j = codeLengths[i] - 1; j >= 0; j--) {
      if (node->literal >= 0)
        return -1;

      if (((code >> j) & 1) == 0) {
        if (node->left == NULL) {
          node->left = alloc_node(context);
          if (node->left == NULL)
            return -1;
        }
        node = node->left;
      }
      else {
        if (node->right == NULL) {
          node->right = alloc_node(context);
          if (node->right == NULL)
            return -1;
        }
        node = node->right;
      }
    }

    if (node->literal >= 0 || node->left != NULL || node->right != NULL)
      return -1;
    node->literal = i;
  }

  return 0;
}


static int decode_symbol(struct bitreader *reader, struct InflateNode *root, int *symbol) {

  struct InflateNode *node;
  int steps, bit;

  if (root == NULL)
    return INFLATE_BAD_CODE;

  node = root;
  steps = 0;
  while (node->literal < 0) {
    if (steps >= INFLATE_MAX_BITS)
      return INFLATE_BAD_CODE;
    if (bit_get(reader, &bit) != 0)
      return INFLATE_TRUNCATED;
    if (bit)
      node = node->right;
    else
      node = node->left;
    if (node == NULL)
      return INFLATE_BAD_CODE;
    steps++;
  }

  *symbol = node->literal;
  return INFLATE_OK;
}


static int sink_put(struct sink *out, int value) {

  unsigned char byte;

  if (out->n >= out->cap)
    return -1;

  byte = (unsigned char)value;
  if (out->vram == NULL) {
    out->bytes[out->n] = byte;
  }
  else if ((out->n & 1U) == 0) {
    out->vram[out->n >> 1] = byte;
  }
  else {
    out->vram[out->n >> 1] =
      (unsigned short)((byte << 8) | out->vram[out->n >> 1]);
  }
  out->n++;
  return 0;
}


static int sink_get(struct sink *out, int index, int *value) {

  unsigned short word;

  if (index >= out->n)
    return -1;

  if (out->vram == NULL) {
    *value = out->bytes[index];
    return 0;
  }

  word = out->vram[index >> 1];
  if (index & 1U)
    *value = (word >> 8) & 0xFF;
  else
    *value = word & 0xFF;
  return 0;
}


static unsigned long read_u32le(const unsigned char *p) {

  return (unsigned long)p[0]
       | ((unsigned long)p[1] << 8)
       | ((unsigned long)p[2] << 16)
       | ((unsigned long)p[3] << 24);
}


static int inflate_stream(const unsigned char *data, int dataSize,
                          struct sink *out, int *outputSize,
                          struct InflateContext *context) {

  struct bitreader reader;
  unsigned long claimedLong;
  int claimed;
  int j, n, m, b, e, length, distance, codesN, bPrevious, rc, sym;

  if (outputSize != NULL)
    *outputSize = 0;

  if (data == NULL || context == NULL || out == NULL || dataSize < 0)
    return INFLATE_BAD_ARGUMENT;
  if (out->bytes == NULL && out->vram == NULL && out->cap != 0)
    return INFLATE_BAD_ARGUMENT;

  if (dataSize < 4)
    return INFLATE_TRUNCATED;
  if (data[0] != 'D' || data[1] != 'E' || data[2] != 'F' || data[3] != 'c')
    return INFLATE_WRONG_HEADER;
  if (dataSize < 9)
    return INFLATE_TRUNCATED;

  claimedLong = read_u32le(data + 4);
  if (claimedLong > (unsigned long)INT_MAX)
    return INFLATE_OUTPUT_OVERFLOW;
  claimed = (int)claimedLong;
  if (claimed > out->cap)
    return INFLATE_OUTPUT_OVERFLOW;
  out->cap = claimed;
  out->n = 0;

  codesN = data[8];
  if (codesN < 4 || codesN > INFLATE_CODELEN_SYMBOLS)
    return INFLATE_BAD_CODE;

  reader.data = data;
  reader.size = dataSize;
  reader.index = 9;
  reader.bit = 7;

  if (bit_getn(&reader, 3, &m) != 0)
    return INFLATE_TRUNCATED;
  if (m < 1 || m > INFLATE_MAX_BITS)
    return INFLATE_BAD_CODE;

  for (j = 0; j < codesN; j++) {
    context->codeLengthCombined[j] = 0;
    for (n = 0; n < m; n++) {
      if (bit_get(&reader, &b) != 0)
        return INFLATE_TRUNCATED;
      context->codeLengthCombined[j] = (context->codeLengthCombined[j] << 1) | b;
    }
    if (context->codeLengthCombined[j] > INFLATE_MAX_BITS)
      return INFLATE_BAD_CODE;
  }
  for (j = codesN; j < INFLATE_CODELEN_SYMBOLS; j++)
    context->codeLengthCombined[j] = 0;

  if (assign_canonical(context, codesN, context->codeLengthCombined,
                       context->codeCombined) != 0)
    return INFLATE_BAD_CODE;

  context->treeNodesFreeCurrent = 0;
  if (construct_tree(context, &context->treeCombined, context->codeCombined,
                     context->codeLengthCombined, codesN) != 0)
    return INFLATE_BAD_CODE;

  for (j = 0; j < INFLATE_LIT_SYMBOLS; j++)
    context->codeLengthLiterals[j] = 0;
  for (j = 0; j < INFLATE_DIST_SYMBOLS; j++)
    context->codeLengthDistances[j] = 0;

  j = 0;
  bPrevious = 0;
  while (j < INFLATE_LIT_SYMBOLS + INFLATE_DIST_SYMBOLS) {
    int count, value;

    rc = decode_symbol(&reader, context->treeCombined, &sym);
    if (rc != INFLATE_OK)
      return rc;
    if (sym < 0 || sym >= codesN)
      return INFLATE_BAD_CODE;

    if (sym <= codesN - 4) {
      if (sym > INFLATE_MAX_BITS)
        return INFLATE_BAD_CODE;
      value = sym;
      count = 1;
    }
    else if (sym == codesN - 3) {
      if (bit_getn(&reader, 2, &e) != 0)
        return INFLATE_TRUNCATED;
      count = e + 3;
      value = bPrevious;
    }
    else if (sym == codesN - 2) {
      if (bit_getn(&reader, 3, &e) != 0)
        return INFLATE_TRUNCATED;
      count = e + 3;
      value = 0;
    }
    else if (sym == codesN - 1) {
      if (bit_getn(&reader, 7, &e) != 0)
        return INFLATE_TRUNCATED;
      count = e + 11;
      value = 0;
    }
    else
      return INFLATE_BAD_CODE;

    if (value < 0 || value > INFLATE_MAX_BITS)
      return INFLATE_BAD_CODE;
    if (count < 1 || j + count > INFLATE_LIT_SYMBOLS + INFLATE_DIST_SYMBOLS)
      return INFLATE_BAD_CODE;

    while (count > 0) {
      if (j < INFLATE_LIT_SYMBOLS)
        context->codeLengthLiterals[j] = value;
      else
        context->codeLengthDistances[j - INFLATE_LIT_SYMBOLS] = value;
      j++;
      count--;
    }
    bPrevious = value;
  }

  if (assign_canonical(context, INFLATE_LIT_SYMBOLS, context->codeLengthLiterals,
                       context->codeLiterals) != 0 ||
      assign_canonical(context, INFLATE_DIST_SYMBOLS, context->codeLengthDistances,
                       context->codeDistances) != 0)
    return INFLATE_BAD_CODE;

  context->treeNodesFreeCurrent = 0;
  if (construct_tree(context, &context->treeLiterals, context->codeLiterals,
                     context->codeLengthLiterals, INFLATE_LIT_SYMBOLS) != 0 ||
      construct_tree(context, &context->treeDistances, context->codeDistances,
                     context->codeLengthDistances, INFLATE_DIST_SYMBOLS) != 0)
    return INFLATE_BAD_CODE;

  for (;;) {
    int src;
    int t;

    rc = decode_symbol(&reader, context->treeLiterals, &sym);
    if (rc != INFLATE_OK)
      return rc;

    if (sym < 256) {
      if (sink_put(out, sym) != 0)
        return INFLATE_SIZE_MISMATCH;
      continue;
    }

    if (sym == 256)
      break;

    if (sym < 257 || sym > 285)
      return INFLATE_BAD_CODE;

    b = sym - 257;
    length = baseValueLengths[b];
    if (extraBitsLengths[b] > 0) {
      if (bit_getn(&reader, extraBitsLengths[b], &e) != 0)
        return INFLATE_TRUNCATED;
      length += e;
    }
    if (length < 3 || length > 258)
      return INFLATE_BAD_CODE;

    rc = decode_symbol(&reader, context->treeDistances, &sym);
    if (rc != INFLATE_OK)
      return rc;
    if (sym < 0 || sym > 29)
      return INFLATE_BAD_CODE;

    distance = baseValueDistances[sym];
    if (extraBitsDistances[sym] > 0) {
      if (bit_getn(&reader, extraBitsDistances[sym], &e) != 0)
        return INFLATE_TRUNCATED;
      distance += e;
    }
    if (distance < 1 || distance > 32768 || distance > out->n)
      return INFLATE_BAD_DISTANCE;
    if (out->n > out->cap || length > out->cap - out->n)
      return INFLATE_SIZE_MISMATCH;

    src = out->n - distance;
    for (t = 0; t < length; t++) {
      int byte;

      if (sink_get(out, src + t, &byte) != 0)
        return INFLATE_BAD_DISTANCE;
      if (sink_put(out, byte) != 0)
        return INFLATE_SIZE_MISMATCH;
    }
  }

  if (out->n != claimed)
    return INFLATE_SIZE_MISMATCH;

  if (outputSize != NULL)
    *outputSize = out->n;
  return INFLATE_OK;
}


int inflate(const unsigned char *data, int dataSize,
            unsigned char *output, int outputCapacity,
            int *outputSize, struct InflateContext *context) {

  struct sink out;

  if (outputCapacity < 0 || (output == NULL && outputCapacity != 0))
    return INFLATE_BAD_ARGUMENT;

  out.bytes = output;
  out.vram = NULL;
  out.n = 0;
  out.cap = outputCapacity;
  return inflate_stream(data, dataSize, &out, outputSize, context);
}


int inflate_vram16(const unsigned char *data, int dataSize,
                   volatile unsigned short *output, int outputBytes,
                   int *outputSize, struct InflateContext *context) {

  struct sink out;

  if (outputBytes < 0 || (output == NULL && outputBytes != 0))
    return INFLATE_BAD_ARGUMENT;
  if (output != NULL && (((unsigned int)(unsigned long)output) & 1U) != 0)
    return INFLATE_BAD_ARGUMENT;

  out.bytes = NULL;
  out.vram = output;
  out.n = 0;
  out.cap = outputBytes;
  return inflate_stream(data, dataSize, &out, outputSize, context);
}


const char *inflate_strerror(int error) {

  switch (error) {
  case INFLATE_OK:
    return "ok";
  case INFLATE_WRONG_HEADER:
    return "wrong header";
  case INFLATE_TRUNCATED:
    return "truncated input";
  case INFLATE_BAD_CODE:
    return "bad Huffman code";
  case INFLATE_OUTPUT_OVERFLOW:
    return "output buffer too small";
  case INFLATE_BAD_DISTANCE:
    return "bad length or distance";
  case INFLATE_SIZE_MISMATCH:
    return "uncompressed size does not match";
  case INFLATE_BAD_ARGUMENT:
    return "bad argument";
  default:
    return "unknown inflate error";
  }
}
