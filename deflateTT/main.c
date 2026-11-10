/*
 * deflateTT, compresses data that can be uncompressed using inflateTT.
 * Literals, lengths and distances follow RFC 1951. The container is the
 * DEFc format described in the repository README: MSB-first bits, one
 * block, and a custom code-length header.
 *
 * Programmed by Ville Helin <vhelin#iki.fi> in 2007.
 * Hash-chained LZ77, 15-bit Huffman limits and I/O checks added in 2026.
 *
 * This code is under GNU General Public Licence (GPL), version 2, June 1991.
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "deflate.h"
#include "huffman.h"


#define LZ_MAX_DIST 32768
#define LZ_MAX_LEN  258
#define LZ_MIN_LEN  3
#define LZ_TOO_FAR  4096
#define LZ_LAZY_LIMIT LZ_MAX_LEN
#define LZ_HASH_BITS 15
#define LZ_HASH_SIZE (1 << LZ_HASH_BITS)

#define LIT_SYMBOLS  286
#define DIST_SYMBOLS 30
#define LENGTH_TOTAL (LIT_SYMBOLS + DIST_SYMBOLS)


static int freqLiterals[LIT_SYMBOLS];
static int freqDistances[DIST_SYMBOLS];
static int freqCombined[119];

static int codeLengthLiterals[LIT_SYMBOLS];
static int codeLengthDistances[DIST_SYMBOLS];
static int codeLengthCombined[119];
static int codeLengths[LENGTH_TOTAL];

static int codeLiterals[LIT_SYMBOLS];
static int codeDistances[DIST_SYMBOLS];
static int codeCombined[119];


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

/* First length that needs the next length code. Index 0 is length 11. */
static const int extraBitsLimitsLength[21] = {
  11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258
};
static const int extraBitsLimitsDistance[27] = {
  5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577, 32769
};


struct bitbuf {
  unsigned char *data;
  int len;
  int cap;
  unsigned int acc;
  int nbits;
};


static int bitbuf_init(struct bitbuf *b, int cap) {

  if (cap < 64)
    cap = 64;
  b->data = (unsigned char *)malloc(cap);
  if (b->data == NULL)
    return -1;
  b->len = 0;
  b->cap = cap;
  b->acc = 0;
  b->nbits = 0;
  return 0;
}


static int bitbuf_grow(struct bitbuf *b, int need) {

  int cap;
  unsigned char *p;

  if (need < 0 || b->len < 0 || b->len > INT_MAX - need)
    return -1;
  if (b->len + need <= b->cap)
    return 0;

  cap = b->cap > 0 ? b->cap : 64;
  while (b->len + need > cap) {
    if (cap > INT_MAX / 2)
      return -1;
    cap *= 2;
  }

  p = (unsigned char *)realloc(b->data, cap);
  if (p == NULL)
    return -1;
  b->data = p;
  b->cap = cap;
  return 0;
}


static int bitbuf_byte(struct bitbuf *b, unsigned char value) {

  if (bitbuf_grow(b, 1) != 0)
    return -1;
  b->data[b->len++] = value;
  return 0;
}


/* MSB first, matching the inflater. A length of 0 writes nothing. */
static int bitbuf_bits(struct bitbuf *b, unsigned int code, int codeLength) {

  int i;

  if (codeLength < 0 || codeLength > 15)
    return -1;

  for (i = 0; i < codeLength; i++) {
    b->acc = (b->acc << 1) | ((code >> (codeLength - i - 1)) & 1U);
    b->nbits++;
    if (b->nbits == 8) {
      if (bitbuf_byte(b, (unsigned char)b->acc) != 0)
        return -1;
      b->acc = 0;
      b->nbits = 0;
    }
  }
  return 0;
}


static int bitbuf_finish(struct bitbuf *b) {

  if (b->nbits == 0)
    return 0;
  return bitbuf_byte(b, (unsigned char)(b->acc << (8 - b->nbits)));
}


static unsigned lz_hash(const unsigned char *p) {

  return (unsigned)(((p[0] << 10) ^ (p[1] << 5) ^ p[2]) & (LZ_HASH_SIZE - 1));
}


static void lz_insert(const unsigned char *data, int fileSize, int *head, int *prev, int i) {

  int h;

  if (i + LZ_MIN_LEN <= fileSize) {
    h = (int)lz_hash(data + i);
    prev[i] = head[h];
    head[h] = i;
  }
}


/* Longest match at i. The chain is newest-first, so an equal length keeps
   the nearer copy. A length-3 match farther than LZ_TOO_FAR costs more
   than three literals, so it is not reported. */
static void lz_find(const unsigned char *data, int fileSize, const int *head,
                    const int *prev, int i, int *bestLenOut, int *bestDistOut) {

  int h, j, n, dist, chain, bestLen, bestDist;

  bestLen = 0;
  bestDist = 0;

  if (i + LZ_MIN_LEN <= fileSize) {
    h = (int)lz_hash(data + i);
    chain = 0;
    for (j = head[h]; j >= 0; j = prev[j]) {
      chain++;
      if (chain > fileSize)
        break;
      dist = i - j;
      if (dist <= 0 || dist > LZ_MAX_DIST)
        break;
      if (i + bestLen < fileSize && data[j + bestLen] != data[i + bestLen])
        continue;
      n = 0;
      while (n < LZ_MAX_LEN && i + n < fileSize && data[j + n] == data[i + n])
        n++;
      if (n > bestLen) {
        bestLen = n;
        bestDist = dist;
        if (bestLen == LZ_MAX_LEN || i + bestLen == fileSize)
          break;
      }
    }
  }

  if (bestLen == LZ_MIN_LEN && bestDist > LZ_TOO_FAR) {
    bestLen = 0;
    bestDist = 0;
  }

  *bestLenOut = bestLen;
  *bestDistOut = bestDist;
}


/* Lazy LZ77. A 3-byte hash chain replaces the old window scan.
   Any match is skipped when the next byte starts a longer one. */
static int lz77_compress(const unsigned char *data, int fileSize, int *lz77,
                         int lz77Cap, int *lz77SizeOut, int *matches, int *dupBytes) {

  int *head, *prev;
  int i, k, bestLen, bestDist, nextLen, nextDist, lz77Size;
  int matchCount, duplicateBytes;

  head = NULL;
  prev = NULL;
  matchCount = 0;
  duplicateBytes = 0;
  lz77Size = 0;

  if (fileSize < 0 || lz77 == NULL || lz77Cap < 1)
    return -1;

  if (fileSize > 0) {
    head = (int *)malloc(sizeof(int) * LZ_HASH_SIZE);
    prev = (int *)malloc(sizeof(int) * fileSize);
    if (head == NULL || prev == NULL) {
      free(head);
      free(prev);
      return -1;
    }
    for (i = 0; i < LZ_HASH_SIZE; i++)
      head[i] = -1;
  }

  i = 0;
  while (i < fileSize) {
    lz_find(data, fileSize, head, prev, i, &bestLen, &bestDist);
    lz_insert(data, fileSize, head, prev, i);

    if (bestLen >= LZ_MIN_LEN && bestLen <= LZ_LAZY_LIMIT &&
        i + 1 + LZ_MIN_LEN <= fileSize) {
      lz_find(data, fileSize, head, prev, i + 1, &nextLen, &nextDist);
      if (nextLen > bestLen && nextDist > 0) {
        if (lz77Size + 1 >= lz77Cap) {
          free(head);
          free(prev);
          return -1;
        }
        lz77[lz77Size++] = data[i];
        i++;
        continue;
      }
    }

    if (bestLen >= LZ_MIN_LEN) {
      if (lz77Size + 2 >= lz77Cap) {
        free(head);
        free(prev);
        return -1;
      }
      lz77[lz77Size++] = bestLen + 254;
      lz77[lz77Size++] = bestDist;
      matchCount++;
      duplicateBytes += bestLen;

      for (k = 1; k < bestLen; k++)
        lz_insert(data, fileSize, head, prev, i + k);
      i += bestLen;
    }
    else {
      if (lz77Size + 1 >= lz77Cap) {
        free(head);
        free(prev);
        return -1;
      }
      lz77[lz77Size++] = data[i];
      i++;
    }
  }

  if (lz77Size + 1 > lz77Cap) {
    free(head);
    free(prev);
    return -1;
  }
  lz77[lz77Size++] = 256;

  free(head);
  free(prev);
  *lz77SizeOut = lz77Size;
  if (matches != NULL)
    *matches = matchCount;
  if (dupBytes != NULL)
    *dupBytes = duplicateBytes;
  return 0;
}


static int preprocess_lz77(int *lz77, int lz77Size) {

  int i, j, n;

  for (i = 0; i < lz77Size; i++) {
    if (lz77[i] < 257)
      continue;

    n = lz77[i] - 254;
    if (n > 10) {
      for (j = 0; j < 20; j++) {
        if (n < extraBitsLimitsLength[j + 1]) {
          lz77[i] = (265 + j) | ((n - extraBitsLimitsLength[j]) << 16);
          break;
        }
      }
      if (j == 20) {
        if (n == 258)
          lz77[i] = 285;
        else {
          fprintf(stderr, "main(): Unsupported length %d in LZ77 preprocess.\n", n);
          return -1;
        }
      }
    }

    i++;
    if (i >= lz77Size) {
      fprintf(stderr, "main(): LZ77 length is missing its distance.\n");
      return -1;
    }

    n = lz77[i];
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
          lz77[i] = (4 + j) | ((n - extraBitsLimitsDistance[j]) << 16);
          break;
        }
      }
      if (j == 26) {
        fprintf(stderr, "main(): Unsupported distance %d in LZ77 preprocess.\n", n);
        return -1;
      }
    }
  }

  return 0;
}


static int rle_code_lengths(int *lengths, int codeLengthMax, int *outCount) {

  int i, j, k, m;

  i = 0;
  j = 0;
  while (i < LENGTH_TOTAL) {
    m = lengths[i];

    if (m < 0 || m > 15) {
      fprintf(stderr, "main(): Got a code length of %d bits.\n", m);
      return -1;
    }

    if (m != 0) {
      if (j >= LENGTH_TOTAL)
        return -1;
      lengths[j++] = m;

      k = 0;
      while (i + k < LENGTH_TOTAL) {
        if (lengths[i + k] != m)
          break;
        if (k == 7)
          break;
        k++;
      }

      i++;
      k--;

      if (k < 3)
        continue;

      if (j >= LENGTH_TOTAL)
        return -1;
      lengths[j++] = (codeLengthMax + 1) | ((k - 3) << 16);
      i += k;
    }
    else {
      k = 1;
      i++;
      while (i < LENGTH_TOTAL && k < 138) {
        if (lengths[i] != 0)
          break;
        k++;
        i++;
      }

      if (k < 3) {
        while (k > 0) {
          if (j >= LENGTH_TOTAL)
            return -1;
          lengths[j++] = 0;
          k--;
        }
      }
      else if (k < 11) {
        if (j >= LENGTH_TOTAL)
          return -1;
        lengths[j++] = (codeLengthMax + 2) | ((k - 3) << 16);
      }
      else {
        if (j >= LENGTH_TOTAL)
          return -1;
        lengths[j++] = (codeLengthMax + 3) | ((k - 11) << 16);
      }
    }
  }

  *outCount = j;
  return 0;
}


int deflate_bytes(const unsigned char *input, int inputSize,
                  unsigned char **output, int *outputSize,
                  int *matches, int *duplicateBytes) {

  int *lz77;
  int lz77Size, i, j, k, m, n;
  int codeLengthMax, codesN, codeLengthsN;
  int matchCount, dupBytes, allocCount;
  struct bitbuf bits;

  lz77 = NULL;
  bits.data = NULL;
  matchCount = 0;
  dupBytes = 0;

  if (output == NULL || outputSize == NULL)
    return -1;
  *output = NULL;
  *outputSize = 0;

  if (inputSize < 0 || (inputSize > 0 && input == NULL)) {
    fprintf(stderr, "main(): Bad input.\n");
    return -1;
  }

  if (inputSize > 0 &&
      (unsigned long)inputSize > ((unsigned long)-1 / sizeof(int)) - 1UL) {
    fprintf(stderr, "main(): Input is too large.\n");
    return -1;
  }

  if (inputSize == INT_MAX) {
    fprintf(stderr, "main(): Input is too large.\n");
    return -1;
  }

  allocCount = inputSize + 1;
  lz77 = (int *)malloc(sizeof(int) * allocCount);
  if (lz77 == NULL) {
    fprintf(stderr, "main(): Out of memory error.\n");
    return -1;
  }

  if (lz77_compress(input, inputSize, lz77, allocCount, &lz77Size,
                    &matchCount, &dupBytes) != 0) {
    fprintf(stderr, "main(): LZ77 failed.\n");
    free(lz77);
    return -1;
  }

  if (matches != NULL)
    *matches = matchCount;
  if (duplicateBytes != NULL)
    *duplicateBytes = dupBytes;

  if (preprocess_lz77(lz77, lz77Size) != 0) {
    free(lz77);
    return -1;
  }

  for (i = 0; i < LIT_SYMBOLS; i++)
    freqLiterals[i] = 0;
  for (i = 0; i < DIST_SYMBOLS; i++)
    freqDistances[i] = 0;

  i = 0;
  while (i < lz77Size) {
    n = lz77[i++] & 0xFFFF;
    if (n < 257) {
      freqLiterals[n]++;
    }
    else {
      if (n > 285) {
        fprintf(stderr, "main(): Bad length symbol %d.\n", n);
        free(lz77);
        return -1;
      }
      freqLiterals[n]++;
      if (i >= lz77Size) {
        fprintf(stderr, "main(): Truncated LZ77 pair.\n");
        free(lz77);
        return -1;
      }
      n = lz77[i++] & 0xFFFF;
      if (n >= DIST_SYMBOLS) {
        fprintf(stderr, "main(): Bad distance symbol %d.\n", n);
        free(lz77);
        return -1;
      }
      freqDistances[n]++;
    }
  }

  if (huffmanBuildCodes(freqLiterals, LIT_SYMBOLS, codeLengthLiterals, codeLiterals) != 0 ||
      huffmanBuildCodes(freqDistances, DIST_SYMBOLS, codeLengthDistances, codeDistances) != 0) {
    fprintf(stderr, "main(): Failed to build Huffman codes.\n");
    free(lz77);
    return -1;
  }

  j = 0;
  for (i = 0; i < LIT_SYMBOLS; i++)
    codeLengths[j++] = codeLengthLiterals[i];
  for (i = 0; i < DIST_SYMBOLS; i++)
    codeLengths[j++] = codeLengthDistances[i];

  codeLengthMax = 0;
  for (i = 0; i < LENGTH_TOTAL; i++) {
    if (codeLengthMax < codeLengths[i])
      codeLengthMax = codeLengths[i];
  }
  codesN = codeLengthMax + 4;
  if (codesN < 4 || codesN > 119) {
    fprintf(stderr, "main(): Code-length alphabet size %d is unusable.\n", codesN);
    free(lz77);
    return -1;
  }

  if (rle_code_lengths(codeLengths, codeLengthMax, &codeLengthsN) != 0) {
    free(lz77);
    return -1;
  }

  for (i = 0; i < codesN; i++)
    freqCombined[i] = 0;

  for (i = 0; i < codeLengthsN; i++) {
    m = codeLengths[i] & 0xFFFF;
    if (m < 0 || m >= codesN) {
      fprintf(stderr, "main(): Code-length symbol %d is out of range.\n", m);
      free(lz77);
      return -1;
    }
    freqCombined[m]++;
  }

  if (huffmanBuildCodes(freqCombined, codesN, codeLengthCombined, codeCombined) != 0) {
    fprintf(stderr, "main(): Failed to build code-length Huffman codes.\n");
    free(lz77);
    return -1;
  }

  if (inputSize > INT_MAX - 64 || bitbuf_init(&bits, inputSize + 64) != 0) {
    fprintf(stderr, "main(): Out of memory error.\n");
    free(lz77);
    return -1;
  }

  if (bitbuf_byte(&bits, 'D') != 0 || bitbuf_byte(&bits, 'E') != 0 ||
      bitbuf_byte(&bits, 'F') != 0 || bitbuf_byte(&bits, 'c') != 0 ||
      bitbuf_byte(&bits, (unsigned char)(inputSize & 0xFF)) != 0 ||
      bitbuf_byte(&bits, (unsigned char)((inputSize >> 8) & 0xFF)) != 0 ||
      bitbuf_byte(&bits, (unsigned char)((inputSize >> 16) & 0xFF)) != 0 ||
      bitbuf_byte(&bits, (unsigned char)((inputSize >> 24) & 0xFF)) != 0 ||
      bitbuf_byte(&bits, (unsigned char)codesN) != 0) {
    fprintf(stderr, "main(): Out of memory error.\n");
    free(bits.data);
    free(lz77);
    return -1;
  }

  n = 0;
  for (i = 0; i < codesN; i++) {
    if (n < codeLengthCombined[i])
      n = codeLengthCombined[i];
  }

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

  if (bitbuf_bits(&bits, (unsigned int)k, 3) != 0) {
    free(bits.data);
    free(lz77);
    return -1;
  }

  for (i = 0; i < codesN; i++) {
    if (codeLengthCombined[i] < 0 || codeLengthCombined[i] >= (1 << k)) {
      fprintf(stderr, "main(): Combined code length entry %d does not fit in %d bits.\n",
              codeLengthCombined[i], k);
      free(bits.data);
      free(lz77);
      return -1;
    }
    if (bitbuf_bits(&bits, (unsigned int)codeLengthCombined[i], k) != 0) {
      free(bits.data);
      free(lz77);
      return -1;
    }
  }

  for (i = 0; i < codeLengthsN; i++) {
    n = codeLengths[i] & 0xFFFF;
    if (n < 0 || n >= codesN || codeLengthCombined[n] <= 0) {
      fprintf(stderr, "main(): Cannot emit code-length symbol %d.\n", n);
      free(bits.data);
      free(lz77);
      return -1;
    }
    if (bitbuf_bits(&bits, (unsigned int)codeCombined[n], codeLengthCombined[n]) != 0) {
      free(bits.data);
      free(lz77);
      return -1;
    }

    if (n > codeLengthMax) {
      int extra, extraLen;

      extra = codeLengths[i] >> 16;
      if (n == codeLengthMax + 1)
        extraLen = 2;
      else if (n == codeLengthMax + 2)
        extraLen = 3;
      else if (n == codeLengthMax + 3)
        extraLen = 7;
      else {
        fprintf(stderr, "main(): Internal error, combined code %d is not supported.\n", n);
        free(bits.data);
        free(lz77);
        return -1;
      }
      if (extra < 0 || extra >= (1 << extraLen) ||
          bitbuf_bits(&bits, (unsigned int)extra, extraLen) != 0) {
        free(bits.data);
        free(lz77);
        return -1;
      }
    }
  }

  i = 0;
  while (i < lz77Size) {
    int extra;

    n = lz77[i] & 0xFFFF;
    if (n < 0 || n >= LIT_SYMBOLS || codeLengthLiterals[n] <= 0 ||
        codeLiterals[n] < 0) {
      fprintf(stderr, "main(): Cannot emit literal symbol %d.\n", n);
      free(bits.data);
      free(lz77);
      return -1;
    }
    if (bitbuf_bits(&bits, (unsigned int)codeLiterals[n], codeLengthLiterals[n]) != 0) {
      free(bits.data);
      free(lz77);
      return -1;
    }

    if (n > 256) {
      if (n > 285) {
        fprintf(stderr, "main(): Bad length symbol %d.\n", n);
        free(bits.data);
        free(lz77);
        return -1;
      }
      k = extraBitsLengths[n - 257];
      if (k > 0) {
        extra = lz77[i] >> 16;
        if (extra < 0 || extra >= (1 << k) ||
            bitbuf_bits(&bits, (unsigned int)extra, k) != 0) {
          free(bits.data);
          free(lz77);
          return -1;
        }
      }
      i++;

      if (i >= lz77Size) {
        fprintf(stderr, "main(): Truncated LZ77 pair while writing.\n");
        free(bits.data);
        free(lz77);
        return -1;
      }

      n = lz77[i] & 0xFFFF;
      if (n < 0 || n >= DIST_SYMBOLS || codeLengthDistances[n] <= 0 ||
          codeDistances[n] < 0) {
        fprintf(stderr, "main(): Cannot emit distance symbol %d.\n", n);
        free(bits.data);
        free(lz77);
        return -1;
      }
      if (bitbuf_bits(&bits, (unsigned int)codeDistances[n], codeLengthDistances[n]) != 0) {
        free(bits.data);
        free(lz77);
        return -1;
      }

      k = extraBitsDistances[n];
      if (k > 0) {
        extra = lz77[i] >> 16;
        if (extra < 0 || extra >= (1 << k) ||
            bitbuf_bits(&bits, (unsigned int)extra, k) != 0) {
          free(bits.data);
          free(lz77);
          return -1;
        }
      }
      i++;
    }
    else
      i++;
  }

  if (bitbuf_finish(&bits) != 0) {
    free(bits.data);
    free(lz77);
    return -1;
  }

  free(lz77);

  *output = bits.data;
  *outputSize = bits.len;
  return 0;
}


#ifndef DEFLATE_LIBRARY
int main(int argc, char *argv[]) {

  int fileSize, outSize, matches, dupBytes;
  long fileSizeLong;
  unsigned char *data, *compressed;
  FILE *f;

  data = NULL;
  compressed = NULL;

  if (argc != 3) {
    fprintf(stderr, "deflateTT v2.2 Written by Ville Helin 2007\n");
    fprintf(stderr, "USAGE: %s <INPUT FILE> <OUTPUT FILE>\n", argv[0]);
    return 1;
  }

  f = fopen(argv[1], "rb");
  if (f == NULL) {
    fprintf(stderr, "main(): Could not open file \"%s\" for reading.\n", argv[1]);
    return 1;
  }

  if (fseek(f, 0, SEEK_END) != 0) {
    fprintf(stderr, "main(): Could not seek in \"%s\".\n", argv[1]);
    fclose(f);
    return 1;
  }
  fileSizeLong = ftell(f);
  if (fileSizeLong < 0 || fileSizeLong > (long)INT_MAX) {
    fprintf(stderr, "main(): Unsupported file size in \"%s\".\n", argv[1]);
    fclose(f);
    return 1;
  }
  if (fseek(f, 0, SEEK_SET) != 0) {
    fprintf(stderr, "main(): Could not seek in \"%s\".\n", argv[1]);
    fclose(f);
    return 1;
  }
  fileSize = (int)fileSizeLong;

  data = (unsigned char *)malloc(fileSize > 0 ? fileSize : 1);
  if (data == NULL) {
    fprintf(stderr, "main(): Out of memory error.\n");
    fclose(f);
    return 1;
  }

  if (fileSize > 0 && (int)fread(data, 1, fileSize, f) != fileSize) {
    fprintf(stderr, "main(): Could not read \"%s\".\n", argv[1]);
    free(data);
    fclose(f);
    return 1;
  }
  fclose(f);

  if (deflate_bytes(data, fileSize, &compressed, &outSize, &matches, &dupBytes) != 0) {
    free(data);
    return 1;
  }

  fprintf(stderr, "main(): LZ77: %d utilized matches | %d duplicate bytes.\n",
          matches, dupBytes);

  f = fopen(argv[2], "wb");
  if (f == NULL) {
    fprintf(stderr, "main(): Could not open file \"%s\" for writing.\n", argv[2]);
    free(data);
    free(compressed);
    return 1;
  }

  if ((int)fwrite(compressed, 1, outSize, f) != outSize) {
    fprintf(stderr, "main(): Could not write \"%s\".\n", argv[2]);
    fclose(f);
    free(data);
    free(compressed);
    return 1;
  }
  fclose(f);

  if (fileSize == 0) {
    fprintf(stderr, "main(): Original size = 0B, deflated size = %dB.\n", outSize);
  }
  else {
    fprintf(stderr,
            "main(): Original size = %dB, deflated size = %dB -> Got rid of %.2f%%.\n",
            fileSize, outSize, 100.0 - (outSize * 100.0 / fileSize));
  }

  free(data);
  free(compressed);
  return 0;
}
#endif
