/*
 * deflateTT, compresses data that can be uncompressed using inflateTT.
 * Literals, lengths and distances follow RFC 1951. The container is the
 * DEFd format described in the repository README: MSB-first bits, one
 * or more Huffman or stored blocks, and a custom code-length header.
 * The window crosses blocks.
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


/* MSB first, matching the inflater. A length of 0 writes nothing.
   The buffered bits are the low nbits of acc. */
static int bitbuf_bits(struct bitbuf *b, unsigned int code, int codeLength) {

  if (codeLength < 0 || codeLength > 15)
    return -1;
  if (codeLength == 0)
    return 0;

  b->acc = (b->acc << codeLength) | (code & ((1U << codeLength) - 1U));
  b->nbits += codeLength;
  while (b->nbits >= 8) {
    b->nbits -= 8;
    if (bitbuf_byte(b, (unsigned char)(b->acc >> b->nbits)) != 0)
      return -1;
  }
  if (b->nbits == 0)
    b->acc = 0;
  else
    b->acc &= (1U << b->nbits) - 1U;
  return 0;
}


static int bitbuf_finish(struct bitbuf *b) {

  int rc;

  if (b->nbits == 0)
    return 0;
  rc = bitbuf_byte(b, (unsigned char)(b->acc << (8 - b->nbits)));
  b->acc = 0;
  b->nbits = 0;
  return rc;
}


static int bitbuf_bytes(struct bitbuf *b, const unsigned char *bytes, int n) {

  int i;

  if (b->nbits != 0 || n < 0)
    return -1;
  if (bitbuf_grow(b, n) != 0)
    return -1;
  i = 0;
  while (i + 4 <= n) {
    b->data[b->len] = bytes[i];
    b->data[b->len + 1] = bytes[i + 1];
    b->data[b->len + 2] = bytes[i + 2];
    b->data[b->len + 3] = bytes[i + 3];
    b->len += 4;
    i += 4;
  }
  while (i < n)
    b->data[b->len++] = bytes[i++];
  return 0;
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

  int h, j, n, dist, chain, bestLen, bestDist, limit;

  bestLen = 0;
  bestDist = 0;
  limit = 0;

  if (i + LZ_MIN_LEN <= fileSize) {
    limit = fileSize - i;
    if (limit > LZ_MAX_LEN)
      limit = LZ_MAX_LEN;
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
      while (n + 4 <= limit &&
             data[j + n] == data[i + n] &&
             data[j + n + 1] == data[i + n + 1] &&
             data[j + n + 2] == data[i + n + 2] &&
             data[j + n + 3] == data[i + n + 3])
        n += 4;
      while (n < limit && data[j + n] == data[i + n])
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
   A match shorter than 258 bytes is skipped when the next byte starts
   a longer one. A match of 258 bytes cannot be beaten by the next byte. */
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

    if (bestLen >= LZ_MIN_LEN && bestLen < LZ_MAX_LEN &&
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


/* Output bytes covered by one preprocessed LZ77 token. Distances are skipped. */
static int span_bytes(const int *lz77, int lz77Size, int i, int *step) {

  int sym, extra, j;

  if (i < 0 || i >= lz77Size)
    return -1;

  sym = lz77[i] & 0xFFFF;
  extra = lz77[i] >> 16;
  if (extra < 0)
    return -1;

  if (sym < 256) {
    *step = 1;
    return 1;
  }
  if (sym == 256) {
    *step = 1;
    return 0;
  }
  if (sym > 285 || i + 1 >= lz77Size)
    return -1;

  *step = 2;
  if (sym <= 264)
    return sym - 254;
  if (sym == 285)
    return 258;
  j = sym - 265;
  if (j < 0 || j >= 20)
    return -1;
  return extraBitsLimitsLength[j] + extra;
}


/* Running size of a block that is not written out. bytes counts finished
   bytes. nbits is the unfinished tail, 0..7, matching bitbuf_finish. */
struct bitcount {
  int bytes;
  int nbits;
};


static int count_bits(struct bitcount *count, int n) {

  int room, take;

  if (n < 0)
    return -1;
  while (n > 0) {
    room = 8 - count->nbits;
    take = n < room ? n : room;
    count->nbits += take;
    n -= take;
    if (count->nbits == 8) {
      if (count->bytes == INT_MAX)
        return -1;
      count->bytes++;
      count->nbits = 0;
    }
  }
  return 0;
}


static int count_size(const struct bitcount *count) {

  if (count->nbits == 0)
    return count->bytes;
  if (count->bytes == INT_MAX)
    return -1;
  return count->bytes + 1;
}


static int put_bits(struct bitbuf *bits, struct bitcount *count,
                    unsigned int code, int codeLength) {

  if (bits == NULL)
    return count_bits(count, codeLength);
  return bitbuf_bits(bits, code, codeLength);
}


static int put_byte(struct bitbuf *bits, struct bitcount *count, unsigned char value) {

  if (bits == NULL) {
    if (count->nbits != 0)
      return -1;
    return count_bits(count, 8);
  }
  return bitbuf_byte(bits, value);
}


/* Huffman-encode tokens [begin, end) and an end-of-block symbol.
   The block starts with the code-length alphabet size, on a byte boundary.
   A NULL bit buffer counts the same bytes without writing them. */
static int emit_huffman_range(struct bitbuf *bits, struct bitcount *counted,
                              const int *lz77, int begin, int end) {

  int i, j, k, m, n;
  int codeLengthMax, codesN, codeLengthsN;

  if (begin < 0 || end < begin)
    return -1;
  if (bits == NULL) {
    if (counted == NULL)
      return -1;
  }
  else if (bits->nbits != 0)
    return -1;

  for (i = 0; i < LIT_SYMBOLS; i++)
    freqLiterals[i] = 0;
  for (i = 0; i < DIST_SYMBOLS; i++)
    freqDistances[i] = 0;
  freqLiterals[256]++;

  i = begin;
  while (i < end) {
    n = lz77[i++] & 0xFFFF;
    if (n < 257) {
      freqLiterals[n]++;
    }
    else {
      if (n > 285) {
        fprintf(stderr, "main(): Bad length symbol %d.\n", n);
        return -1;
      }
      freqLiterals[n]++;
      if (i >= end) {
        fprintf(stderr, "main(): Truncated LZ77 pair.\n");
        return -1;
      }
      n = lz77[i++] & 0xFFFF;
      if (n >= DIST_SYMBOLS) {
        fprintf(stderr, "main(): Bad distance symbol %d.\n", n);
        return -1;
      }
      freqDistances[n]++;
    }
  }

  if (huffmanBuildCodes(freqLiterals, LIT_SYMBOLS, codeLengthLiterals, codeLiterals) != 0 ||
      huffmanBuildCodes(freqDistances, DIST_SYMBOLS, codeLengthDistances, codeDistances) != 0) {
    fprintf(stderr, "main(): Failed to build Huffman codes.\n");
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
    return -1;
  }

  if (rle_code_lengths(codeLengths, codeLengthMax, &codeLengthsN) != 0)
    return -1;

  for (i = 0; i < codesN; i++)
    freqCombined[i] = 0;

  for (i = 0; i < codeLengthsN; i++) {
    m = codeLengths[i] & 0xFFFF;
    if (m < 0 || m >= codesN) {
      fprintf(stderr, "main(): Code-length symbol %d is out of range.\n", m);
      return -1;
    }
    freqCombined[m]++;
  }

  if (huffmanBuildCodes(freqCombined, codesN, codeLengthCombined, codeCombined) != 0) {
    fprintf(stderr, "main(): Failed to build code-length Huffman codes.\n");
    return -1;
  }

  if (put_byte(bits, counted, (unsigned char)codesN) != 0)
    return -1;

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

  if (put_bits(bits, counted, (unsigned int)k, 3) != 0)
    return -1;

  for (i = 0; i < codesN; i++) {
    if (codeLengthCombined[i] < 0 || codeLengthCombined[i] >= (1 << k)) {
      fprintf(stderr, "main(): Combined code length entry %d does not fit in %d bits.\n",
              codeLengthCombined[i], k);
      return -1;
    }
    if (put_bits(bits, counted, (unsigned int)codeLengthCombined[i], k) != 0)
      return -1;
  }

  for (i = 0; i < codeLengthsN; i++) {
    int extra, extraLen;

    n = codeLengths[i] & 0xFFFF;
    if (n < 0 || n >= codesN || codeLengthCombined[n] <= 0) {
      fprintf(stderr, "main(): Cannot emit code-length symbol %d.\n", n);
      return -1;
    }
    if (put_bits(bits, counted, (unsigned int)codeCombined[n], codeLengthCombined[n]) != 0)
      return -1;

    if (n > codeLengthMax) {
      extra = codeLengths[i] >> 16;
      if (n == codeLengthMax + 1)
        extraLen = 2;
      else if (n == codeLengthMax + 2)
        extraLen = 3;
      else if (n == codeLengthMax + 3)
        extraLen = 7;
      else {
        fprintf(stderr, "main(): Internal error, combined code %d is not supported.\n", n);
        return -1;
      }
      if (extra < 0 || extra >= (1 << extraLen) ||
          put_bits(bits, counted, (unsigned int)extra, extraLen) != 0)
        return -1;
    }
  }

  i = begin;
  while (i < end) {
    int extra;

    n = lz77[i] & 0xFFFF;
    if (n < 0 || n >= LIT_SYMBOLS || codeLengthLiterals[n] <= 0 ||
        codeLiterals[n] < 0) {
      fprintf(stderr, "main(): Cannot emit literal symbol %d.\n", n);
      return -1;
    }
    if (put_bits(bits, counted, (unsigned int)codeLiterals[n], codeLengthLiterals[n]) != 0)
      return -1;

    if (n > 256) {
      if (n > 285) {
        fprintf(stderr, "main(): Bad length symbol %d.\n", n);
        return -1;
      }
      k = extraBitsLengths[n - 257];
      if (k > 0) {
        extra = lz77[i] >> 16;
        if (extra < 0 || extra >= (1 << k) ||
            put_bits(bits, counted, (unsigned int)extra, k) != 0)
          return -1;
      }
      i++;

      if (i >= end) {
        fprintf(stderr, "main(): Truncated LZ77 pair while writing.\n");
        return -1;
      }

      n = lz77[i] & 0xFFFF;
      if (n < 0 || n >= DIST_SYMBOLS || codeLengthDistances[n] <= 0 ||
          codeDistances[n] < 0) {
        fprintf(stderr, "main(): Cannot emit distance symbol %d.\n", n);
        return -1;
      }
      if (put_bits(bits, counted, (unsigned int)codeDistances[n], codeLengthDistances[n]) != 0)
        return -1;

      k = extraBitsDistances[n];
      if (k > 0) {
        extra = lz77[i] >> 16;
        if (extra < 0 || extra >= (1 << k) ||
            put_bits(bits, counted, (unsigned int)extra, k) != 0)
          return -1;
      }
      i++;
    }
    else
      i++;
  }

  n = 256;
  if (codeLengthLiterals[n] <= 0 || codeLiterals[n] < 0) {
    fprintf(stderr, "main(): Cannot emit literal symbol %d.\n", n);
    return -1;
  }
  if (put_bits(bits, counted, (unsigned int)codeLiterals[n], codeLengthLiterals[n]) != 0)
    return -1;
  if (bits == NULL)
    return 0;
  return bitbuf_finish(bits);
}


static int emit_stored(struct bitbuf *bits, const unsigned char *raw, int n) {

  if (bits->nbits != 0 || n < 1)
    return -1;
  if (bitbuf_byte(bits, 0) != 0 ||
      bitbuf_byte(bits, (unsigned char)(n & 0xFF)) != 0 ||
      bitbuf_byte(bits, (unsigned char)((n >> 8) & 0xFF)) != 0 ||
      bitbuf_byte(bits, (unsigned char)((n >> 16) & 0xFF)) != 0 ||
      bitbuf_byte(bits, (unsigned char)((n >> 24) & 0xFF)) != 0)
    return -1;
  return bitbuf_bytes(bits, raw, n);
}


static int huffman_cost(const int *lz77, int tok0, int tok1) {

  struct bitcount counted;
  int rc;

  counted.bytes = 0;
  counted.nbits = 0;
  rc = emit_huffman_range(NULL, &counted, lz77, tok0, tok1);
  if (rc != 0)
    return -1;
  return count_size(&counted);
}


static int build_bounds(const int *lz77, int lz77Size, int inputSize,
                        int *tokAt, int *byteAt, int *nboundsOut) {

  int i, step, bytes, pos, nbounds;

  nbounds = 0;
  tokAt[nbounds] = 0;
  byteAt[nbounds] = 0;
  nbounds++;

  i = 0;
  pos = 0;
  while (i < lz77Size && (lz77[i] & 0xFFFF) != 256) {
    bytes = span_bytes(lz77, lz77Size, i, &step);
    if (bytes < 0 || step < 1 || pos > inputSize - bytes)
      return -1;
    i += step;
    pos += bytes;
    if (nbounds >= lz77Size)
      return -1;
    tokAt[nbounds] = i;
    byteAt[nbounds] = pos;
    nbounds++;
  }

  if (i >= lz77Size || (lz77[i] & 0xFFFF) != 256 || pos != inputSize) {
    fprintf(stderr, "main(): LZ77 token stream does not cover the input.\n");
    return -1;
  }
  *nboundsOut = nbounds;
  return 0;
}


static int histogram_score(const unsigned char *data, int n, int pos, int window) {

  int a[256], b[256];
  int i, d, score;

  if (window < 1 || pos < window || pos > n - window)
    return 0;
  for (i = 0; i < 256; i++) {
    a[i] = 0;
    b[i] = 0;
  }
  for (i = pos - window; i < pos; i++)
    a[data[i]]++;
  for (i = pos; i < pos + window; i++)
    b[data[i]]++;
  score = 0;
  for (i = 0; i < 256; i++) {
    d = a[i] - b[i];
    if (d < 0)
      d = -d;
    score += d;
  }
  return score;
}


/* Up to two byte positions where the neighbouring histograms disagree.
   A shift has to be stronger than one count per symbol in the window,
   and the two positions have to be at least one window apart. */
static void find_splits(const unsigned char *data, int n, int *s0, int *s1, int *nfound) {

  int window, step, pos, score, best1, best2, score1, score2, gap, have;

  *nfound = 0;
  *s0 = 0;
  *s1 = 0;
  if (n < 4096)
    return;

  window = 1024;
  step = 256;
  best1 = 0;
  best2 = 0;
  score1 = 0;
  score2 = 0;
  have = 0;

  for (pos = window; pos <= n - window; pos += step) {
    score = histogram_score(data, n, pos, window);
    if (score <= window)
      continue;
    if (have == 0 || score > score1) {
      if (have >= 1) {
        gap = pos - best1;
        if (gap < 0)
          gap = -gap;
        if (gap >= window) {
          best2 = best1;
          score2 = score1;
          have = 2;
        }
      }
      best1 = pos;
      score1 = score;
      if (have == 0)
        have = 1;
      else if (have == 2) {
        gap = best2 - best1;
        if (gap < 0)
          gap = -gap;
        if (gap < window)
          have = 1;
      }
    }
    else {
      gap = pos - best1;
      if (gap < 0)
        gap = -gap;
      if (gap >= window && score > score2) {
        best2 = pos;
        score2 = score;
        have = 2;
      }
    }
  }

  if (have >= 1) {
    *s0 = best1;
    *nfound = 1;
  }
  if (have >= 2) {
    *s1 = best2;
    *nfound = 2;
  }
}


static int snap_bound(const int *byteAt, int nbounds, int target) {

  int i, best, dist, bestd;

  if (nbounds < 3)
    return -1;
  best = -1;
  bestd = 0;
  for (i = 1; i < nbounds - 1; i++) {
    dist = byteAt[i] - target;
    if (dist < 0)
      dist = -dist;
    if (best < 0 || dist < bestd) {
      best = i;
      bestd = dist;
    }
  }
  return best;
}


struct blockseg {
  int tok0;
  int tok1;
  int byte0;
  int byte1;
  int stored;
};


static int segment_choice(const int *lz77, int tok0, int tok1,
                          int byte0, int byte1, int *storedOut, int *costOut) {

  int raw, storedCost, huff;

  raw = byte1 - byte0;
  if (raw < 0 || tok0 > tok1)
    return -1;
  huff = huffman_cost(lz77, tok0, tok1);
  if (huff < 0)
    return -1;
  if (raw > 0 && raw <= INT_MAX - 5) {
    storedCost = 5 + raw;
    if (storedCost < huff) {
      *storedOut = 1;
      *costOut = storedCost;
      return 0;
    }
  }
  *storedOut = 0;
  *costOut = huff;
  return 0;
}


static int consider_blocks(const int *lz77, const int *tokAt, const int *byteAt,
                           int i0, int i1, int i2, int i3, int nparts,
                           struct blockseg *best, int *bestN, int *bestCost) {

  struct blockseg part[3];
  int idx[4];
  int p, total, stored, partCost, raw;

  idx[0] = i0;
  idx[1] = i1;
  idx[2] = i2;
  idx[3] = i3;
  total = 0;
  for (p = 0; p < nparts; p++) {
    raw = byteAt[idx[p + 1]] - byteAt[idx[p]];
    if (raw == 0 && nparts > 1)
      return 0;
    if (segment_choice(lz77, tokAt[idx[p]], tokAt[idx[p + 1]],
                       byteAt[idx[p]], byteAt[idx[p + 1]],
                       &stored, &partCost) != 0)
      return -1;
    if (partCost < 0 || partCost > INT_MAX - total)
      return -1;
    total += partCost;
    part[p].tok0 = tokAt[idx[p]];
    part[p].tok1 = tokAt[idx[p + 1]];
    part[p].byte0 = byteAt[idx[p]];
    part[p].byte1 = byteAt[idx[p + 1]];
    part[p].stored = stored;
  }
  if (total < *bestCost) {
    *bestCost = total;
    *bestN = nparts;
    for (p = 0; p < nparts; p++)
      best[p] = part[p];
  }
  return 0;
}


static int plan_blocks(const unsigned char *input, int inputSize,
                       const int *lz77, const int *tokAt, const int *byteAt,
                       int nbounds, struct blockseg *plan, int *planCount) {

  struct blockseg best[3];
  int bestN, bestCost, nf, s0, s1, a, b, lo, hi;

  bestN = 0;
  bestCost = INT_MAX;
  a = -1;
  b = -1;
  if (consider_blocks(lz77, tokAt, byteAt, 0, nbounds - 1, 0, 0, 1,
                      best, &bestN, &bestCost) != 0)
    return -1;

  find_splits(input, inputSize, &s0, &s1, &nf);
  if (nf >= 1) {
    a = snap_bound(byteAt, nbounds, s0);
    if (a > 0 &&
        consider_blocks(lz77, tokAt, byteAt, 0, a, nbounds - 1, 0, 2,
                        best, &bestN, &bestCost) != 0)
      return -1;
  }
  if (nf >= 2) {
    b = snap_bound(byteAt, nbounds, s1);
    if (b > 0 && b != a &&
        consider_blocks(lz77, tokAt, byteAt, 0, b, nbounds - 1, 0, 2,
                        best, &bestN, &bestCost) != 0)
      return -1;
    if (a > 0 && b > 0 && a != b) {
      if (a < b) {
        lo = a;
        hi = b;
      }
      else {
        lo = b;
        hi = a;
      }
      if (consider_blocks(lz77, tokAt, byteAt, 0, lo, hi, nbounds - 1, 3,
                          best, &bestN, &bestCost) != 0)
        return -1;
    }
  }

  if (bestN < 1)
    return -1;
  for (lo = 0; lo < bestN; lo++)
    plan[lo] = best[lo];
  *planCount = bestN;
  return 0;
}


int deflate_bytes(const unsigned char *input, int inputSize,
                  unsigned char **output, int *outputSize,
                  int *matches, int *duplicateBytes) {

  int *lz77;
  int *tokAt;
  int *byteAt;
  struct blockseg segs[3];
  int lz77Size, i, nbounds, nsegs;
  int matchCount, dupBytes, allocCount;
  struct bitbuf bits;

  lz77 = NULL;
  tokAt = NULL;
  byteAt = NULL;
  bits.data = NULL;
  matchCount = 0;
  dupBytes = 0;
  nbounds = 0;
  nsegs = 0;

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

  tokAt = (int *)malloc(sizeof(int) * lz77Size);
  byteAt = (int *)malloc(sizeof(int) * lz77Size);
  if (tokAt == NULL || byteAt == NULL) {
    fprintf(stderr, "main(): Out of memory error.\n");
    free(tokAt);
    free(byteAt);
    free(lz77);
    return -1;
  }

  if (build_bounds(lz77, lz77Size, inputSize, tokAt, byteAt, &nbounds) != 0 ||
      plan_blocks(input, inputSize, lz77, tokAt, byteAt, nbounds, segs, &nsegs) != 0) {
    free(tokAt);
    free(byteAt);
    free(lz77);
    return -1;
  }

  if (inputSize > INT_MAX - 64 || bitbuf_init(&bits, inputSize + 64) != 0) {
    fprintf(stderr, "main(): Out of memory error.\n");
    free(tokAt);
    free(byteAt);
    free(lz77);
    return -1;
  }

  if (bitbuf_byte(&bits, 'D') != 0 || bitbuf_byte(&bits, 'E') != 0 ||
      bitbuf_byte(&bits, 'F') != 0 || bitbuf_byte(&bits, 'd') != 0 ||
      bitbuf_byte(&bits, (unsigned char)(inputSize & 0xFF)) != 0 ||
      bitbuf_byte(&bits, (unsigned char)((inputSize >> 8) & 0xFF)) != 0 ||
      bitbuf_byte(&bits, (unsigned char)((inputSize >> 16) & 0xFF)) != 0 ||
      bitbuf_byte(&bits, (unsigned char)((inputSize >> 24) & 0xFF)) != 0) {
    fprintf(stderr, "main(): Out of memory error.\n");
    free(bits.data);
    free(tokAt);
    free(byteAt);
    free(lz77);
    return -1;
  }

  for (i = 0; i < nsegs; i++) {
    if (segs[i].stored) {
      if (emit_stored(&bits, input + segs[i].byte0,
                      segs[i].byte1 - segs[i].byte0) != 0) {
        free(bits.data);
        free(tokAt);
        free(byteAt);
        free(lz77);
        return -1;
      }
    }
    else if (emit_huffman_range(&bits, NULL, lz77, segs[i].tok0, segs[i].tok1) != 0) {
      free(bits.data);
      free(tokAt);
      free(byteAt);
      free(lz77);
      return -1;
    }
  }

  free(tokAt);
  free(byteAt);
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
    fprintf(stderr, "deflateTT v2.3 Written by Ville Helin 2007\n");
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
