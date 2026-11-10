/*
 * Round-trip and error tests for deflateTT / inflateTT-MP.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "deflate.h"
#include "huffman.h"
#include "inflate.h"

static int g_fail;

static void fail(const char *name, const char *why) {

  fprintf(stderr, "FAIL %s: %s\n", name, why);
  g_fail++;
}


static int codes_well_formed(const int *lengths, const int *codes, int n) {

  int i, j;

  for (i = 0; i < n; i++) {
    if (lengths[i] < 0 || lengths[i] > 15)
      return 0;
    if (lengths[i] == 0)
      continue;
    if (codes[i] < 0 || codes[i] >= (1 << lengths[i]))
      return 0;
    for (j = i + 1; j < n; j++) {
      if (lengths[j] == lengths[i] && codes[j] <= codes[i])
        return 0;
    }
  }
  return 1;
}


static int kraft_ok(const int *lengths, int n) {

  int count[16];
  int i, bits, left;

  for (i = 0; i <= 15; i++)
    count[i] = 0;
  for (i = 0; i < n; i++) {
    if (lengths[i] < 0 || lengths[i] > 15)
      return 0;
    count[lengths[i]]++;
  }
  left = 1;
  for (bits = 1; bits <= 15; bits++) {
    left <<= 1;
    left -= count[bits];
    if (left < 0)
      return 0;
  }
  return 1;
}


static void check_huffman(void) {

  int freq[40], lengths[40], codes[40];
  int i, heavy, light, rc;

  for (i = 0; i < 40; i++)
    freq[i] = 0;

  rc = huffmanBuildCodes(NULL, 4, lengths, codes);
  if (rc == 0)
    fail("huffman null", "accepted a NULL frequency table");

  rc = huffmanBuildCodes(freq, 0, lengths, codes);
  if (rc == 0)
    fail("huffman n=0", "accepted an empty alphabet");

  freq[0] = -1;
  rc = huffmanBuildCodes(freq, 4, lengths, codes);
  if (rc == 0)
    fail("huffman negative", "accepted a negative frequency");
  freq[0] = 0;

  rc = huffmanBuildCodes(freq, 4, lengths, codes);
  if (rc != 0 || lengths[0] != 1 || lengths[1] != 1 ||
      lengths[2] != 0 || lengths[3] != 0 || codes[0] != 0 || codes[1] != 1)
    fail("huffman empty", "expected two dummy length-1 codes");
  else if (!kraft_ok(lengths, 4))
    fail("huffman empty", "Kraft inequality failed");

  freq[4] = 9;
  rc = huffmanBuildCodes(freq, 8, lengths, codes);
  if (rc != 0 || lengths[0] != 1 || lengths[4] != 1 || codes[0] != 0 || codes[4] != 1)
    fail("huffman one symbol", "expected the used symbol and one dummy at length 1");
  else {
    for (i = 0; i < 8; i++) {
      if (i != 0 && i != 4 && lengths[i] != 0)
        fail("huffman one symbol", "unexpected extra code");
    }
  }
  freq[4] = 0;

  freq[0] = 3;
  rc = huffmanBuildCodes(freq, 1, lengths, codes);
  if (rc != 0 || lengths[0] != 1 || codes[0] != 0)
    fail("huffman n=1", "expected a single length-1 code");
  freq[0] = 0;

  freq[0] = 1;
  freq[1] = 1;
  for (i = 2; i < 32; i++)
    freq[i] = freq[i - 1] + freq[i - 2];
  rc = huffmanBuildCodes(freq, 32, lengths, codes);
  if (rc != 0)
    fail("huffman fibonacci", "length-limited build failed");
  else if (!kraft_ok(lengths, 32))
    fail("huffman fibonacci", "Kraft inequality failed");
  else {
    heavy = 0;
    light = 0;
    for (i = 1; i < 32; i++) {
      if (freq[i] > freq[heavy])
        heavy = i;
      if (freq[i] < freq[light])
        light = i;
    }
    if (lengths[heavy] > lengths[light])
      fail("huffman fibonacci", "a rarer symbol got a shorter code");
    for (i = 0; i < 32; i++) {
      int j;
      if (lengths[i] <= 0)
        fail("huffman fibonacci", "a used symbol was given length 0");
      for (j = i + 1; j < 32; j++) {
        if (lengths[i] == lengths[j] && codes[i] == codes[j])
          fail("huffman fibonacci", "duplicate canonical code");
      }
    }
    if (!codes_well_formed(lengths, codes, 32))
      fail("huffman fibonacci", "codes are not canonical");
  }
}


static void check_huffman_limits(void) {

  int *freq, *lengths, *codes;
  int rc, i;

  freq = (int *)malloc(sizeof(int) * 4096);
  lengths = (int *)malloc(sizeof(int) * 4096);
  codes = (int *)malloc(sizeof(int) * 4096);
  if (freq == NULL || lengths == NULL || codes == NULL) {
    fail("huffman limits", "out of memory");
    free(freq);
    free(lengths);
    free(codes);
    return;
  }

  for (i = 0; i < 4096; i++)
    freq[i] = 0;

  rc = huffmanBuildCodes(freq, 4097, lengths, codes);
  if (rc == 0)
    fail("huffman n=4097", "accepted an oversized alphabet");

  rc = huffmanBuildCodes(freq, 4, NULL, codes);
  if (rc == 0)
    fail("huffman null lengths", "accepted NULL code lengths");

  rc = huffmanBuildCodes(freq, 4, lengths, NULL);
  if (rc == 0)
    fail("huffman null codes", "accepted NULL codes");

  rc = huffmanBuildCodes(freq, 4096, lengths, codes);
  if (rc != 0 || lengths[0] != 1 || lengths[1] != 1 || lengths[4095] != 0 ||
      codes[0] != 0 || codes[1] != 1)
    fail("huffman n=4096", "expected two dummy length-1 codes");
  else if (!kraft_ok(lengths, 4096) || !codes_well_formed(lengths, codes, 4096))
    fail("huffman n=4096", "codes are not a complete 15-bit set");

  freq[0] = 1;
  freq[1] = 1;
  freq[2] = 1;
  freq[3] = 1;
  rc = huffmanBuildCodes(freq, 4, lengths, codes);
  if (rc != 0 || lengths[0] != 2 || lengths[1] != 2 || lengths[2] != 2 || lengths[3] != 2 ||
      codes[0] != 0 || codes[1] != 1 || codes[2] != 2 || codes[3] != 3)
    fail("huffman equal", "expected codes 0..3 at length 2");

  freq[0] = 100;
  freq[1] = 0;
  freq[2] = 1;
  freq[3] = 0;
  rc = huffmanBuildCodes(freq, 4, lengths, codes);
  if (rc != 0 || lengths[1] != 0 || lengths[3] != 0 || lengths[0] < 1 || lengths[2] < 1 ||
      lengths[0] > lengths[2] || codes[0] == codes[2])
    fail("huffman skewed", "a heavier symbol got a longer code");
  else if (!codes_well_formed(lengths, codes, 4))
    fail("huffman skewed", "codes are not canonical");

  free(freq);
  free(lengths);
  free(codes);
}


static void check_roundtrip(const char *name, const unsigned char *in, int n) {

  unsigned char *comp;
  unsigned char *out;
  unsigned short *vram;
  struct InflateContext *ctx;
  int compSize, rc, matches, stored;
  int produced, i;

  comp = NULL;
  out = NULL;
  vram = NULL;
  ctx = NULL;
  compSize = 0;

  rc = deflate_bytes(in, n, &comp, &compSize, &matches, NULL);
  if (rc != 0 || comp == NULL || compSize < 9) {
    fail(name, "deflate failed");
    free(comp);
    return;
  }
  if (comp[0] != 'D' || comp[1] != 'E' || comp[2] != 'F' || comp[3] != 'd')
    fail(name, "missing DEFd magic");
  stored = (int)comp[4] | ((int)comp[5] << 8) |
           ((int)comp[6] << 16) | ((int)comp[7] << 24);
  if (stored != n)
    fail(name, "header size does not match");

  ctx = (struct InflateContext *)malloc(sizeof(struct InflateContext));
  out = (unsigned char *)malloc(n > 0 ? n : 1);
  if (ctx == NULL || out == NULL) {
    fail(name, "out of memory");
    free(comp);
    free(out);
    free(ctx);
    return;
  }

  produced = 1;
  rc = inflate(comp, compSize, out, n, &produced, ctx);
  if (rc != INFLATE_OK || produced != n ||
      (n > 0 && memcmp(in, out, n) != 0))
    fail(name, "byte inflate mismatch");

  rc = inflate(comp, compSize, out, n, NULL, ctx);
  if (rc != INFLATE_OK || (n > 0 && memcmp(in, out, n) != 0))
    fail(name, "inflate with no size output failed");

  produced = 1;
  rc = inflate(comp, compSize, out, n, &produced, ctx);
  if (rc != INFLATE_OK || produced != n)
    fail(name, "reused context failed");

  vram = (unsigned short *)malloc((n / 2 + 2) * sizeof(unsigned short));
  if (vram == NULL) {
    fail(name, "out of memory");
  }
  else {
    produced = 1;
    rc = inflate_vram16(comp, compSize, vram, n, &produced, ctx);
    if (rc != INFLATE_OK || produced != n)
      fail(name, "vram inflate failed");
    else {
      for (i = 0; i < n; i++) {
        unsigned short word;
        unsigned char byte;

        word = vram[i >> 1];
        byte = (i & 1U) ? (unsigned char)(word >> 8) : (unsigned char)(word & 0xFF);
        if (byte != in[i]) {
          fail(name, "vram byte mismatch");
          break;
        }
      }
      if ((n & 1) != 0 && (vram[n >> 1] & 0xFF00) != 0)
        fail(name, "odd length left a stale high byte in VRAM");
    }
  }

  if (n > 0) {
    rc = inflate(comp, compSize, out, n - 1, &produced, ctx);
    if (rc != INFLATE_OUTPUT_OVERFLOW)
      fail(name, "small buffer was accepted");
  }

  free(vram);
  free(out);
  free(ctx);
  free(comp);
  (void)matches;
}


static void check_errors(void) {

  unsigned char *comp;
  unsigned char *out;
  unsigned char raw[8];
  struct InflateContext *ctx;
  int compSize, rc;
  int produced;
  unsigned char *misaligned;

  comp = NULL;
  compSize = 0;
  raw[0] = 1;
  raw[1] = 2;
  raw[2] = 3;
  raw[3] = 4;
  raw[4] = 5;
  raw[5] = 6;
  raw[6] = 7;
  raw[7] = 8;

  if (deflate_bytes(raw, 8, &comp, &compSize, NULL, NULL) != 0) {
    fail("errors", "could not compress the fixture");
    return;
  }

  ctx = (struct InflateContext *)malloc(sizeof(struct InflateContext));
  out = (unsigned char *)malloc(8);
  if (ctx == NULL || out == NULL) {
    fail("errors", "out of memory");
    free(comp);
    free(ctx);
    free(out);
    return;
  }

  rc = inflate(NULL, 8, out, 8, &produced, ctx);
  if (rc != INFLATE_BAD_ARGUMENT)
    fail("null data", "expected bad argument");

  rc = inflate(comp, compSize, out, 8, &produced, NULL);
  if (rc != INFLATE_BAD_ARGUMENT)
    fail("null context", "expected bad argument");

  rc = inflate(comp, 0, out, 8, &produced, ctx);
  if (rc != INFLATE_TRUNCATED)
    fail("empty input", "expected truncated");

  rc = inflate((const unsigned char *)"DEFd", 4, out, 8, &produced, ctx);
  if (rc != INFLATE_TRUNCATED)
    fail("four byte header", "expected truncated");

  rc = inflate((const unsigned char *)"XXXX\0\0\0\0", 8, out, 8, &produced, ctx);
  if (rc != INFLATE_WRONG_HEADER)
    fail("bad magic", "expected wrong header");

  if (compSize > 16) {
    rc = inflate(comp, compSize / 2, out, 8, &produced, ctx);
    if (rc == INFLATE_OK)
      fail("truncated body", "accepted a short buffer");
  }
  else
    fail("truncated body", "fixture was unexpectedly tiny");

  comp[4] = 9;
  comp[5] = 0;
  comp[6] = 0;
  comp[7] = 0;
  rc = inflate(comp, compSize, out, 9, &produced, ctx);
  if (rc != INFLATE_SIZE_MISMATCH)
    fail("size plus one", "expected a size mismatch");

  comp[4] = 0xFF;
  comp[5] = 0xFF;
  comp[6] = 0xFF;
  comp[7] = 0x7F;
  rc = inflate(comp, compSize, out, 8, &produced, ctx);
  if (rc != INFLATE_OUTPUT_OVERFLOW)
    fail("huge size", "expected output overflow");

  misaligned = (unsigned char *)malloc(8);
  if (misaligned != NULL) {
    rc = inflate_vram16(comp, compSize,
                        (volatile unsigned short *)(misaligned + 1),
                        8, &produced, ctx);
    if (rc != INFLATE_BAD_ARGUMENT)
      fail("unaligned vram", "expected bad argument");
    free(misaligned);
  }

  free(out);
  free(ctx);
  free(comp);
}


static void expect_matches(const char *name, const unsigned char *in, int n,
                           int wantMatches, int wantDup) {

  unsigned char *comp;
  char why[80];
  int compSize, matches, dups, rc;

  comp = NULL;
  compSize = 0;
  matches = -1;
  dups = -1;
  rc = deflate_bytes(in, n, &comp, &compSize, &matches, &dups);
  if (rc != 0 || matches != wantMatches || dups != wantDup) {
    sprintf(why, "got %d matches and %d duplicate bytes", matches, dups);
    fail(name, why);
  }
  free(comp);
}


static void check_deflate_args(void) {

  unsigned char raw[4];
  unsigned char *comp, *again;
  int compSize, againSize, matches, dups, rc;

  raw[0] = 'a';
  raw[1] = 'b';
  raw[2] = 'c';
  raw[3] = 'd';
  comp = NULL;
  again = NULL;
  compSize = 0;
  againSize = 0;

  if (deflate_bytes(raw, 4, NULL, &compSize, NULL, NULL) == 0)
    fail("deflate null output", "accepted a NULL output pointer");
  if (deflate_bytes(raw, 4, &comp, NULL, NULL, NULL) == 0)
    fail("deflate null size", "accepted a NULL size pointer");
  if (deflate_bytes(raw, -1, &comp, &compSize, NULL, NULL) == 0)
    fail("deflate negative", "accepted a negative size");
  if (deflate_bytes(NULL, 4, &comp, &compSize, NULL, NULL) == 0)
    fail("deflate null input", "accepted a NULL input");

  matches = -1;
  dups = -1;
  rc = deflate_bytes(NULL, 0, &comp, &compSize, &matches, &dups);
  if (rc != 0 || comp == NULL || matches != 0 || dups != 0 || compSize < 9)
    fail("deflate empty", "empty input failed");
  else if (comp[4] != 0 || comp[5] != 0 || comp[6] != 0 || comp[7] != 0)
    fail("deflate empty", "empty file stored a non-zero size");

  free(comp);
  comp = NULL;
  rc = deflate_bytes(raw, 4, &comp, &compSize, NULL, NULL);
  if (rc != 0 || comp == NULL)
    fail("deflate repeat", "compression failed");
  else {
    rc = deflate_bytes(raw, 4, &again, &againSize, NULL, NULL);
    if (rc != 0 || again == NULL || againSize != compSize ||
        memcmp(comp, again, againSize) != 0)
      fail("deflate repeat", "repeated compression was not identical");
  }

  free(comp);
  free(again);
}


static void fill_unique3(unsigned char *buf, int n) {

  int i, k;

  k = 0;
  for (i = 0; i + 2 < n; i += 3) {
    buf[i] = (unsigned char)(k & 255);
    buf[i + 1] = (unsigned char)((k >> 8) & 255);
    buf[i + 2] = 255;
    k++;
  }
  while (i < n) {
    buf[i] = 0;
    i++;
  }
}


static void check_window(unsigned char *buf) {

  unsigned char *comp, *out;
  struct InflateContext *ctx;
  char why[80];
  int compSize, matches, dups, produced, rc;

  comp = NULL;
  out = NULL;
  ctx = (struct InflateContext *)malloc(sizeof(struct InflateContext));
  if (ctx == NULL) {
    fail("window", "out of memory");
    return;
  }

  fill_unique3(buf, 32772);
  buf[32768] = buf[0];
  buf[32769] = buf[1];
  buf[32770] = buf[2];
  matches = -1;
  dups = -1;
  rc = deflate_bytes(buf, 32772, &comp, &compSize, &matches, &dups);
  if (rc != 0 || matches != 0 || dups != 0) {
    sprintf(why, "got %d matches and %d duplicate bytes", matches, dups);
    fail("far length 3", why);
  }
  out = (unsigned char *)malloc(32772);
  produced = 0;
  if (out != NULL && comp != NULL) {
    rc = inflate(comp, compSize, out, 32772, &produced, ctx);
    if (rc != INFLATE_OK || produced != 32772 || memcmp(buf, out, 32772) != 0)
      fail("far length 3", "round-trip mismatch");
  }
  free(comp);
  comp = NULL;
  free(out);

  fill_unique3(buf, 33027);
  memcpy(buf + 32768, buf, 258);
  matches = -1;
  dups = -1;
  rc = deflate_bytes(buf, 33027, &comp, &compSize, &matches, &dups);
  if (rc != 0 || matches != 1 || dups != 258) {
    sprintf(why, "got %d matches and %d duplicate bytes", matches, dups);
    fail("length 258 at distance 32768", why);
  }
  out = (unsigned char *)malloc(33027);
  produced = 0;
  if (out != NULL && comp != NULL) {
    rc = inflate(comp, compSize, out, 33027, &produced, ctx);
    if (rc != INFLATE_OK || produced != 33027 || memcmp(buf, out, 33027) != 0)
      fail("length 258 at distance 32768", "round-trip mismatch");
  }
  free(comp);
  comp = NULL;
  free(out);

  fill_unique3(buf, 32772);
  buf[32769] = buf[0];
  buf[32770] = buf[1];
  buf[32771] = buf[2];
  matches = -1;
  dups = -1;
  rc = deflate_bytes(buf, 32772, &comp, &compSize, &matches, &dups);
  if (rc != 0 || matches != 0 || dups != 0) {
    sprintf(why, "got %d matches and %d duplicate bytes", matches, dups);
    fail("distance 32769", why);
  }

  free(comp);
  free(ctx);
}


static void pack_init(unsigned char *bytes, int cap, int *nbits) {

  int i;

  for (i = 0; i < cap; i++)
    bytes[i] = 0;
  *nbits = 0;
}


static void pack_bits(unsigned char *bytes, int *nbits, unsigned int value, int n) {

  int i;

  for (i = n - 1; i >= 0; i--) {
    if (((value >> i) & 1U) != 0)
      bytes[*nbits >> 3] = (unsigned char)(bytes[*nbits >> 3] | (1U << (7 - (*nbits & 7))));
    (*nbits)++;
  }
}


static int build_bad_distance(unsigned char *data) {

  unsigned char bits[16];
  int nbits, nbytes, i;

  pack_init(bits, 16, &nbits);
  pack_bits(bits, &nbits, 1, 3);
  pack_bits(bits, &nbits, 0, 1);
  pack_bits(bits, &nbits, 1, 1);
  pack_bits(bits, &nbits, 0, 1);
  pack_bits(bits, &nbits, 0, 1);
  pack_bits(bits, &nbits, 1, 1);
  pack_bits(bits, &nbits, 1, 1);
  pack_bits(bits, &nbits, 127, 7);
  pack_bits(bits, &nbits, 1, 1);
  pack_bits(bits, &nbits, 107, 7);
  pack_bits(bits, &nbits, 0, 1);
  pack_bits(bits, &nbits, 0, 1);
  pack_bits(bits, &nbits, 1, 1);
  pack_bits(bits, &nbits, 17, 7);
  pack_bits(bits, &nbits, 0, 1);
  pack_bits(bits, &nbits, 1, 1);
  pack_bits(bits, &nbits, 18, 7);
  pack_bits(bits, &nbits, 1, 1);
  pack_bits(bits, &nbits, 0, 1);

  data[0] = 'D';
  data[1] = 'E';
  data[2] = 'F';
  data[3] = 'd';
  data[4] = 0;
  data[5] = 0;
  data[6] = 0;
  data[7] = 0;
  data[8] = 5;
  nbytes = (nbits + 7) / 8;
  for (i = 0; i < nbytes; i++)
    data[9 + i] = bits[i];
  return 9 + nbytes;
}


static void check_decoder_edges(void) {

  unsigned char *comp, *other, *out;
  unsigned char raw[8];
  unsigned char crafted[32];
  unsigned short vram[4];
  struct InflateContext *ctx;
  int compSize, otherSize, produced, rc, n;
  const char *text;

  comp = NULL;
  other = NULL;
  out = (unsigned char *)malloc(16);
  ctx = (struct InflateContext *)malloc(sizeof(struct InflateContext));
  if (out == NULL || ctx == NULL) {
    fail("decoder edges", "out of memory");
    free(out);
    free(ctx);
    return;
  }

  raw[0] = 'A';
  raw[1] = 'B';
  raw[2] = 'C';
  raw[3] = 'D';
  raw[4] = 'E';
  raw[5] = 'F';
  raw[6] = 'G';
  raw[7] = 'H';
  if (deflate_bytes(raw, 4, &comp, &compSize, NULL, NULL) != 0 ||
      deflate_bytes(raw, 0, &other, &otherSize, NULL, NULL) != 0) {
    fail("decoder edges", "could not compress a fixture");
    free(comp);
    free(other);
    free(out);
    free(ctx);
    return;
  }

  produced = 7;
  rc = inflate(other, otherSize, NULL, 0, &produced, ctx);
  if (rc != INFLATE_OK || produced != 0)
    fail("empty null output", "expected an empty success");

  rc = inflate(other, otherSize, NULL, 1, &produced, ctx);
  if (rc != INFLATE_BAD_ARGUMENT)
    fail("null output", "expected bad argument");

  rc = inflate(comp, -1, out, 4, &produced, ctx);
  if (rc != INFLATE_BAD_ARGUMENT)
    fail("negative input", "expected bad argument");

  rc = inflate(comp, compSize, out, -1, &produced, ctx);
  if (rc != INFLATE_BAD_ARGUMENT)
    fail("negative capacity", "expected bad argument");

  rc = inflate_vram16(comp, compSize, NULL, 4, &produced, ctx);
  if (rc != INFLATE_BAD_ARGUMENT)
    fail("null vram", "expected bad argument");

  rc = inflate_vram16(comp, -5, vram, 4, &produced, ctx);
  if (rc != INFLATE_BAD_ARGUMENT)
    fail("negative vram input", "expected bad argument");

  produced = 7;
  rc = inflate_vram16(other, otherSize, vram, 0, NULL, ctx);
  if (rc != INFLATE_OK)
    fail("empty vram", "expected success");

  rc = inflate(comp, 8, out, 4, &produced, ctx);
  if (rc != INFLATE_TRUNCATED)
    fail("eight byte header", "expected truncated input");

  comp[4] = 0;
  comp[5] = 0;
  comp[6] = 0;
  comp[7] = 0;
  rc = inflate(comp, compSize, out, 4, &produced, ctx);
  if (rc != INFLATE_SIZE_MISMATCH)
    fail("claimed zero", "expected a size mismatch");

  free(comp);
  comp = NULL;
  if (deflate_bytes(raw, 4, &comp, &compSize, NULL, NULL) != 0) {
    fail("decoder edges", "could not recompress the fixture");
    free(comp);
    free(other);
    free(out);
    free(ctx);
    return;
  }
  comp[4] = 3;
  comp[5] = 0;
  comp[6] = 0;
  comp[7] = 0;
  rc = inflate(comp, compSize, out, 4, &produced, ctx);
  if (rc != INFLATE_SIZE_MISMATCH)
    fail("size minus one", "expected a size mismatch");

  comp[4] = 0xFF;
  comp[5] = 0xFF;
  comp[6] = 0xFF;
  comp[7] = 0xFF;
  rc = inflate(comp, compSize, out, 4, &produced, ctx);
  if (rc != INFLATE_OUTPUT_OVERFLOW)
    fail("size above int", "expected output overflow");

  comp[4] = 4;
  comp[5] = 0;
  comp[6] = 0;
  comp[7] = 0;
  comp[8] = 0;
  rc = inflate(comp, 9, out, 4, &produced, ctx);
  if (rc != INFLATE_TRUNCATED)
    fail("stored header", "expected truncated input");

  comp[8] = 3;
  rc = inflate(comp, compSize, out, 4, &produced, ctx);
  if (rc != INFLATE_BAD_CODE)
    fail("codesN 3", "expected a bad code");

  comp[8] = 120;
  rc = inflate(comp, compSize, out, 4, &produced, ctx);
  if (rc != INFLATE_BAD_CODE)
    fail("codesN 120", "expected a bad code");

  crafted[0] = 'D';
  crafted[1] = 'E';
  crafted[2] = 'F';
  crafted[3] = 'd';
  crafted[4] = 0;
  crafted[5] = 0;
  crafted[6] = 0;
  crafted[7] = 0;
  crafted[8] = 4;
  crafted[9] = 0;
  rc = inflate(crafted, 10, out, 0, &produced, ctx);
  if (rc != INFLATE_BAD_CODE)
    fail("width 0", "expected a bad code");

  pack_init(crafted + 9, 8, &n);
  pack_bits(crafted + 9, &n, 1, 3);
  pack_bits(crafted + 9, &n, 1, 1);
  pack_bits(crafted + 9, &n, 1, 1);
  pack_bits(crafted + 9, &n, 1, 1);
  pack_bits(crafted + 9, &n, 1, 1);
  rc = inflate(crafted, 9 + (n + 7) / 8, out, 0, &produced, ctx);
  if (rc != INFLATE_BAD_CODE)
    fail("oversubscribed", "expected a bad code");

  pack_init(crafted + 9, 8, &n);
  pack_bits(crafted + 9, &n, 6, 3);
  pack_bits(crafted + 9, &n, 32, 6);
  rc = inflate(crafted, 9 + (n + 7) / 8, out, 0, &produced, ctx);
  if (rc != INFLATE_BAD_CODE)
    fail("length 32", "expected a bad code");

  n = build_bad_distance(crafted);
  rc = inflate(crafted, n, NULL, 0, &produced, ctx);
  if (rc != INFLATE_BAD_DISTANCE)
    fail("distance before output", "expected a bad distance");

  free(comp);
  free(other);
  comp = NULL;
  other = NULL;
  if (deflate_bytes((const unsigned char *)"Hello", 5, &comp, &compSize, NULL, NULL) != 0 ||
      deflate_bytes((const unsigned char *)"AAAA", 4, &other, &otherSize, NULL, NULL) != 0) {
    fail("reused context", "could not compress the fixtures");
  }
  else {
    produced = 0;
    rc = inflate(comp, compSize, out, 5, &produced, ctx);
    if (rc != INFLATE_OK || produced != 5 || memcmp(out, "Hello", 5) != 0)
      fail("reused context", "first file failed");
    produced = 0;
    rc = inflate(other, otherSize, out, 4, &produced, ctx);
    if (rc != INFLATE_OK || produced != 4 || memcmp(out, "AAAA", 4) != 0)
      fail("reused context", "second file failed");
    produced = 0;
    rc = inflate_vram16(other, otherSize, vram, 4, NULL, ctx);
    if (rc != INFLATE_OK)
      fail("reused context", "vram pass failed");
    else if ((vram[0] & 0xFF) != 'A' || (vram[0] >> 8) != 'A' ||
             (vram[1] & 0xFF) != 'A' || (vram[1] >> 8) != 'A')
      fail("reused context", "vram bytes mismatch");
  }

  text = inflate_strerror(INFLATE_OK);
  if (text == NULL || strcmp(text, "ok") != 0)
    fail("strerror", "missing ok");
  if (inflate_strerror(INFLATE_WRONG_HEADER) == NULL ||
      inflate_strerror(INFLATE_TRUNCATED) == NULL ||
      inflate_strerror(INFLATE_BAD_CODE) == NULL ||
      inflate_strerror(INFLATE_OUTPUT_OVERFLOW) == NULL ||
      inflate_strerror(INFLATE_BAD_DISTANCE) == NULL ||
      inflate_strerror(INFLATE_SIZE_MISMATCH) == NULL ||
      inflate_strerror(INFLATE_BAD_ARGUMENT) == NULL)
    fail("strerror", "missing a known error");
  text = inflate_strerror(99);
  if (text == NULL || strcmp(text, "unknown inflate error") != 0)
    fail("strerror", "unknown code was not reported");

  free(comp);
  free(other);
  free(out);
  free(ctx);
}


/* One epoch is 4 pattern bytes, 256 four-byte traps, 2049 zeros and
   256 five-byte blocks. The length-3 trap for each lead byte sits
   3073..3328 bytes back: 13 extra distance bits, still under the
   matcher's 4096 cutoff, so this is the lazy decision and not the
   far length-3 rejection. The next byte is a length-4 copy of the
   pattern. A greedy matcher keeps the length-3 trap and stores the
   last two pattern bytes as literals (one epoch 1061 bytes and 3584
   duplicate bytes; sixteen epochs 20022 bytes and 57419 duplicate
   bytes). Lazy emits the lead byte and keeps the length-4 match. */
static int fill_lazy_traps(unsigned char *buf, int epochs) {

  int e, m, n, base;
  unsigned char pat[4];

  n = 0;
  for (e = 0; e < epochs; e++) {
    base = 1 + (e * 4) % 200;
    pat[0] = (unsigned char)(200 + ((base + 0) % 50));
    pat[1] = (unsigned char)(200 + ((base + 1) % 50));
    pat[2] = (unsigned char)(200 + ((base + 2) % 50));
    pat[3] = (unsigned char)(200 + ((base + 3) % 50));
    memcpy(buf + n, pat, 4);
    n += 4;
    for (m = 0; m < 256; m++) {
      buf[n++] = (unsigned char)m;
      buf[n++] = pat[0];
      buf[n++] = pat[1];
      buf[n++] = 0;
    }
    memset(buf + n, 0, 2049);
    n += 2049;
    for (m = 0; m < 256; m++) {
      buf[n++] = (unsigned char)m;
      memcpy(buf + n, pat, 4);
      n += 4;
    }
  }
  return n;
}


static void check_lazy_gain(void) {

  unsigned char *buf;
  unsigned char *comp;
  char why[80];
  int n, compSize, matches, dups, rc;

  buf = (unsigned char *)malloc(4357);
  if (buf == NULL) {
    fail("lazy distant traps", "out of memory");
    return;
  }
  n = fill_lazy_traps(buf, 1);
  comp = NULL;
  compSize = 0;
  matches = -1;
  dups = -1;
  rc = deflate_bytes(buf, n, &comp, &compSize, &matches, &dups);
  if (n != 4357 || rc != 0 || matches != 519 || dups != 3839) {
    sprintf(why, "got %d matches and %d duplicate bytes", matches, dups);
    fail("lazy distant traps", why);
  }
  else if (compSize > 980) {
    sprintf(why, "compressed to %d bytes", compSize);
    fail("lazy distant traps", why);
  }
  free(comp);
  check_roundtrip("lazy distant traps", buf, n);
  free(buf);

  buf = (unsigned char *)malloc(16 * 4357);
  if (buf == NULL) {
    fail("lazy distant traps x16", "out of memory");
    return;
  }
  n = fill_lazy_traps(buf, 16);
  comp = NULL;
  compSize = 0;
  matches = -1;
  dups = -1;
  rc = deflate_bytes(buf, n, &comp, &compSize, &matches, &dups);
  /* 16 epochs keep 8304 matches. Cross-epoch copies that add 30
     duplicate bytes do not pay, so the count is 16 * 3839. */
  if (n != 16 * 4357 || rc != 0 || matches != 8304 || dups != 61424) {
    sprintf(why, "got %d matches and %d duplicate bytes", matches, dups);
    fail("lazy distant traps x16", why);
  }
  else if (compSize > 15000) {
    sprintf(why, "compressed to %d bytes", compSize);
    fail("lazy distant traps x16", why);
  }
  free(comp);
  check_roundtrip("lazy distant traps x16", buf, n);
  free(buf);
}


/* Each lead byte has a length-20 match, and the next byte starts a
   length-80 copy. Length 20 is past the old 16-byte lookahead cutoff,
   which stores this file as 1382 bytes and 606 matches. Looking ahead
   for every match stores it as 1128 bytes and 406 matches. */
static int fill_long_lazy(unsigned char *buf) {

  int i, m, n;
  unsigned char pat[80];
  unsigned char stop;

  for (i = 0; i < 80; i++)
    pat[i] = (unsigned char)((i * 17 + 3) & 255);
  stop = (unsigned char)((pat[19] != 0xFF) ? 0xFF : 0xFE);
  n = 0;
  memcpy(buf + n, pat, 80);
  n += 80;
  for (m = 0; m < 200; m++) {
    buf[n++] = (unsigned char)(m & 255);
    memcpy(buf + n, pat, 19);
    n += 19;
    buf[n++] = stop;
  }
  memset(buf + n, 0x11, 1500);
  n += 1500;
  for (m = 0; m < 200; m++) {
    buf[n++] = (unsigned char)(m & 255);
    memcpy(buf + n, pat, 80);
    n += 80;
  }
  return n;
}


static void check_long_lazy(void) {

  unsigned char *buf;
  unsigned char *comp;
  char why[80];
  int n, compSize, matches, dups, rc;

  buf = (unsigned char *)malloc(21980);
  if (buf == NULL) {
    fail("lazy long match", "out of memory");
    return;
  }
  n = fill_long_lazy(buf);
  comp = NULL;
  compSize = 0;
  matches = -1;
  dups = -1;
  rc = deflate_bytes(buf, n, &comp, &compSize, &matches, &dups);
  if (n != 21980 || rc != 0 || matches != 406 || dups != 21499) {
    sprintf(why, "got %d matches and %d duplicate bytes", matches, dups);
    fail("lazy long match", why);
  }
  else if (compSize > 1200) {
    sprintf(why, "compressed to %d bytes", compSize);
    fail("lazy long match", why);
  }
  free(comp);
  check_roundtrip("lazy long match", buf, n);
  free(buf);
}


static void fill_lcg(unsigned char *buf, int n, unsigned long seed);


static void put_u32le(unsigned char *p, int n) {

  p[0] = (unsigned char)(n & 0xFF);
  p[1] = (unsigned char)((n >> 8) & 0xFF);
  p[2] = (unsigned char)((n >> 16) & 0xFF);
  p[3] = (unsigned char)((n >> 24) & 0xFF);
}


static void check_blocks(void) {

  unsigned char *buf;
  unsigned char *comp;
  unsigned char raw[24];
  unsigned char out[8];
  struct InflateContext *ctx;
  char why[80];
  int i, compSize, rc, produced;

  buf = (unsigned char *)malloc(16000);
  ctx = (struct InflateContext *)malloc(sizeof(struct InflateContext));
  if (buf == NULL || ctx == NULL) {
    fail("blocks", "out of memory");
    free(buf);
    free(ctx);
    return;
  }

  fill_lcg(buf, 2000, 7UL);
  comp = NULL;
  compSize = 0;
  rc = deflate_bytes(buf, 2000, &comp, &compSize, NULL, NULL);
  if (rc != 0 || comp == NULL || compSize > 2013) {
    sprintf(why, "compressed to %d bytes", compSize);
    fail("stored lcg", why);
  }
  free(comp);
  check_roundtrip("stored lcg", buf, 2000);

  for (i = 0; i < 8000; i++)
    buf[i] = (unsigned char)(i % 20);
  for (i = 8000; i < 16000; i++)
    buf[i] = (unsigned char)(200 + (i % 20));
  memcpy(buf + 8000, buf + 7800, 120);
  comp = NULL;
  compSize = 0;
  rc = deflate_bytes(buf, 16000, &comp, &compSize, NULL, NULL);
  if (rc != 0 || comp == NULL || compSize > 10000) {
    sprintf(why, "compressed to %d bytes", compSize);
    fail("split alphabets", why);
  }
  free(comp);
  check_roundtrip("split alphabets", buf, 16000);

  memset(raw, 0, 24);
  raw[0] = 'D';
  raw[1] = 'E';
  raw[2] = 'F';
  raw[3] = 'd';
  raw[8] = 0;
  put_u32le(raw + 4, 2);
  put_u32le(raw + 9, 2);
  raw[13] = 'A';
  raw[14] = 'B';
  produced = 0;
  rc = inflate(raw, 15, out, 2, &produced, ctx);
  if (rc != INFLATE_OK || produced != 2 || out[0] != 'A' || out[1] != 'B')
    fail("stored block", "did not return the stored bytes");

  put_u32le(raw + 4, 4);
  put_u32le(raw + 9, 0);
  rc = inflate(raw, 13, out, 4, &produced, ctx);
  if (rc != INFLATE_BAD_CODE)
    fail("stored length 0", "expected a bad code");

  put_u32le(raw + 9, 5);
  rc = inflate(raw, 13, out, 4, &produced, ctx);
  if (rc != INFLATE_SIZE_MISMATCH)
    fail("stored length past end", "expected a size mismatch");

  put_u32le(raw + 9, 2);
  raw[13] = 'A';
  raw[14] = 'B';
  rc = inflate(raw, 15, out, 4, &produced, ctx);
  if (rc != INFLATE_SIZE_MISMATCH)
    fail("stored then eof", "expected a size mismatch");

  raw[15] = 2;
  rc = inflate(raw, 16, out, 4, &produced, ctx);
  if (rc != INFLATE_BAD_CODE)
    fail("bad block type", "expected a bad code");

  free(buf);
  free(ctx);
}


static void fill_lcg(unsigned char *buf, int n, unsigned long seed) {

  int i;
  unsigned long s;

  s = seed;
  for (i = 0; i < n; i++) {
    s = s * 1103515245UL + 12345UL;
    buf[i] = (unsigned char)(s >> 16);
  }
}


int main(void) {

  unsigned char *buf;
  unsigned char hello[13];
  int i;

  g_fail = 0;
  memcpy(hello, "Hello, world!", 13);

  check_huffman();
  check_huffman_limits();
  check_deflate_args();
  check_roundtrip("empty", hello, 0);
  check_roundtrip("hello", hello, 13);
  check_roundtrip("A", (const unsigned char *)"A", 1);
  check_roundtrip("AA", (const unsigned char *)"AA", 2);
  check_roundtrip("AAA", (const unsigned char *)"AAA", 3);
  check_roundtrip("AAAA", (const unsigned char *)"AAAA", 4);
  check_roundtrip("abc", (const unsigned char *)"abc", 3);
  check_roundtrip("abcabc", (const unsigned char *)"abcabc", 6);
  check_roundtrip("aaaabaaa", (const unsigned char *)"aaaabaaa", 8);
  expect_matches("hello matches", hello, 13, 0, 0);
  expect_matches("AAA matches", (const unsigned char *)"AAA", 3, 0, 0);
  /* A length-3 match does not make AAAA smaller than a stored block.
     aaabaaaab is smaller as literals than as one match. */
  expect_matches("AAAA matches", (const unsigned char *)"AAAA", 4, 0, 0);
  expect_matches("abcabc matches", (const unsigned char *)"abcabc", 6, 1, 3);
  expect_matches("lazy aaabaaaab", (const unsigned char *)"aaabaaaab", 9, 0, 0);
  check_roundtrip("lazy aaabaaaab", (const unsigned char *)"aaabaaaab", 9);

  buf = (unsigned char *)malloc(40000);
  if (buf == NULL) {
    fail("buffers", "out of memory");
    return 1;
  }

  for (i = 0; i < 300; i++)
    buf[i] = 'Q';
  expect_matches("300 Q matches", buf, 300, 2, 299);
  check_roundtrip("300 Q", buf, 300);

  for (i = 0; i < 256; i++)
    buf[i] = (unsigned char)i;
  check_roundtrip("bytes 0..255", buf, 256);

  for (i = 0; i < 258; i++)
    buf[i] = (unsigned char)(i & 255);
  check_roundtrip("258 sequential", buf, 258);

  for (i = 0; i < 259; i++)
    buf[i] = 'Z';
  expect_matches("259 Z matches", buf, 259, 1, 258);
  check_roundtrip("259 Z", buf, 259);

  for (i = 0; i < 260; i++)
    buf[i] = 'Z';
  expect_matches("260 Z matches", buf, 260, 1, 258);
  check_roundtrip("260 Z", buf, 260);

  fill_lcg(buf, 4096, 1UL);
  check_roundtrip("lcg 4096", buf, 4096);

  fill_lcg(buf, 37568, 99UL);
  check_roundtrip("lcg 37568", buf, 37568);

  for (i = 0; i < 40000; i++)
    buf[i] = (unsigned char)(i & 255);
  memcpy(buf + 32768, buf, 100);
  check_roundtrip("distance 32768", buf, 40000);

  check_window(buf);
  check_lazy_gain();
  check_long_lazy();
  check_blocks();
  check_errors();
  check_decoder_edges();

  free(buf);

  if (g_fail != 0) {
    fprintf(stderr, "%d failure(s)\n", g_fail);
    return 1;
  }
  fprintf(stderr, "all tests passed\n");
  return 0;
}
