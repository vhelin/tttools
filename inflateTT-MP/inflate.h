#ifndef INFLATE_H
#define INFLATE_H

#ifdef __cplusplus
extern "C" {
#endif

/* The compressor emits codes of at most 15 bits. The decoder accepts
   codes up to 31 bits so files written by older deflateTT builds,
   which did not limit code length, still decompress. */
#define INFLATE_MAX_BITS 31

#define INFLATE_LIT_SYMBOLS 286
#define INFLATE_DIST_SYMBOLS 30
#define INFLATE_CODELEN_SYMBOLS 119
#define INFLATE_MAX_TREE_NODES 1024

struct InflateNode {
  struct InflateNode *left;
  struct InflateNode *right;
  int literal;
};

/* About 17 KB on a 32-bit target and about 30 KB on a 64-bit target,
   because the node pool is embedded. Allocate it with malloc or keep
   it in static storage. Do not put it on a small stack. */
struct InflateContext {
  int tmpCount[INFLATE_MAX_BITS + 1];
  int nextCode[INFLATE_MAX_BITS + 1];

  int codeLengthLiterals[INFLATE_LIT_SYMBOLS];
  int codeLengthDistances[INFLATE_DIST_SYMBOLS];
  int codeLengthCombined[INFLATE_CODELEN_SYMBOLS];

  int codeLiterals[INFLATE_LIT_SYMBOLS];
  int codeDistances[INFLATE_DIST_SYMBOLS];
  int codeCombined[INFLATE_CODELEN_SYMBOLS];

  struct InflateNode *treeLiterals;
  struct InflateNode *treeDistances;
  struct InflateNode *treeCombined;

  int treeNodesFreeCurrent;
  struct InflateNode treeNodesFree[INFLATE_MAX_TREE_NODES];
};

#define INFLATE_OK              0
#define INFLATE_WRONG_HEADER    1
#define INFLATE_TRUNCATED       2
#define INFLATE_BAD_CODE        3
#define INFLATE_OUTPUT_OVERFLOW 4
#define INFLATE_BAD_DISTANCE    5
#define INFLATE_SIZE_MISMATCH   6
#define INFLATE_BAD_ARGUMENT    7

/* Decompress one DEFd buffer.
   outputCapacity must be at least the uncompressed size stored in the
   header. On success *outputSize (when non-NULL) is that many bytes.
   The output buffer may be partially written when the return is not
   INFLATE_OK. A NULL output is accepted only when outputCapacity is 0
   and the uncompressed size is 0. */
int inflate(const unsigned char *data, int dataSize,
            unsigned char *output, int outputCapacity,
            int *outputSize, struct InflateContext *context);

/* Same decoder, storing bytes as little-endian 16-bit writes.
   output must be 2-byte aligned. outputBytes is the capacity in bytes.
   Lengths are int, the same type as the rest of the byte counts.
   inflate_ds.c wraps this with one static context and exports it as
   inflate(). Compile that file instead of inflate.c on the DS. */
int inflate_vram16(const unsigned char *data, int dataSize,
                   volatile unsigned short *output, int outputBytes,
                   int *outputSize, struct InflateContext *context);

const char *inflate_strerror(int error);

#ifdef __cplusplus
}
#endif

#endif
