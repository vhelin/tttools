
#ifndef INFLATE_H
#define INFLATE_H

#ifdef __cplusplus
extern "C" {
#endif

/* the number of bits we can have in a code */
#define HUFFMAN_CODE_MAX_BITS 32

struct InflateNode {
	struct InflateNode *left;
	struct InflateNode *right;
	int literal;
};

/* the inflate context */
struct InflateContext {
	/* tmp buffers for huffmanRecreateCodes() */
	int tmpCount[HUFFMAN_CODE_MAX_BITS];
	int nextCode[HUFFMAN_CODE_MAX_BITS+1];

	/* code lengths */
	int codeLengthLiterals[286];
	int codeLengthDistances[30];
	int codeLengthCombined[119];

	/* codes */
	int codeLiterals[286];
	int codeDistances[30];
	int codeCombined[119];

	/* the huffman trees */
	struct InflateNode *treeLiterals;
	struct InflateNode *treeDistances;
	struct InflateNode *treeCombined;

	/* the free tree nodes */
	int treeNodesFreeCurrent;
	struct InflateNode treeNodesFree[1024];
};

/* the return values */
#define INFLATE_OK           0
#define INFLATE_WRONG_HEADER 1

int inflate(unsigned char *data, unsigned char *output, struct InflateContext *context);

#ifdef __cplusplus
}
#endif

#endif
