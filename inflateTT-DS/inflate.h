
#ifndef INFLATE_H
#define INFLATE_H

struct node {
	struct node *left;
	struct node *right;
	s32 literal;
};

/* inflate() decompresses the deflateTT'ed data to output, with 16bit writes */
void inflate(u8 *data, vu16 *output);

#endif
