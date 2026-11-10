#ifndef DEFLATE_H
#define DEFLATE_H

/* Compress input into a newly malloc'd DEFc buffer. *output is set on success and must be freed by the caller.
   matches and duplicateBytes may be NULL. Returns 0 on success and -1 on failure. */
int deflate_bytes(const unsigned char *input, int inputSize, unsigned char **output, int *outputSize, int *matches, int *duplicateBytes);

#endif
