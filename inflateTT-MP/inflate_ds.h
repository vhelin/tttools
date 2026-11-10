#ifndef INFLATE_DS_H
#define INFLATE_DS_H

/* Nintendo DS entry point. Compile inflate_ds.c instead of inflate.c.
   Linking both objects defines inflate() twice.

   output must be 2-byte aligned. outputBytes is the destination capacity
   in bytes and must be at least the uncompressed size stored in the
   DEFc header. One static context is used, so this is not re-entrant.
   A 0 return is success. Any other value is an INFLATE_* code from
   inflate.h. Include that header only from a file that does not also
   include this one: both headers declare inflate(). */

int inflate(const unsigned char *data, int dataSize, volatile unsigned short *output, int outputBytes);

#endif
