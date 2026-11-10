/*
 * inflateTT, decompresses data that has been compressed using deflateTT.
 *
 * Programmed by Ville Helin <vhelin#iki.fi> in 2007.
 * File I/O around inflateTT-MP, with size checks, added in 2026.
 *
 * This code is under GNU General Public Licence (GPL), version 2, June 1991.
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "inflate.h"

#define INFLATE_CLI_MAX_OUTPUT (256UL * 1024UL * 1024UL)

int main(int argc, char *argv[]) {

  unsigned char *data, *output;
  struct InflateContext *context;
  FILE *f;
  long fileSizeLong;
  int fileSize, claimed, produced;
  unsigned long claimedLong;
  int rc;

  data = NULL;
  output = NULL;
  context = NULL;
  f = NULL;

  if (argc != 3) {
    fprintf(stderr, "inflateTT v2.1 Written by Ville Helin 2007\n");
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
  f = NULL;

  if (fileSize < 4 || data[0] != 'D' || data[1] != 'E' || data[2] != 'F' || data[3] != 'd') {
    fprintf(stderr, "main(): File \"%s\" doesn't start with \"DEFd\".\n", argv[1]);
    free(data);
    return 1;
  }
  if (fileSize < 8) {
    fprintf(stderr, "main(): %s.\n", inflate_strerror(INFLATE_TRUNCATED));
    free(data);
    return 1;
  }

  claimedLong = (unsigned long)data[4]
              | ((unsigned long)data[5] << 8)
              | ((unsigned long)data[6] << 16)
              | ((unsigned long)data[7] << 24);
  if (claimedLong > INFLATE_CLI_MAX_OUTPUT) {
    fprintf(stderr, "main(): Uncompressed size %lu is unreasonably large.\n", claimedLong);
    free(data);
    return 1;
  }
  claimed = (int)claimedLong;

  output = (unsigned char *)malloc(claimed > 0 ? claimed : 1);
  context = (struct InflateContext *)malloc(sizeof(struct InflateContext));
  if (output == NULL || context == NULL) {
    fprintf(stderr, "main(): Out of memory error.\n");
    free(data);
    free(output);
    free(context);
    return 1;
  }

  produced = 0;
  rc = inflate(data, fileSize, output, claimed, &produced, context);
  if (rc != INFLATE_OK) {
    fprintf(stderr, "main(): %s.\n", inflate_strerror(rc));
    free(data);
    free(output);
    free(context);
    return 1;
  }

  fprintf(stderr, "main(): Original size = %d, uncompressed size = %d.\n",
          claimed, produced);

  f = fopen(argv[2], "wb");
  if (f == NULL) {
    fprintf(stderr, "main(): Could not open file \"%s\" for writing.\n", argv[2]);
    free(data);
    free(output);
    free(context);
    return 1;
  }
  if (produced > 0 && (int)fwrite(output, 1, produced, f) != produced) {
    fprintf(stderr, "main(): Could not write \"%s\".\n", argv[2]);
    fclose(f);
    free(data);
    free(output);
    free(context);
    return 1;
  }
  fclose(f);

  free(data);
  free(output);
  free(context);
  return 0;
}
