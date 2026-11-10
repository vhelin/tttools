/*
 * Nintendo DS entry point for the inflateTT-MP decoder.
 * 16-bit stores, so the destination can be VRAM.
 *
 * Programmed by Ville Helin <vhelin#iki.fi> in 2007.
 *
 * This code is under GNU Lesser General Public Licence (LGPL), version 2.1,
 * February 1999.
 *
 * Compile this file instead of inflate.c. The included decoder's
 * byte-oriented inflate() is renamed, and this file exports inflate()
 * as the 16-bit writer with one static context.
 */

#define inflate inflate_bytes
#include "inflate.c"
#undef inflate

#include "inflate_ds.h"

static struct InflateContext dsContext;

int inflate(const unsigned char *data, int dataSize, volatile unsigned short *output, int outputBytes) {

  return inflate_vram16(data, dataSize, output, outputBytes, NULL, &dsContext);
}
