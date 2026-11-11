inflateTT (GPL v2) and deflateTT (GPL v2) are command-line tools that compress
and decompress data. They were used in Explosion of an Aeroplane to compress
graphics. inflateTT-MP (LGPL v2.1) is the portable library. It can write
bytes, or 16-bit words so a Nintendo DS destination can be VRAM.

Literal bytes, length codes 257..285 and distance codes 0..29, including the
extra-bit tables, match RFC 1951. The container does not. A DEFd file is not
zlib, gzip, or a raw deflate block. The format is described below.

Build and test with gcc:

```
make
make test
```

Cygwin's `run` looks up `tests.exe` before the tests directory, so this also works:

```
run -wait tests
```

deflateTT reads a raw file and writes DEFd. inflateTT does the reverse.

```
deflateTT input.raw output.def
inflateTT output.def output.raw
```

The compressor may write several blocks. Each block is Huffman-coded
or stored uncompressed, whichever is smaller, and a match can refer to
a byte written by an earlier block. Matches are LZ77 over a 32768-byte
window, found with a 3-byte hash chain and chosen by bit cost, so a
shorter nearer copy is kept when it encodes in fewer bits. A split is
kept only when the blocks are smaller than one block.
Code lengths are limited to 15 bits. An alphabet with fewer than two used
symbols is given two length-1 codes, so a run of identical bytes and a file
with no repeated 3-byte sequence both round-trip.

## inflateTT-MP

```c
int inflate(const unsigned char *data, int dataSize,
            unsigned char *output, int outputCapacity,
            int *outputSize, struct InflateContext *context);
```

`outputCapacity` must be at least the uncompressed size stored in the header.
The return value is `INFLATE_OK` or another `INFLATE_*` code from `inflate.h`.
`InflateContext` holds the Huffman node pool and the 9-bit decode tables.
Measured size is 40,504 bytes when pointers are 64 bits. The same layout is
28,200 bytes when pointers are 32 bits. Allocate it with `malloc`, or keep
it in static storage.

On the DS, compile `inflateTT-MP/inflate_ds.c` instead of `inflate.c`. Do not
link both: each defines `inflate()`. The DS entry point writes 16-bit words
and keeps one static context, so it is not re-entrant:

```c
int inflate(const unsigned char *data, int dataSize,
            volatile unsigned short *output, int outputBytes);
```

`output` must be 2-byte aligned. `outputBytes` is the capacity in bytes. A
program that calls `inflate_vram16` and supplies its own context compiles
`inflate.c` instead, and does not use `inflate_ds.c`.

## DEFd format

Offset 0, 4 bytes: `D` `E` `F` `d`.
Offset 4, uint32 little-endian: uncompressed size.
The file is one or more blocks. Each block starts on a byte boundary.
Offset 8, and the first byte of every later block, is a type byte T.

T = 0 is a stored block. The next four bytes are a little-endian length
of at least 1, then that many raw bytes. The length must not run past
the uncompressed size in the header.

T in 4..119 is a Huffman block. T is the code-length alphabet size N.
The encoder sets N to (longest literal or distance code) + 4. The
following bits are packed MSB first (bit 7 of the next byte is first):

- 3 bits: width W of the following code lengths. The encoder uses 2..7.
- N values of W bits each: code lengths for the code-length alphabet.
- A Huffman stream of those symbols, expanding to exactly 316 lengths:
  286 literal/length lengths, then 30 distance lengths.
  Symbol <= N-4 is a literal length.
  Symbol N-3 repeats the previous length. Two extra bits, count = extra+3
  (3..6 extra copies, in addition to the length already emitted).
  Symbol N-2 is a run of zeros. Three extra bits, count = extra+3 (3..10).
  Symbol N-1 is a run of zeros. Seven extra bits, count = extra+11 (11..138).
- A Huffman stream of literal/length symbols, then distances.
  0..255 are bytes, 256 ends this block, 257..285 are match lengths.
  Extra length and distance bits match RFC 1951.
  A match distance must point at a byte already written, including a
  byte written by an earlier block in this file.

If the bytes written so far are still short of the uncompressed size,
the decoder moves to the next byte and reads another block. Padding
bits in the last byte of a Huffman block are ignored. When the count
matches the header, the file ends. A file with one Huffman block is
the same bitstream the previous DEFc container used after its magic.

Canonical codes are assigned from the lengths as in RFC 1951 section 3.2.2.
The encoder will not emit a code longer than 15 bits. The decoder rejects
a longer code.

## Compared with zlib

Measured on 10 Nov 2026 on an AMD Ryzen 9 7950X, gcc `-O2`, zlib 1.3.2.
Sizes are the whole container. Every output was inflated and compared
with the input. `deflate_bytes` times are the mean of 3 calls after one
warmup for files of 60 KiB or more, and 5 calls for smaller files. zlib
was called three times as often. Each timed call allocates its output
buffer. The best zlib size is the smallest of levels 1 through 9, which
is not always level 9.

Hello is the 13 bytes of `Hello, world!`. Zeros and the single-byte run
are 64 KiB. The sentence is `The quick brown fox jumps over the lazy dog. `
repeated, and the edited file XORs every 53rd byte with `0x5A`. Two
periods is 32 KiB of `i mod 20`, then 32 KiB of `200 + (i mod 20)`, with
120 bytes from the start copied onto offset 32704. Two alphabets is an
LCG restricted to 0..39 and then to 180..255. Text then noise is 32 KiB
of the sentence and 32 KiB of LCG noise. Period 100 repeats a 100-byte
pattern. The far copy plants 4096 bytes from the start at offset 60000,
past the 32768-byte window, so both DEFd and zlib store it. Lazy traps
are 16 epochs of the distant length-3 fixture. Long matches are the
21980-byte lookahead fixture. Sources is one snapshot of the READMEs plus
`deflateTT` main/huffman, `inflateTT-MP` inflate.c/h, `inflateTT` main,
and `tests/roundtrip.c` (125,833 bytes).

deflateTT is a little smaller than the best zlib level on the repetitive
files, and 4.6% smaller on that source snapshot (27,063 against 28,373).
Hello is the size loss, 26 bytes against 21. Incompressible input is
stored by both. DEFd's stored framing is 13 bytes and zlib's is 26.

| File | Raw | DEFd | zlib -1 | zlib -6 | zlib -9 | Best zlib |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| Hello | 13 | 26 | 21 | 21 | 21 | 21 (level 1) |
| Zeros | 65,536 | 82 | 307 | 84 | 84 | 84 (level 4) |
| One byte | 65,536 | 83 | 308 | 85 | 85 | 85 (level 4) |
| Sentence | 65,536 | 256 | 456 | 258 | 258 | 258 (level 4) |
| Edited sentence | 65,536 | 628 | 1,209 | 647 | 647 | 647 (level 6) |
| Two periods | 65,536 | 228 | 461 | 232 | 232 | 232 (level 6) |
| Two alphabets | 65,536 | 48,948 | 51,234 | 50,028 | 50,028 | 50,028 (level 4) |
| Text then noise | 65,536 | 32,934 | 33,194 | 33,089 | 33,089 | 33,089 (level 4) |
| Noise | 65,536 | 65,549 | 65,562 | 65,562 | 65,562 | 65,562 (level 1) |
| Period 100 | 65,536 | 391 | 631 | 391 | 391 | 391 (level 4) |
| Copy past 32 KiB | 65,536 | 65,549 | 65,562 | 65,562 | 65,562 | 65,562 (level 1) |
| Lazy traps | 69,712 | 13,654 | 13,783 | 13,767 | 13,767 | 13,751 (level 4) |
| Long matches | 21,980 | 1,073 | 1,144 | 1,376 | 1,133 | 1,133 (level 7) |
| Sources | 125,833 | 27,063 | 35,850 | 28,488 | 28,374 | 28,373 (level 8) |

Compression is much slower than zlib. A 64 KiB run of zeros takes 44.7 ms
in deflateTT and 0.127 ms at zlib level 9. The source snapshot takes 131 ms
(0.96 MB/s) against 5.85 ms at zlib level 9 (21.5 MB/s) and 0.69 ms at
level 1. A megabyte here is 1,000,000 bytes. Noise is the fast end of
deflateTT, 21.4 ms against 0.83 ms.

| File | deflateTT | zlib -1 | zlib -6 | zlib -9 |
| --- | ---: | ---: | ---: | ---: |
| Zeros | 44.7 | 0.053 | 0.128 | 0.127 |
| One byte | 43.1 | 0.053 | 0.128 | 0.127 |
| Sentence | 43.8 | 0.055 | 0.166 | 0.167 |
| Edited sentence | 57.8 | 0.060 | 0.199 | 0.200 |
| Two periods | 46.8 | 0.054 | 0.126 | 0.135 |
| Two alphabets | 20.8 | 0.911 | 1.057 | 1.057 |
| Text then noise | 34.4 | 0.345 | 0.433 | 0.442 |
| Noise | 21.4 | 0.911 | 0.838 | 0.833 |
| Period 100 | 44.9 | 0.057 | 0.128 | 0.129 |
| Copy past 32 KiB | 21.0 | 0.785 | 0.834 | 0.833 |
| Lazy traps | 181 | 0.253 | 1.84 | 2.56 |
| Long matches | 14.3 | 0.026 | 0.099 | 0.322 |
| Sources | 131 | 0.693 | 2.50 | 5.85 |

Times are milliseconds per call.

inflateTT-MP's byte decoder was timed on the same machine. The context and
the output buffer were allocated once, outside the timer. zlib used one
`z_stream` and `inflateReset` between calls. Each figure is uncompressed
megabytes per second over a loop of at least 300 ms. The files are 128 KiB
versions of the same generators, except the source snapshot above and a
1 MiB noise file. Zipf draws each byte independently with weight
`256 / (k + 1)`. The 4096-byte period is `(i * 17 + 3) & 255`, repeated.
The far copy plants 8192 bytes from the start at offset 100000.

On a repeated sentence the speeds are 2,068 MB/s for DEFd and 2,656 MB/s
for zlib level 9, at almost the same size (447 and 449 bytes). A stored
block is the other way around: 7,194 MB/s against 4,359 MB/s, and the DEFd
file is 13 bytes smaller. Zipf is the slow decode, 132 MB/s against 262 MB/s,
and DEFd is the smaller file (100,378 against 105,627).

| File | Raw | DEFd | MB/s | zlib -1 | MB/s | zlib -9 | MB/s |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Zeros | 131,072 | 147 | 789 | 595 | 2,594 | 149 | 658 |
| One byte | 131,072 | 148 | 789 | 596 | 2,632 | 151 | 656 |
| Sentence | 131,072 | 447 | 2,068 | 837 | 2,603 | 449 | 2,656 |
| Period 4096 | 131,072 | 838 | 1,861 | 1,358 | 2,527 | 840 | 2,566 |
| Period 100 | 131,072 | 612 | 2,006 | 1,074 | 2,568 | 613 | 2,446 |
| Zipf | 131,072 | 100,378 | 132 | 107,700 | 250 | 105,627 | 262 |
| Noise | 131,072 | 131,085 | 7,194 | 131,118 | 4,369 | 131,118 | 4,359 |
| Copy past 32 KiB | 131,072 | 131,085 | 7,166 | 131,118 | 4,374 | 131,118 | 4,366 |
| Sources | 125,833 | 27,063 | 361 | 35,850 | 690 | 28,374 | 1,176 |
| Noise, 1 MiB | 1,048,576 | 1,048,589 | 6,214 | 1,048,902 | 4,301 | 1,048,902 | 4,295 |

The sources were originally hosted at <https://www.villehelin.com/nintendoDS.html>
and moved to GitHub in 2021.
