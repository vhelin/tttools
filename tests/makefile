CC = gcc
CFLAGS = -Wall -Wextra -O2 -ansi -pedantic -I../deflateTT -I../inflateTT-MP

.PHONY: all clean

all: ../tests.exe ds.o
	../tests.exe

tests.exe: ../tests.exe

../tests.exe: roundtrip.c ../deflateTT/main.c ../deflateTT/huffman.c ../inflateTT-MP/inflate.c \
           ../deflateTT/deflate.h ../deflateTT/huffman.h ../inflateTT-MP/inflate.h
	$(CC) $(CFLAGS) -DDEFLATE_LIBRARY -o ../tests.exe roundtrip.c \
	  ../deflateTT/main.c ../deflateTT/huffman.c ../inflateTT-MP/inflate.c

ds.o: ../inflateTT-MP/inflate_ds.c ../inflateTT-MP/inflate_ds.h \
      ../inflateTT-MP/inflate.c ../inflateTT-MP/inflate.h
	$(CC) $(CFLAGS) -c ../inflateTT-MP/inflate_ds.c -o ds.o

clean:
	rm -f ../tests.exe ds.o roundtrip roundtrip.exe
