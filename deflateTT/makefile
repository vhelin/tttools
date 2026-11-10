CC = gcc
CFLAGS = -Wall -Wextra -O2 -ansi -pedantic

.PHONY: all clean install

all: deflateTT

deflateTT: main.o huffman.o
	$(CC) -o deflateTT main.o huffman.o

main.o: main.c deflate.h huffman.h
	$(CC) $(CFLAGS) -c main.c

huffman.o: huffman.c huffman.h
	$(CC) $(CFLAGS) -c huffman.c

clean:
	rm -f main.o huffman.o deflateTT deflateTT.exe

install: deflateTT
	cp deflateTT /usr/local/bin/
