CC = gcc
CFLAGS = -Wall -Wextra -O2 -ansi -pedantic -c

.PHONY: all clean

all: inflate.o inflate_ds.o

inflate.o: inflate.c inflate.h
	$(CC) $(CFLAGS) inflate.c

inflate_ds.o: inflate_ds.c inflate_ds.h inflate.c inflate.h
	$(CC) $(CFLAGS) inflate_ds.c

clean:
	rm -f inflate.o inflate_ds.o
