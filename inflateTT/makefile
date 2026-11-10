CC = gcc
CFLAGS = -Wall -Wextra -O2 -ansi -pedantic
INFLATE = ../inflateTT-MP

.PHONY: all clean install

all: inflateTT

inflateTT: main.o inflate.o
	$(CC) -o inflateTT main.o inflate.o

main.o: main.c $(INFLATE)/inflate.h
	$(CC) $(CFLAGS) -I$(INFLATE) -c main.c

inflate.o: $(INFLATE)/inflate.c $(INFLATE)/inflate.h
	$(CC) $(CFLAGS) -I$(INFLATE) -c $(INFLATE)/inflate.c -o inflate.o

clean:
	rm -f main.o inflate.o inflateTT inflateTT.exe

install: inflateTT
	cp inflateTT /usr/local/bin/
