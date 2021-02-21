
CC = gcc
LD = gcc

CFLAGS = -Wall -c -g -O0 -ansi -pedantic
LDFLAGS = 

CFILES = main.c
HFILES = main.h
OFILES = main.o
EXECUT = inflateTT


all: $(OFILES) makefile
	$(LD) $(LDFLAGS) $(OFILES) -o $(EXECUT) -lm

main.o: main.c defines.h
	$(CC) $(CFLAGS) main.c


$(OFILES): $(HFILES)


clean:
	rm -f $(OFILES) core *~ $(EXECUT) gmon.out

nice:
	rm -f *~ gmon.out

install:
	make; cp $(EXECUT) /usr/local/bin/
