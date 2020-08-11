CC = mpicc 
CFLAGS = -O3 -finline-functions -fomit-frame-pointer -funroll-loops -Wall -pedantic -lm -D DEBUG 

HEADERS = $(wildcard *.h)
SOURCES = $(wildcard *.c)
OBJECTS= $(SOURCES:.c=.o)

exe = dot

all: .depend $(OBJECTS) $(exe)

dot: $(OBJECTS)
	$(CC) $(OBJECTS) -lm -o $(exe) 

.PHONY: clean 

clean:
	@-rm $(OBJECTS)
	@-rm -f *~ .depend $(exe) 

.depend:
	$(CC) -MM $(SOURCES) > .depend

ifeq ($(wildcard .depend), .depend)
include .depend
endif
# DO NOT DELETE
