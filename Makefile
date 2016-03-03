# Project: FHT
CC   = gcc
CFLAGS = -Wall -Werror -pedantic

OBJ  = FHT.o lalgebra.o
LINKOBJ  = FHT.o lalgebra.o
LIBS = -lm -lfftw3
BIN  = FHT

RM = rm -f


.PHONY: all clean

all: FHT

clean:
	${RM} $(OBJ) $(BIN)

$(BIN): $(OBJ)
	$(CC) $(LINKOBJ) -o $(BIN) $(LIBS)

FHT.o: FHT.c
	$(CC) -c FHT.c -o FHT.o $(CFLAGS)

lalgebra.o: lalgebra.c
	$(CC) -c lalgebra.c -o lalgebra.o $(CFLAGS)
