# Project: FHT
CC   = gcc
CFLAGS =

OBJ  = FHT.o lalgebra.o testingstuff.o
LINKOBJ  = FHT.o lalgebra.o testingstuff.o
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

testingstuff.o: testingstuff.c
	$(CC) -c testingstuff.c -o testingstuff.o $(CFLAGS)
