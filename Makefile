CC   = gcc
CFLAGS = -Wall -Werror -pedantic

OBJ  = FHT.o lalgebra.o
LIBS = -lm -lfftw3
BIN  = FHT

RM = rm -f


.PHONY: all clean

all: FHT

clean:
	${RM} $(OBJ) $(BIN)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

$(BIN): $(OBJ)
	$(CC) $(OBJ) -o $(BIN) $(LIBS)
