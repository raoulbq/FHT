CC   = gcc
CFLAGS = -Wall -Werror -pedantic

OBJ  = fht.o lalgebra.o
LIBS = -lm -lfftw3
BIN  = fht

RM = rm -f


.PHONY: all clean

all: fht

clean:
	${RM} $(OBJ) $(BIN)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

$(BIN): $(OBJ)
	$(CC) $(OBJ) -o $(BIN) $(LIBS)
