CC = gcc
CFLAGS = -Wall
EXEC = a.out
LIBS = -lm

SRC = $(wildcard *.c)
OBJ = $(SRC:.c=.o)

all: $(EXEC)

%.o: %.c
	$(CC) $(CFLAGS) -c $<

$(EXEC): $(OBJ)
	$(CC) -o $@ $^ $(LIBS)

clean:
	rm -f $(OBJ) $(EXEC)
