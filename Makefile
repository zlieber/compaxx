
all: compaxx

SRC := compaxx.c extra.c test.c

OBJ := ${SRC:.c=.o}

compaxx: $(OBJ)
	gcc -o $@ $^ -lm

%.o: %.c %.h
	gcc -o $@ -g -c $<

test: compaxx
	./compaxx

clean:
	rm -f *.o

