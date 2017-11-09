CC = gcc
CFLAGS = -Wall -O3
DEPS = hash_table.h 
OBJ = kmer.o hash_table.o

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

kmer: $(OBJ) 
	$(CC) $(CFLAGS) -o $@ $^ -lpthread
clean:
	rm *.o
