#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <arpa/inet.h>
#define MAX_LOAD_FACTOR .66
/* This isn't the growth factor that causes big muscles (or possibly cancer), this is how much our stack expands when it reaches its capacity */
#define GROWTH_FACTOR 1.5
typedef struct _kmer_pointer{
  /* the sequence_number tells us which sequence the kmer is in. The amino_acid_index tells us where the kmer is within the sequence */
  size_t sequence_number;
  size_t amino_acid_index;
  /* location_pointer points to the exact start of the kmer */
  char* location_pointer;
  /* sequence_start_pointer points to the start of the sequence that contains the kmer */
  char* sequence_start_pointer;
} KMerPointer;

typedef struct _stack{
  KMerPointer* elements;
  size_t capacity;
  size_t current_size;
} Stack;

Stack* init_stack(size_t);


int add_to_stack(Stack*, KMerPointer);

typedef struct node{
  char* sequence;
  int num_characters;
  unsigned long count;
  //the sequences that contain the kmer
  Stack* sequences_contain_kmer;
  struct node *nextNode;
} Node;

typedef struct _linkedList {
    Node* start;
    Node* end;
} LinkedList;



typedef struct _hashTable{
    unsigned long num_buckets;
    unsigned long num_entries;
    LinkedList* lists;
} HashTable;
Node* getNode(HashTable*, char*);
typedef struct _fileIndex{
  off_t hash_table_area;
  off_t linked_list_area;
  off_t stack_area;
  unsigned int fragment_size;
} FileIndex;


int hash(char* sequence, int num_characters);



HashTable* doubleSize(HashTable* table);


HashTable* make_table(unsigned long num_buckets);


void print_and_free_table(HashTable*, FILE*);






void add_node(Node* node, HashTable* table);




void increment_count(char*, int, size_t, size_t, char*, char*, HashTable*);

HashTable* doubleSize(HashTable* table);

