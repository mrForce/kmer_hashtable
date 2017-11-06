#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <arpa/inet.h>
#define MAX_LOAD_FACTOR .66
/* This isn't the growth factor that causes big muscles (or possibly cancer), this is how much our stack expands when it reaches its capacity */
#define GROWTH_FACTOR 1.5
typedef struct _kmer_pointer{
  /* the sequence_number tells us which sequence the kmer is in. The amino_acid_index tells us where the kmer is within the sequence */
  uint32_t sequence_number;
  uint32_t amino_acid_index;
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

typedef struct _fileIndex{
  off_t hash_table_area;
  off_t linked_list_area;
  off_t stack_area;
  unsigned int fragment_size;
} FileIndex;


typedef struct _serial_stack{
  KMerPointer* elements;
  uint32_t num_elements;
} SerializedStack;



typedef struct _serial_node{
  char* sequence;
  off_t sequences_stack;
} SerializedNode;


typedef struct _nodeArray{
  //location of the array within the file
  off_t nodes;
  size_t num_nodes;
} SerializedNodeArray;

typedef struct _serial_hash_table{
  off_t node_arrays[NUM_BUCKETS];
} SerializedHashTable;


SerializedStack read_serialized_stack(off_t, FILE*);

/* 
Returns 0 if successfull, -1 otherwise.

Moves the file pointer right after the stack.
 */
char write_serialized_stack(SerializedStack, FILE*);

//returns 1 if successfull, 0 otherwise
int saveToDisk(HashTable*, FILE*);
int hash(char* sequence, int num_characters);



HashTable* doubleSize(HashTable* table);


HashTable* make_table(unsigned long num_buckets);


void print_and_free_table(HashTable* table);






void add_node(Node* node, HashTable* table);




void increment_count(char* sequence, int num_characters, size_t sequence_index, size_t amino_acid_index, HashTable* table);

HashTable* doubleSize(HashTable* table);

