
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#define MAX_LOAD_FACTOR .66

typedef struct node{
    char* sequence;
    int num_characters;
    unsigned long count;
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

int hash(char* sequence, int num_characters);



HashTable* doubleSize(HashTable* table);


HashTable* make_table(unsigned long num_buckets);


void print_and_free_table(HashTable* table);






void add_node(Node* node, HashTable* table);




void increment_count(char* sequence, int num_characters, HashTable* table);

HashTable* doubleSize(HashTable* table);

