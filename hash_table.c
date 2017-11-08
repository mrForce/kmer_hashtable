#include "hash_table.h"
#include<mcheck.h>
#define INIT_STACK_SIZE 5
/*
  This implements a stack that can only be added to. We do this like a dynamic array
 */
Stack* init_stack(size_t initial_capacity){
  Stack* stack = (Stack*) malloc(sizeof(Stack));
  KMerPointer* elements = (KMerPointer*) malloc(sizeof(KMerPointer)*initial_capacity);
  stack->elements = elements;
  stack->capacity = initial_capacity;
  stack->current_size = 0;
  return stack;
}



void free_stack(Stack* stack){

  free(stack->elements);
  free(stack);
}

//returns 1 if successful, 0 otherwise
int add_to_stack(Stack* stack, KMerPointer kmer_pointer){
  if(stack->capacity == stack->current_size){
    size_t new_capacity = 1.5*stack->capacity;
    KMerPointer* new_stack_elements = realloc(stack->elements, sizeof(KMerPointer)*new_capacity);
    if(new_stack_elements){
      stack->elements = new_stack_elements;
      stack->capacity = new_capacity;
    }else{
      return 0;
    }
  }
  memcpy(&(stack->elements[stack->current_size++]), &(kmer_pointer), sizeof(KMerPointer));
  return 1;
}


int hash(char* sequence, int num_characters){
    int hash = 0;
    int i;
    for(i = 0; i < num_characters; i++){
	//multiply by two
	hash = hash << 1;
	//do an xor with the character
	hash = hash ^ sequence[i];
    }
    return hash;
}

HashTable* make_table(unsigned long num_buckets){
    
    HashTable* table = (HashTable*) malloc(sizeof(HashTable));
    table->num_buckets = num_buckets;
    table->num_entries = 0;
    LinkedList* lists = (LinkedList*) malloc(sizeof(LinkedList)*num_buckets);
    unsigned long i;
    for(i = 0; i < num_buckets; i++){
	lists[i].start = NULL;
	lists[i].end = NULL;
    }
    table->lists = lists;
    return table;
}


/*
  Prints the table, and then frees it.
 */

void print_and_free_table(HashTable* table){
  int i;
  size_t j;
    Node* node;
    Node* nextNode;
    LinkedList* lists = table->lists;
    int z = 0;
    for(i = 0; i < table->num_buckets; i++){
	node = lists[i].start;
	if(node != NULL){
	  z = 0;
	    while(node != NULL){
		printf("K-mer: %s, count: %lu, ", node->sequence, node->count);
		for(j = 0; j < node->sequences_contain_kmer->current_size; j++){
		  printf("(%lu, %lu), ", node->sequences_contain_kmer->elements[j].sequence_number, node->sequences_contain_kmer->elements[j].amino_acid_index);
		}
		printf("\n");
		nextNode = node->nextNode;
		free_stack(node->sequences_contain_kmer);
		free(node->sequence);
		free(node);
		node = nextNode;
		z++;
	    }
	}

    }
    free(table->lists);
    free(table);
}







void add_node(Node* node, HashTable* table){
    node->nextNode = NULL;
    fflush(stdin);
    int node_hash = hash(node->sequence, node->num_characters);
    int bucket = node_hash % table->num_buckets;
    //we have the bucket index, so add to bucket.
    LinkedList* list = &table->lists[bucket];
    node->nextNode = NULL;
    Node* endNode = list->end;
    if(endNode == NULL){
	list->start = node;
    }else{
	endNode->nextNode = node; 
    }
    list->end = node;
}


//don't need to remove from hash table.

//since we're taking a substring of a sequence, 
void increment_count(char* sequence, int num_characters, size_t sequence_index, size_t amino_acid_index, HashTable* table){
    if(((float) table->num_entries)/(table->num_buckets) >= MAX_LOAD_FACTOR){
	//double the size of the table
	table = doubleSize(table);
    }
    int node_hash = hash(sequence, num_characters);
    LinkedList* list = &table->lists[node_hash % table->num_buckets];
    Node* node = list->start;
    int i;
    if(node == NULL){
	node = (Node*) malloc(sizeof(Node));

	node->count = 1;
	node->nextNode = NULL;
	node->num_characters = num_characters;
	node->sequence = (char*) malloc(sizeof(char)*(num_characters + 1));
	if(node->sequence == NULL){
	  printf("could not reserver memory\n");
	}
	for(i = 0; i < num_characters; i++){
	    node->sequence[i] = sequence[i];
	}
	node->sequence[num_characters] = '\0';
	node->sequences_contain_kmer = init_stack(INIT_STACK_SIZE);	
	KMerPointer pointer;
	pointer.amino_acid_index = amino_acid_index;
	pointer.sequence_number = sequence_index;
	if(!add_to_stack(node->sequences_contain_kmer, pointer)){
	  printf("Problem with increment count\n");
	}
	list->end = node;
	list->start = node;
    }else{

	char keepGoing = 1;

	Node* newNode;
	while(keepGoing){
	    if(memcmp(node->sequence, sequence, num_characters*sizeof(char)) == 0){
		//then the same
	      node->count++;
		KMerPointer pointer;
		pointer.amino_acid_index = amino_acid_index;
		pointer.sequence_number = sequence_index;
		if(!add_to_stack(node->sequences_contain_kmer, pointer)){
		  printf("Problem with increment count\n");
		}
		keepGoing = 0;
	    }else{
		if(node->nextNode == NULL){
		    newNode = (Node*) malloc(sizeof(Node));
		    newNode->count = 1;
		    newNode->nextNode = NULL;
		    newNode->num_characters = num_characters;
		    newNode->sequence = (char*) malloc(sizeof(char)*(num_characters + 1));
		    if(newNode->sequence == NULL){
		      printf("could not reserver memory\n");
		    }
		    for(i = 0; i < num_characters; i++){
			newNode->sequence[i] = sequence[i];			
		    }
		    newNode->sequence[num_characters] = '\0';
		    list->end = newNode;
		    
		    newNode->sequences_contain_kmer = init_stack(INIT_STACK_SIZE);
		    KMerPointer pointer;
		    pointer.amino_acid_index = amino_acid_index;
		    pointer.sequence_number = sequence_index;
		    if(!add_to_stack(newNode->sequences_contain_kmer, pointer)){
		      printf("Problem with increment count\n");
		    }
		    keepGoing = 0;
		    node->nextNode = newNode;
		}else{
		    node = node->nextNode;
		}
	    }
	}
    }
}


HashTable* doubleSize(HashTable* table){
    int i;
    Node* node;
    Node* nextNode;

    HashTable* newTable = make_table(table->num_buckets*2);
    LinkedList* lists = table->lists;
    for(i = 0; i < table->num_buckets; i++){
	node = lists[i].start;
	while(node != NULL){
	    add_node(node, newTable);
	    nextNode = node->nextNode;
	    free(node);
	    node = nextNode;
	}
    }
    free(table->lists);
    free(table);
    
    return newTable;
}
/*
HashTable* doubleSize(HashTable* table){
    //make a new table with double the number of buckets
    HashTable* newTable = make_table(table->num_buckets*2);
    int i;
    Node* node;
    LinkedList* list;
    //we need to loop through the buckets, and copy over all of the elements in them.
    for(i = 0; i < table->num_buckets; i++){	
	list = &table->lists[i];
	node = list->start;
	while(node != NULL){
	    add_node(node, newTable);
	    node = node->nextNode;
	}
        free(list);
    }
    free(table->lists);
    free(table);
    return newTable;
}
*/

