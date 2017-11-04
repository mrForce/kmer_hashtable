#include "hash_table.h"



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
    int i;
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
    Node* node;
    Node* nextNode;
    LinkedList* lists = table->lists;
    for(i = 0; i < table->num_buckets; i++){
	node = lists[i].start;
	if(node != NULL){
	    while(node != NULL){
		printf("K-mer: %s, count: %lu \n", node->sequence, node->count);
		nextNode = node->nextNode;
		free(node->sequence);
		free(node);
		node = nextNode;
	    }
	}

    }
    free(table->lists);
    free(table);
}







void add_node(Node* node, HashTable* table){
    node->nextNode = NULL;
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
void increment_count(char* sequence, int num_characters, HashTable* table){
   
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
	for(i = 0; i < num_characters; i++){
	    node->sequence[i] = sequence[i];
	}
	node->sequence[num_characters] = '\0';
	list->end = node;
	list->start = node;
    }else{
	char keepGoing = 1;

	Node* newNode;
	while(keepGoing){
	    if(memcmp(node->sequence, sequence, num_characters*sizeof(char)) == 0){
		//then the same
		node->count++;
		keepGoing = 0;
	    }else{
		if(node->nextNode == NULL){
		    newNode = (Node*) malloc(sizeof(Node));
		    newNode->count = 1;
		    newNode->nextNode = NULL;
		    newNode->num_characters = num_characters;
		    newNode->sequence = (char*) malloc(sizeof(char)*(num_characters + 1));
		    for(i = 0; i < num_characters; i++){
			newNode->sequence[i] = sequence[i];			
		    }
		    node->sequence[num_characters] = '\0';
		    list->end = newNode;
		    node->nextNode = newNode;
		    keepGoing = 0;
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

