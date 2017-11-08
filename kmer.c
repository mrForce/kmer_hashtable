#define _FILE_OFFSET_BITS 64

#include<stdio.h>
#include<stdlib.h>
#include<errno.h>
#include <pthread.h>
#include "hash_table.h"
#include<string.h>
#include<mcheck.h>
#define LINE_LENGTH 256
#define BUFFER_CAPACITY 50

#define NUM_BUCKETS 1021

typedef struct section{
    //inclusive
    int min_sum;
    //inclusive
    int max_sum;
    HashTable* table;
    pthread_mutex_t lock;
    struct section* previous_section;
    struct section* next_section;
} Section;


typedef struct {
    //points to the block of memory that the thread should work with; this contains the sequences, seperated by \0 between each FASTE entry. 
    char* nucleotide_block;
    int k;
  unsigned int* weights;
  size_t first_sequence_index;
  size_t first_amino_acid_index;
    size_t num_sequences;
    Section* sectionList;
    int num_sections;
} KMerPayload;




unsigned int* distribute(unsigned int num_items, int num_buckets);
unsigned long long* distribute_nucleotides(unsigned long long num_nucleotides, int num_threads);
//what each thread does
void* mer_count(void* arg);



int compare_numbers (const void * x, const void * y) {
  //do it in descending order
   return ( *(unsigned long long*)y - *(unsigned long long*)x );
}
/*
  Takes in path to a FASTA file, k, and number of threads.
 */
int main(int argc, char* argv[]){
  mtrace();
    if(argc == 4){
	int num_threads = atoi(argv[3]);
	int k = atoi(argv[2]);
	char* fasta_file_path = argv[1];
	FILE* fasta_file = fopen(fasta_file_path, "r");
	char line[LINE_LENGTH];
	/*
	  I'm going to make two passes through the file.

	  First, I'm going to count the number of FASTA sequences, and the number of nucleotides in each. I will use that information to distribute the sequences to the threads (we want to give each thread a similar number of sequences).

	  We may still need to make a few realloc calls, but I'm not worried about that, since we're only storing integers. I will assume that the length of a FASTA sequence can be stored in an unsigned long int. 

	  Then I will allocate a whole block of memory for each of the threads. 
	  
	  I will then make a second pass through the file, and 
	 */
	unsigned int num_fasta_sequences = 0;
	
	//use initial size of 10 sequences. We'll double the length whenever we need to using realloc.
        unsigned long* sequence_lengths = (unsigned long*) malloc(sizeof(unsigned long)*10);
	//how many sizes can be stored before we need to realloc.
	unsigned int capacity = 10;
	//total number of nucleotides.

	unsigned long long num_nucleotides = 0;
	int i;
	//I'd like to allocate a whole block of memory for the FASTA file, so let's 
	while(fgets(line, sizeof(line), fasta_file) != NULL){
	    //if line starts with '>', then it's the header to a fasta file.
	    if(line[0] == '>'){
		/*
		  If the sequence length is less than k, then don't include it in the number of nucleotides.
		 */
		if(num_fasta_sequences > 0 && sequence_lengths[num_fasta_sequences] < k){
		    num_nucleotides -= sequence_lengths[num_fasta_sequences];
		}
		if(capacity == num_fasta_sequences){
		    //then realloc -- since it's an array of unsigned longs, this shouldn't be too expensive.
		    sequence_lengths = realloc(sequence_lengths, 2*capacity*sizeof(unsigned long));
		    capacity = capacity*2;
		    
		    if(sequence_lengths == NULL){
			//then realloc failed.
			//if there isn't enough memory to do a realloc, then our k-mer counting attempt will surely fail. Let's just exit at this point.
			printf("Realloc failed while parsing the FASTA file. Exiting\n");
			return -1;
		    }
		}




		sequence_lengths[num_fasta_sequences] = 0;
		  
		

	    }else{
		//need to count number of nucleotides in the line
		i = 0;
	        while(line[i] >= 'A' && line[i] <= 'Z'){
		    i++;
		}
	        sequence_lengths[num_fasta_sequences] += i;
		num_nucleotides += i;
		num_fasta_sequences++;
	    }
	}

	//make sure the last sequence is of at least length k.
	if(sequence_lengths[num_fasta_sequences - 1] < k){
	    num_nucleotides -= sequence_lengths[num_fasta_sequences - 1];
	}


	//distribute the nucleotides over the threads
	unsigned long long* nucleotide_distribution = distribute_nucleotides(num_nucleotides, num_threads);
	//copy the distribution to use for when we are actually copying the sequences from the file
	unsigned long long* nucleotide_distribution_copy = (unsigned long long*) malloc(sizeof(unsigned long long)*num_threads);
	memcpy(nucleotide_distribution_copy, nucleotide_distribution, sizeof(unsigned long long)*num_threads);
	unsigned long last_sequence_index = 0;
	unsigned long long partial_nucleotide_index = 0;
	unsigned long long  space_to_reserve;
	size_t num_sequences;


	KMerPayload* payloads = (KMerPayload*) malloc(sizeof(KMerPayload)*num_threads);
	for(i = 0; i < num_threads; i++){
	    space_to_reserve = 0;
	    num_sequences = 0;
	    payloads[i].first_sequence_index = last_sequence_index;
	    payloads[i].first_amino_acid_index = partial_nucleotide_index;
	    while(nucleotide_distribution[i] > 0){
		if(sequence_lengths[last_sequence_index] - partial_nucleotide_index > nucleotide_distribution[i]){
		    //then part of a sequence. Keep last sequence index here.
		    /*
		      We need to reserve an additional char for the last \0

		      We need to reserve k - 1 chars to deal with the fact that we're splitting a sequence between two or more threads.  

		      k - 1 + 1 = k.
		     */
		    space_to_reserve += nucleotide_distribution[i] + k;
		    partial_nucleotide_index += nucleotide_distribution[i];
		    nucleotide_distribution[i] = 0;
		   
		}else{
		    //then this thread finishes off the sequence.
		    
		    space_to_reserve += (sequence_lengths[last_sequence_index] - partial_nucleotide_index + 1);
		    nucleotide_distribution[i] -= (sequence_lengths[last_sequence_index] - partial_nucleotide_index);
		    last_sequence_index++;
		    partial_nucleotide_index = 0;
		}
		num_sequences++;
	    }
	    payloads[i].nucleotide_block = (char*) malloc(sizeof(char)*space_to_reserve);
	    payloads[i].num_sequences = num_sequences;
	}
	  

	printf("Memory allocated for storing sequences\n");
	//now we have the memory allocated -- make another pass through the file
	rewind(fasta_file);

	unsigned int sequence_index = 0;
	//this should be the index of where we can start writing to!
	unsigned long long nucleotide_index = 0;
	unsigned long long temp_nucleotide_index = 0;
	/*
	  Get count of nucleotides in FASTA file
	*/
	
	//set equal to 1 to avoid issues when  a certain peptide doesn't exist
	unsigned long long peptide_count[26];
	for(i = 0; i < 26; i++){
	  peptide_count[i] = 1;
	}
	fgets(line, sizeof(line), fasta_file);
	free(nucleotide_distribution);
	unsigned long long m = 0;
	unsigned long long m_temp = 0;
	for(i = 0; i < num_threads; i++){
	    nucleotide_index = 0;
	    while(nucleotide_distribution_copy[i] > 0){
		
		if(sequence_lengths[sequence_index] >= k){
		  while(nucleotide_distribution_copy[i] > 0 && (line[m] >= 'A' && line[m] <= 'Z')){
					//do the counts here.... increment m, decrement nucleotide_distribution_copy[i]
			//don't forget to increment nucleotide_index.
		    peptide_count[line[m] - 'A']++;
		    payloads[i].nucleotide_block[nucleotide_index] = line[m];
		    nucleotide_distribution_copy[i]--;
		    nucleotide_index++;
		    m++;
		    
		  }
		}
		/*
		  Get the new line. Detect if it's a FASTA header.

		  If it is. then add a \0 to the nucleotide block, and increment the nucleotide_index. 
		  
		 */
		if(nucleotide_distribution_copy[i] == 0){
		  if(line[m] >= 'A' && line[m] <= 'Z'){
			//then stopping in middle of line. Need to move to next payload, and add k-1 part, and \0.
			temp_nucleotide_index = nucleotide_index;
			m_temp = m;
			while(temp_nucleotide_index - nucleotide_index < k - 1 && (line[m] >= 'A' && line[m] <= 'Z')){
			    payloads[i].nucleotide_block[temp_nucleotide_index] = line[m_temp];
			    m_temp++;
			    temp_nucleotide_index++;
			}
			//	m++;
			payloads[i].nucleotide_block[temp_nucleotide_index] = '\0';
		    }else{
			//then at end of the line.
			m = 0;
			fgets(line, sizeof(line), fasta_file);
		        /*
			  If at end of the sequence, increment the sequence index, and move past the FASTA header.
			 */
			if(line[0] == '>'){
			    payloads[i].nucleotide_block[nucleotide_index] = '\0';
			    while(line[0] == '>'){
				//we use a while loop, just incase there are two consective FASTA header lines.
				
				sequence_index++;
				fgets(line, sizeof(line), fasta_file);
			    }

			}
			 
		
		    }
		}else{
		    //we know that we are at the end of the line. Check if we're at the end of the sequence.
		    m = 0;
		    fgets(line, sizeof(line), fasta_file);
		    if(line[0] == '>'){
			//then place in a \0 at the end
			payloads[i].nucleotide_block[nucleotide_index] = '\0';
			nucleotide_index++;
			while(line[0] == '>'){
			    sequence_index++;
			    fgets(line, sizeof(line), fasta_file);
			}
		    }
		}
	    }
	
	}
	free(sequence_lengths);
	free(nucleotide_distribution_copy);
	printf("Copied the sequences from the FASTA file\n");
	fclose(fasta_file);
	//We need to assign weights to the nucleotides, where the most frequent nucleotide has a weight of 1, second most frequent has weight of 2, then 3, then 4
	unsigned long long counts[26];
	for(int i = 0; i < 26; i++){
	  counts[i] = peptide_count[i];
	}
	qsort(counts, 26, sizeof(unsigned long long), compare_numbers);
	unsigned int weights[26];
	for(int i = 0; i < 26; i++){
	  weights[counts[i] - 'A'] = i + 1;
	}
	int num_sections = 2*num_threads;
	/*
	  The idea behind sections (I'm adding this on November 8th, 2017):

	  Basically, each polypeptide has a score associated with it (which is the sum of the weights of the peptides in it)
	  We section up the hash table into sections; each section has a range of sums. 	  
	 */
	unsigned int* section_distribution = distribute(26*k, num_sections);
	unsigned int last_max_sum = k;

	Section* sections = (Section*) malloc(num_sections*sizeof(Section));



	for(i = 0; i < 2*num_threads; i++){
	    sections[i].min_sum = last_max_sum;
	    sections[i].max_sum = section_distribution[i] + last_max_sum;
	    last_max_sum = sections[i].max_sum;
	    sections[i].table = make_table(NUM_BUCKETS);	    
	    pthread_mutex_init(&sections[i].lock, NULL);
	    //sections[i].lock = lock;
	    if(i == 0){
		sections[0].previous_section = NULL;
	    }else{
		sections[i].previous_section = &sections[i - 1];
	        sections[i - 1].next_section = &sections[i];
	    }
	}
	free(section_distribution);
	sections[2*num_threads - 1].next_section = NULL;
	//now give the sections to the payloads
	for(i = 0; i < num_threads; i++){
	    payloads[i].sectionList = sections;
	    payloads[i].k = k;
	    payloads[i].weights = &weights[0];
	    payloads[i].num_sections = 2*num_threads;
	}
	printf("Going to create threads\n");
	pthread_t* threads = (pthread_t*) malloc(sizeof(pthread_t)*num_threads);
	//now that we have the payloads set up, need to make the threads
	for(i = 0; i < num_threads; i++){
	    pthread_create(&threads[i], NULL, mer_count, (void*) (&payloads[i]));
	}
	//wait until the threads finish.
	printf("Waiting for threads to finish...\n");
	for(i = 0; i < num_threads; i++){
	    pthread_join(threads[i], NULL);
	}

	HashTable* tempTable;


	//now take the sections, and print out their hashtables.
	for(i = 0; i < num_sections; i++){
	    tempTable = sections[i].table;
	    print_and_free_table(tempTable);
	    //now need to destroy the mutex

	    pthread_mutex_destroy(&sections[i].lock);
	    
	    }
	
	for( i = 0; i < num_threads; i++){
	    free(payloads[i].nucleotide_block);
	    
	}
	    

	free(threads);

	free(payloads);

	free(sections);

    }else{
	printf("Usage: kmer <FASTA file> <k> <number of threads> \n");
    }
    muntrace();
    return 0;
}

unsigned int* distribute(unsigned int num_items, int num_buckets){
    unsigned int* distribution = (unsigned int*) malloc(sizeof(unsigned int)*num_buckets);
//As a base amount, each gets num_items/num_buckets. 
    unsigned int base_amount = num_items/num_buckets;
    //number left to distribute.
    unsigned int left_over = num_items % base_amount;
    int i;
    for(i = 0; i < num_buckets; i++){
	if(left_over > 0){
	    distribution[i] = base_amount + 1;
	    left_over--;
	}else{
	    distribution[i] = base_amount;
	}
    }
    return distribution;
}

/*
  This is used to balance the nucleotides across the threads.

  The human genome is a little over 3 billion base pairs. Since some creatures have genomes MUCH larger than the human genome, so I will use an unsigned long long
 */

unsigned long long* distribute_nucleotides(unsigned long long num_nucleotides, int num_threads){
    unsigned long long* distribution = (unsigned long long*) malloc(sizeof(unsigned long long)*num_threads);
    unsigned long long base_amount = num_nucleotides/num_threads;
    unsigned long long left_over = num_nucleotides % base_amount;
    int i;
    for(i = 0; i < num_threads; i++){
	if(left_over > 0){
	    distribution[i] = base_amount + 1;
	    left_over--;
	}else{
	    distribution[i] = base_amount;
	}
    }
    return distribution;
}

//returns the sum of the window. If we hit the end of the sequence (that is, a '\0', then we return -1;
int calculateSum(char* window_start, int k, unsigned int* weights){
    int i;
    int sum = 0;
    for(i = 0; i < k; i++){
	if(window_start[i] == '\0'){
	    return -1;
	}else{
	  sum += weights[window_start[i] - 'A'];
	}
	
    }
    return sum;
}
/*
  This gets the section for the given sum
 */
Section* getSection(Section* sectionList, int sum){
    Section* section = sectionList;
    int min_sum;
    int max_sum;
    
    while(section != NULL){
	min_sum = section->min_sum;
	max_sum = section->max_sum;
	if(min_sum <= sum && max_sum >= sum){
	    return section;
	}
	section = section->next_section;
    }
    return section;
}



void* mer_count(void* arg){
    KMerPayload* payload = (KMerPayload*) arg;
    char* nucleotide_block = payload->nucleotide_block;
    int k = payload->k;
    size_t sequence_index = payload->first_sequence_index;
    int beginning = 1;
    if(sequence_index > 0){
      beginning = 0;
    }
    size_t amino_acid_index = payload->first_amino_acid_index;
    size_t num_sequences = payload->num_sequences;
    Section* sectionList = payload->sectionList;
    if(beginning){
      printf("sequence index: %zu, amino acid index: %zu, num sequences: %zu\n", sequence_index, amino_acid_index, num_sequences);
    }
    //we need a window of size k
    char* window = nucleotide_block;
    int window_sum, min_sum, max_sum;


    Section* temp_section;
    
    char *buffer[BUFFER_CAPACITY];
    size_t sequence_indices[BUFFER_CAPACITY];
    size_t amino_acid_indices[BUFFER_CAPACITY];
    int num_elements_in_buffer = 0;
    unsigned int* weights = payload->weights;



    int i, j;
    for(i = 0; i < num_sequences; i++){
      
	window_sum = calculateSum(window, k, weights);
	
	if(window_sum > 0){
	    temp_section = getSection(sectionList, window_sum);
	    min_sum = temp_section->min_sum;
	    max_sum = temp_section->max_sum;

	    while(window[k - 1] != '\0'){	    
	      if(beginning){
		printf("num elements in buffer: %i\n", num_elements_in_buffer);
		printf("window: %s\n", window);
	      }
		if(num_elements_in_buffer == BUFFER_CAPACITY|| window_sum < min_sum || window_sum > max_sum){
		    //commit the buffer!
		    pthread_mutex_lock(&(temp_section->lock));
		    for(j = 0; j < num_elements_in_buffer; j++){
		      printf("Commiting from buffer: %s\n", buffer[j]);
		      increment_count(buffer[j], k, sequence_indices[j], amino_acid_indices[j], temp_section->table);
		    }
		    
		    pthread_mutex_unlock(&(temp_section->lock));
		    num_elements_in_buffer = 0;

		}
		if(window_sum < min_sum || window_sum > max_sum){
		    temp_section = getSection(sectionList, window_sum);
		    min_sum = temp_section->min_sum;
		    max_sum = temp_section->max_sum;
		    //buffer[num_elements_in_buffer] = window;
		    //num_elements_in_buffer++;
		    //window++;
		}else{
		    //if we're in the same section, then move window forward by doing this stuff.
		    buffer[num_elements_in_buffer] = window;
		    sequence_indices[num_elements_in_buffer] = sequence_index + i;
		    amino_acid_indices[num_elements_in_buffer] = amino_acid_index;
		    if(window[0] >= 'A' && window[0] <= 'Z'){
		      window_sum -= weights[window[0] - 'A'];
		    }
		    if(window[k] >= 'A' && window[k] <= 'Z'){
		      window_sum += weights[window[k] - 'A'];
		    }
		    /*
		      Just to clear up a concern that I had, but I don't think I need to worry about:

		      Suppose that window[k] == '\0'. Then the window_sum gets out of whack; but that doesn't matter, because we'll hae to recaclulate the window sum anyway, since we're at a new sequence.
		     */
		    window++;
		    amino_acid_index++;
		    num_elements_in_buffer++;
		}

		    

	    }
	    if(i < num_sequences - 1){
		//if we aren't at the last sequence, be sure to push the window over
		window += strlen(window) + 1;
		amino_acid_index = 0;
	    }
	}else{
	    window += strlen(window) + 1;
	}




    }
    if(num_elements_in_buffer > 0){
	 //commit the buffer!
	 pthread_mutex_lock(&(temp_section->lock));
	 for(i = 0; i < num_elements_in_buffer; i++){
	   printf("Commiting from buffer: %s\n", buffer[j]);
	   increment_count(buffer[i], k, sequence_indices[i], amino_acid_indices[i], temp_section->table);
	 }
	 pthread_mutex_unlock(&(temp_section->lock));
	 num_elements_in_buffer = 0;

    }

    pthread_exit(NULL);
    
    
}
