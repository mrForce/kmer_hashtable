#include<stdio.h>
#include<stdlib.h>
#include<errno.h>
#include <pthread.h>
#include "hash_table.h"
#include<string.h>
#include<mcheck.h>
#define LINE_LENGTH 256
#define BUFFER_CAPACITY 50

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
    unsigned int a_weight;
    unsigned int g_weight;
    unsigned int c_weight;
    unsigned int t_weight;
  size_t first_sequence_index;
    size_t num_sequences;
    Section* sectionList;
    int num_sections;
} KMerPayload;




unsigned int* distribute(unsigned int num_items, int num_buckets);
unsigned long long* distribute_nucleotides(unsigned long long num_nucleotides, int num_threads);
//what each thread does
void* mer_count(void* arg);




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



		num_fasta_sequences++;
		sequence_lengths[num_fasta_sequences] = 0;
		

	    }else{
		//need to count number of nucleotides in the line
		i = 0;
	        while(line[i] == 'A' || line[i] == 'G' || line[i] == 'C' || line[i] == 'T'){
		    i++;
		}
	        sequence_lengths[num_fasta_sequences] += i;
		num_nucleotides += i;
	    }
	}

	//make sure the last sequence is of at least length k.
	if(sequence_lengths[num_fasta_sequences] < k){
	    num_nucleotides -= sequence_lengths[num_fasta_sequences];
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

	unsigned int sequence_index = 1;
	//this should be the index of where we can start writing to!
	unsigned long long nucleotide_index = 0;
	unsigned long long temp_nucleotide_index = 0;
	/*
	  Get count of nucleotides in FASTA file
	*/
	
	//set equal to 1 to avoid issues when  a certain nucleotide doesn't exist
	unsigned long long adenine_count = 1;
	unsigned long long guanine_count = 1;
	unsigned long long thymine_count = 1;
	unsigned long long cytosine_count = 1;
	fgets(line, sizeof(line), fasta_file);
	free(nucleotide_distribution);
	unsigned long long m = 0;
	unsigned long long m_temp = 0;
	for(i = 0; i < num_threads; i++){
	    nucleotide_index = 0;
	    while(nucleotide_distribution_copy[i] > 0){
		
		if(sequence_lengths[sequence_index] >= k){
		    while(nucleotide_distribution_copy[i] > 0 && (line[m] == 'A' || line[m] == 'G' || line[m] == 'C' || line[m] == 'T')){
					//do the counts here.... increment m, decrement nucleotide_distribution_copy[i]
			//don't forget to increment nucleotide_index.
			if(line[m] == 'A'){
			    adenine_count++;
			}else if(line[m] == 'C'){
			    cytosine_count++;
			}else if(line[m] == 'G'){
			    guanine_count++;
			}else if(line[m] == 'T'){
			    thymine_count++;
			}
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
		    if(line[m] == 'A' || line[m] == 'G' || line[m] == 'C' || line[m] == 'T'){
			//then stopping in middle of line. Need to move to next payload, and add k-1 part, and \0.
			temp_nucleotide_index = nucleotide_index;
			m_temp = m;
			while(temp_nucleotide_index - nucleotide_index < k - 1&& (line[m_temp] == 'A' || line[m_temp] == 'G' || line[m_temp] == 'C' || line[m_temp] == 'T')){
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
	unsigned long long counts[4];
	
	counts[0] = adenine_count;
	counts[1] = guanine_count;
	counts[2] = cytosine_count;
	counts[3] = thymine_count;
	
	//now find the maximum 
	int max_index;
	int t;
	unsigned int weights[4];
	for(t = 0; t < 4; t++){
	    max_index = 0;
	    for(i = 0; i < 4; i++){
		if(counts[i] > counts[max_index]){
		    max_index = i;
		}
		
	    }
	    weights[max_index] = t + 1;
	    counts[max_index] = 0;
	}
	int num_sections = 2*num_threads;
	//set up sections. 2*num_threads sections
	unsigned int* section_distribution = distribute(4*k, num_sections);
	unsigned int last_max_sum = k;

	Section* sections = (Section*) malloc(num_sections*sizeof(Section));
	for(i = 0; i < 2*num_threads; i++){
	    sections[i].min_sum = last_max_sum;
	    sections[i].max_sum = section_distribution[i] + last_max_sum;
	    last_max_sum = sections[i].max_sum;
	    sections[i].table = make_table(1021);
	    pthread_mutex_init(&sections[i].lock, NULL);
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
	    payloads[i].a_weight = weights[0];
	    payloads[i].g_weight = weights[1];
	    payloads[i].c_weight = weights[2];
	    payloads[i].t_weight = weights[3];
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
int calculateSum(char* window_start, int k, unsigned int a_weight, unsigned int g_weight, unsigned int c_weight, unsigned int t_weight){
    int i;
    int sum = 0;
    for(i = 0; i < k; i++){
	if(window_start[i] == '\0'){
	    return -1;
	}else
	if(window_start[i] == 'A'){
	    sum += a_weight;
	}else 
	if(window_start[i] == 'G'){
	    sum += g_weight;
	}else 
	if(window_start[i] == 'C'){
	    sum += c_weight;
	}else 
	if(window_start[i] == 'T'){
	    sum += t_weight;
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
    size_t num_sequences = payload->num_sequences;
    Section* sectionList = payload->sectionList;

    //we need a window of size k
    char* window = nucleotide_block;
    int window_sum, min_sum, max_sum;


    Section* temp_section;
    
    char *buffer[BUFFER_CAPACITY];
    size_t sequence_indices[BUFFER_CAPACITY];
    int num_elements_in_buffer = 0;
    unsigned int a_weight = payload->a_weight;
    unsigned int g_weight = payload->g_weight;
    unsigned int c_weight = payload->c_weight;
    unsigned int t_weight = payload->t_weight;



    int i, j;
    for(i = 0; i < num_sequences; i++){
      
	window_sum = calculateSum(window, k, a_weight, g_weight, c_weight, t_weight);
	
	if(window_sum > 0){
	    temp_section = getSection(sectionList, window_sum);
	    min_sum = temp_section->min_sum;
	    max_sum = temp_section->max_sum;

	    while(window[k - 1] != '\0'){
	      if(strlen(window) < 4){
		printf("window: %s\n", window);
	      }
		if(num_elements_in_buffer == BUFFER_CAPACITY|| window_sum < min_sum || window_sum > max_sum){
		    //commit the buffer!
		    pthread_mutex_lock(&(temp_section->lock));
		    for(j = 0; j < num_elements_in_buffer; j++){
		      increment_count(buffer[j], k, sequence_indices[j], 0, temp_section->table);
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
		    if(window[0] == 'A'){
			window_sum -= a_weight;
		    }else if(window[0] == 'G'){
			window_sum -= g_weight;
		    }else if(window[0] == 'C'){
			window_sum -= c_weight;
		    }else if(window[0] == 'T'){
			window_sum -= t_weight;
		    }
		    /*
		      Just to clear up a concern that I had, but I don't think I need to worry about:

		      Suppose that window[k] == '\0'. Then the window_sum gets out of whack; but that doesn't matter, because we'll hae to recaclulate the window sum anyway, since we're at a new sequence.
		     */
		    if(window[k] == 'A'){
			window_sum += a_weight;	
		    }else if(window[k] == 'G'){
			window_sum += g_weight;
		    }else if(window[k] == 'C'){
			window_sum += c_weight;
		    }else if(window[k] == 'T'){
			window_sum += t_weight;
		    }
		    window++;
		    num_elements_in_buffer++;
		}

		    

	    }
	    if(i < num_sequences - 1){
		//if we aren't at the last sequence, be sure to push the window over
		window += strlen(window) + 1;
	    }
	}else{
	    window += strlen(window) + 1;
	}




    }
    if(num_elements_in_buffer > 0){
	 //commit the buffer!
	 pthread_mutex_lock(&(temp_section->lock));
	 for(i = 0; i < num_elements_in_buffer; i++){
	   increment_count(buffer[i], k, sequence_indices[i], 0, temp_section->table);
	 }
	 pthread_mutex_unlock(&(temp_section->lock));
	 num_elements_in_buffer = 0;

    }

    pthread_exit(NULL);
    
    
}
