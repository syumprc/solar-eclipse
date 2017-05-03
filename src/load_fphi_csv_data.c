#include <stdio.h>
#include <pthread.h>
#include "safelib.h"
#include "csv.h"
#define BUFFER_SIZE 100

void append_phenotype_containers( char * id, double * trait, double * covariates, int n_traits, int n_covariates);
void append_ibdids_and_ids( char * id, int ibdid);

int is_empty(char * buffer);

#define N_PTHREADS 5

static int search_indices(char * field, char ** search_array, int * index_array, int index, int   n_found, int n){
	if(n_found == n)
		return 0;
		
	int i = 0; 
	while(i < n){
		if(index_array[i] == -1){
			if(!StringCmp(field, search_array[i], 0)){
				index_array[i] = index;
				return 1;
			}
		}
		i++;
	}
	
	return 0;
}
				
	
static int check_indices(int index, int * index_array, int  n_found, int n){
	if(n_found == n)
		return -1;
		
	int i = 0;
	while(i < n){
		if(index_array[i] == index){
			return i;
		 }
		 i++;
	}
	
	return -1;
	
}
char * clean_string(char * buffer);

static int comma_count(char * buffer){
	if(buffer == NULL)
		return 0;
	

	int count = 0;
	while(*buffer != '\0'){
		
		if(*buffer == ',' || *buffer == '\n')
			count++;
		buffer++;
	}
	
	return count;
}

static char * clean_cstring(char * str, int * n_lines){
	char *newstr = malloc(strlen(str)+1);

	char *write = newstr, *read = str;
	printf("%i \n", *n_lines);
	int lines = 0;
	do {
		if(*read == '\n')
			lines++;
		
		if (isgraph(*read) || *read == '\n')
			*write++ = *read;
			
			
	} while (*read++);	
	printf("%i \n", lines);
	*n_lines = lines;
	return newstr;
	
}		

struct pedindex_data{
	char * line;
	int id_index;
	int ibdid;
	
};

struct phenotype_data{
	
	char * line;
	int n_traits;
	int n_covs; 
	int * trait_indices;
	int * cov_indices;
	int  id_index;
	int thread_number;
};


void * read_pedindex_line(void * data){
	struct pedindex_data * ped_data = struct pedindex_data)data;
	
	
	
}

void * read_phenotype_line(void * data){
	
	struct phenotype_data  * pheno_data = ( struct phenotype_data*)data;
	
	char * token;
	int id_found = 0;
	int found_traits = 0;
	int found_covs = 0;
	int n_trts = pheno_data->n_traits;
	int n_covs = pheno_data->n_covs;
	int cov_array_index = 0;
	int trt_array_index  = 0;
	int index = 0;
	double * trait_values = malloc(sizeof(double)*n_trts);
	double * cov_values;// = mallo(sizeof(double)*n_covs);
	if(n_covs > 0){
		cov_values = (double*)malloc(sizeof(double)*n_covs);
	}	
	char * id;
	do{
		token = NextArg(pheno_data->line, ',', &pheno_data->line);	
		if(pheno_data->trait_indices[trt_array_index] == index){
			if(token == NULL){
				free(pheno_data->line);
				free(trait_values);
				if(n_covs > 0) free(cov_values);
				return NULL;
				}
				trait_values[trt_array_index] = strtod(token, NULL);
				index++;
				found_traits++;
				trt_array_index++;
				if(trt_array_index == n_trts)
					trt_array_index--;
			//	token = NextArg(line, ',', &line);
				continue;
			}
		//	array_index  = check_indices(index, cov_indices, found_covs, n_covs);
			
			if(n_covs != 0){
				if(pheno_data->cov_indices[cov_array_index] == index){
					if(token == NULL){
						free(pheno_data->line);
						free(cov_values);
						free(trait_values);
						return NULL;
					}
				
					if(!StringCmp(token, "M", 0)){
						cov_values[cov_array_index] = 1.0;
					}else if(!StringCmp(token, "F", 0)){
						cov_values[cov_array_index]= 2.0;
					}else{
						cov_values[cov_array_index] = strtod(token, NULL);
					}

					index++;
					found_covs++;
					cov_array_index++;
					
					if(cov_array_index == n_covs)
						cov_array_index--;					
		//		token = NextArg(line, ',', &line);
					continue;
				}	
			}		
			
			if(index == pheno_data->id_index){
				id = strdup(token);
			//	printf("%s \n", token);
				index++;
				id_found = 1;
		//		token = NextArg(line, ',', &line);
				continue;
			}

			
			index++;
		if(token != NULL)
			free(token);		
			
		}while(token != NULL && pheno_data->line != NULL && (found_traits != n_trts||found_covs != n_covs|| !id_found));	


			append_phenotype_containers(id, trait_values, cov_values, n_trts, n_covs);	


	
}


int create_id_list(const char * phenotype_filename, const char ** trts, const char ** covs, int n_trts, int n_covs){

	int * trait_index = malloc(sizeof(int)*n_trts);
	

		
	int * covariate_index = malloc(sizeof(int)*n_covs);

	printf("%s \n", phenotype_filename);

	FILE * phenotype_file = fopen(phenotype_filename, "rb");
	
	fseeko(phenotype_file, 0, SEEK_END);

	long fsize = ftello(phenotype_file);
	fseeko(phenotype_file, 0, SEEK_SET); 

	char *read_buffer = malloc(fsize );
	int bytes_read =fread(read_buffer, 1, fsize, phenotype_file);

	fclose(phenotype_file);

	int * n_lines = malloc(sizeof(int));
	*n_lines = 0; 
	char * buffer_ptr  = clean_cstring(read_buffer, n_lines);
	printf("%i\n", *n_lines);

	int * trait_indices = malloc(sizeof(int)*n_trts);
	
	int i = 0;
	while(i < n_trts) trait_indices[i++] = -1;

	
	
	
	int * cov_indices;
	if(n_covs > 0){
		cov_indices = malloc(sizeof(int)*n_covs);
		i = 0;
		while(i < n_covs) cov_indices[i++] = -1;
	}

	int id_index = -1;
	
	int found_covs = 0;
	int found_traits = 0;
	int id_found  = 0;
	char * 	line = NextArg(buffer_ptr, '\n', &buffer_ptr);
	char * token;// = NextArg(line, ',', &line);
	int found;
	int index = 0;
	int field_count = 0;
	printf("Checking pheno fields\n");

	do{
		token = NextArg(line, ',', &line);
		
		if(!id_found){
			if(!StringCmp(token, "ID", 0)){
				id_index = index;
				id_found = 1;
				index++;
		//		token = NextArg(line, ',', &line);
				continue;
			}
		}
		
		if(found_traits != n_trts){
			found = search_indices(token, trts, trait_indices, index, found_traits, n_trts);
			if(found){
				found_traits++;
				index++;
			//	token = NextArg(line, ',', &line);
				continue;
			}
			
		}
		
		if(found_covs != n_covs){
			found = search_indices(token, covs, cov_indices, index, found_covs, n_covs);
			if(found){
				found_covs++;
				index++;
			//	token = NextArg(line, ',', &line);
				continue;
			}
			
		}
		
		index++;
	}while(line != NULL  && (!id_found || found_traits != n_trts || found_covs != n_covs));
	
	if(found_traits != n_trts||found_covs != n_covs|| !id_found){
		i = 0;
		while(i < n_trts){
			
			if(trait_indices[i] == -1){
				printf("%s %i\n", trts[i], i);
			}
			i++;
		}
		
		
		
		printf("%i %i", found_traits, id_found);
		printf("Failed to find a trait, covariate, or id field when parsing fields of currently loaded phenotype file\n");
		return 1;
	}
	
	
	
//	double * trait_values = (double*)malloc(sizeof(double)*n_trts);
//	double * cov_values;
/*	if(n_covs > 0){
		cov_values = (double*)malloc(sizeof(double)*n_covs);
	}*/
	char * id;
	int array_index;
	int trt_array_index = 0;
	int cov_array_index = 0;
	int skip_row;
	int length_count = 0;
	printf("Reading entire phenotype file\n");

		
/*		
	pthread_t * threads = (pthread_t *)malloc(sizeof(pthread_t)*n_threads);
	struct phenotype_data * line_data = (struct phenotype_data*)malloc(sizeof(struct phenotype_data)*n_threads);
	int line_idx = 1;
	int pthread_index = 0;
	while(line_idx < *n_lines){
		pthread_index = 0;
		while(buffer_ptr != NULL && pthread_index < n_threads && line_idx < *n_lines){
		
			token = NextArg(buffer_ptr, '\n', &buffer_ptr);	
			
			line_data[pthread_index].line = strdup(token);
			free(token);
			line_data[pthread_index].trait_indices = trait_indices;
			line_data[pthread_index].cov_indices = cov_indices;
			line_data[pthread_index].n_traits = n_trts;
			line_data[pthread_index].n_covs = n_covs;
			line_data[pthread_index].id_index = id_index;
			
			
			pthread_create(&threads[pthread_index], NULL, read_phenotype_line, &line_data[pthread_index]);
			line_idx++;
			pthread_index++;
		}
		
		
		if(line_idx != *n_lines){
			pthread_index = 0;
			while(pthread_index < n_threads) pthread_join(threads[pthread_index++], NULL);
		}else{
			int remaining_threads = (*n_lines - 1) % n_threads;
			pthread_index = 0;
			while(pthread_index < remaining_threads) pthread_join(threads[pthread_index++], NULL);
		}
		
	}*/
	free(threads);
	free(line_data);
	free(read_buffer);
	

	while(line != NULL && StringCmp(line, 0, 0)){

		index = 0;
	//	token = NextArg(line, ',', &line);
		found_traits = 0;
		found_covs = 0;
		id_found = 0;
		skip_row = 0;
		trt_array_index = 0;
		cov_array_index = 0;
		do{
			token = NextArg(line, ',', &line);	
			if(trait_indices[trt_array_index] == index){
				if(token == NULL){
					skip_row = 1;
					break;
				}
				trait_values[trt_array_index] = strtod(token, NULL);
				index++;
				found_traits++;
				trt_array_index++;
				if(trt_array_index == n_trts)
					trt_array_index--;
			//	token = NextArg(line, ',', &line);
				continue;
			}
		//	array_index  = check_indices(index, cov_indices, found_covs, n_covs);
			
			if(n_covs != 0){
				if(cov_indices[cov_array_index] == index){
					if(token == NULL){
						skip_row = 1;
						break;
					}
				
					if(!StringCmp(token, "M", 0)){
						cov_values[array_index] = 1.0;
					}else if(!StringCmp(token, "F", 0)){
						cov_values[array_index]= 2.0;
					}else{
						cov_values[array_index] = strtod(token, NULL);
					}

					index++;
					found_covs++;
					cov_array_index++;
					
					if(cov_array_index == n_covs)
						cov_array_index--;					
		//		token = NextArg(line, ',', &line);
					continue;
				}	
			}		
			
			if(index == id_index){
				id = (char*)malloc(strlen(token) + 1);
				strcpy(id, token);
			//	printf("%s \n", token);
				index++;
				id_found = 1;
		//		token = NextArg(line, ',', &line);
				continue;
			}

			
			index++;
			
					
			
		}while(line != NULL && (found_traits != n_trts||found_covs != n_covs|| !id_found));
	//	printf("%s fus\n", id );
		if(skip_row == 0 && id != NULL)
			append_phenotype_containers(id, trait_values, cov_values, n_trts, n_covs);
			
		line = NextArg(buffer_ptr,'\n', &buffer_ptr);
		
	}




	FILE * pedindex_file = fopen("pedindex.csv", "rb");

	fseeko(pedindex_file, 0, SEEK_END);

	fsize = ftello(pedindex_file);

	fseeko(pedindex_file, 0, SEEK_SET);  

	char * pedindex_buffer = malloc(fsize );	
	
	
	bytes_read =fread(pedindex_buffer, 1, fsize, pedindex_file);
	fclose(pedindex_file);	
	printf("cleaning pedindex file\n");
	*n_lines = 0;
	buffer_ptr = clean_cstring(pedindex_buffer, n_lines);
	line = NextArg(buffer_ptr,'\n', &buffer_ptr);
	id_found  = 0;
	index = 0;
//	token = NextArg(line, ',', &line);
	int ibdid_index = -1;
	field_count = 0;
	printf("searching pedindex fields\n");

	do{
		token = NextArg(line, ',', &line);		
		
		if(!id_found){
			if(!StringCmp(token, "ID", 0)){
				id_index = index;
				id_found = 1;
				index++;
			//	token = NextArg(line, ',', &line);
				continue;
			}			
		}
		
		index++;
	}while(line != NULL  && (!id_found));
	line = NextArg(buffer_ptr,'\n', &buffer_ptr);	
	int ibdid = 1;
	int n_threads;
	if(*n_lines - 1 < N_PTHREADS)
		n_threads = *n_lines - 1;
	else
		n_threads = N_PTHREADS;
	while(line != NULL && StringCmp(line, 0, 0)){
		index = 0;
		
	
	//	if(token == NULL)
	//		break;


		id_found = 0;
		do{
			
			token = NextArg (line, ',', &line);

			if(index == id_index && token != NULL){
				id = (char*)malloc(strlen(token) + 1);
				id_found = 1;
				strcpy(id, token);
			}

			index++;

		}while(line != NULL && (!id_found) && StringCmp(line, 0, 0));
		//printf("%s \n", id);
		if(id != NULL)
			append_ibdids_and_ids(id, ibdid);
		id = NULL;
		ibdid++;
	
		line = NextArg(buffer_ptr,'\n', &buffer_ptr);		
	}
	free(pedindex_buffer);
	
	return 0;
}




