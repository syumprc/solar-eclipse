#include<plinkio.h>
#include<stdio.h>
#include "safelib.h"

#define SQUARED(value) ((value)*(value))


void write_matrix_to_file(char * output_filename, struct pio_file_t * input_file,  double * GRM, size_t n_subjects);

void create_relationship_matrix(struct pio_file_t * input_file, const char * basename, double * GRM, int per_chromosome, size_t n_subjects , size_t n_snps){

	
	// * int_buffer = new unsigned int[n_subjects];




	int * n_snps_matrix= Calloc(n_subjects*n_subjects, sizeof(int));
	snp_t * buffer = Malloc(n_subjects*sizeof(snp_t));
	//memset(n_snps_matrix, 0, sizeof(int)*n_subjects*n_subjects);


	struct pio_locus_t * locus;
	pio_status_t status;


	double freq;


	unsigned int homo_count;
	unsigned int heter_count;
	unsigned int subject_count;
	locus = bim_get_locus(&input_file->bim_file, 0);

	size_t row, col, current_subject;

	unsigned char  chromosome = locus->chromosome;
	
	
	int snp_count=0;
	
	for(size_t snp = 0 ; snp < n_snps; ++snp){


		
		status = pio_next_row(input_file,buffer);
		locus = bim_get_locus(&input_file->bim_file, snp);
		if(per_chromosome == 1 && locus->chromosome != chromosome) {
			
			for(size_t row = 0 ; row < n_subjects ; ++row) \
					for(size_t col = row;  col < n_subjects; ++col) GRM[row*n_subjects + col]/= (n_snps_matrix[row*n_subjects + col] + snp_count);
					
			char str_buffer[200];
			sprintf(str_buffer, "%s.chr%i.csv", basename, chromosome);
			write_matrix_to_file(str_buffer, input_file, GRM, n_subjects);					
			
			memset(n_snps_matrix, 0, sizeof(int)*n_subjects*n_subjects); 
			memset(GRM, 0, sizeof(double)*n_subjects*n_subjects); 
			snp_count = 0;
			chromosome = locus->chromosome;
			
		}
		
	
		if(status != PIO_OK){
			continue;
		}

		
		if(!StringCmp(locus->allele1, "0", 0) || !StringCmp(locus->allele2, "0", 0)){
			continue;
		}
		snp_count++;
		freq = 0.0;


		subject_count = n_subjects;
		homo_count = 0;
		heter_count = 0;







		snp_t * it_buffer = &buffer[0];

		for(size_t subject = 0 ; subject < n_subjects ; ++subject){

			switch(*it_buffer++){
			case 1:
				heter_count++;
				break;
			case 2:
				homo_count++;
				break;
			case 3:
				subject_count--;

				for(size_t col = subject; col < n_subjects ; ++col) n_snps_matrix[subject*n_subjects + col]--;

				for(size_t row = 0; row < subject; ++row) n_snps_matrix[row*n_subjects + subject]--;




				break;

			}




		}

		freq = (double)(homo_count*2.0 + heter_count)/(double)(2.0*subject_count);



		snp_t * row_value = &buffer[0];
		for(size_t row = 0; row < n_subjects; ++row, ++row_value){
			if(*row_value == 3) continue;
			GRM[row*n_subjects + row] += SQUARED((*row_value - 2.0*freq))/(freq*2.0*(1.0 - freq));
			snp_t * col_value = &buffer[n_subjects - 1];
			double factor = (*row_value - 2.0*freq)/(freq*(1.0 - freq)*2.0);
			for(size_t col = n_subjects - 1 ; col > row ; --col, --col_value) GRM[row*n_subjects + col] += factor*((*col_value != 3 ? *col_value : 2.0*freq)  - \
							2.0*freq);
		}

		

	}
	free(buffer);
	
	for(size_t row = 0 ; row < n_subjects ; ++row) \
	for(size_t col = row;  col < n_subjects; ++col) GRM[row*n_subjects + col]/= (n_snps_matrix[row*n_subjects + col] + snp_count);
	if(per_chromosome == 0){
		char str_buffer[200];
		sprintf(str_buffer, "%s.csv", basename);
		write_matrix_to_file(str_buffer, input_file, GRM, n_subjects);
	}else{
		char str_buffer[200];
		sprintf(str_buffer, "%s.chr%i.csv", basename, chromosome);
		write_matrix_to_file(str_buffer, input_file, GRM, n_subjects);
	}

	
	free(n_snps_matrix);

}




void create_relationship_matrix_with_header(struct pio_file_t * input_file, double * GRM, const char * basename,
		const char * header_list[], int n_headers, int per_chromosome, size_t n_subjects , size_t n_snps){

	
	// * int_buffer = new unsigned int[n_subjects];




	int * n_snps_matrix= Calloc(n_subjects*n_subjects, sizeof(int));
	snp_t * buffer = Malloc(n_subjects*sizeof(snp_t));
	//memset(n_snps_matrix, 0, sizeof(int)*n_subjects*n_subjects);


	struct pio_locus_t * locus;
	pio_status_t status;


	double freq;


	unsigned int homo_count;	 
	unsigned int heter_count;
	unsigned int subject_count;
	locus = bim_get_locus(&input_file->bim_file, 0);

	size_t row, col, current_subject;

	unsigned char  chromosome = locus->chromosome;
	
	
	int snp_count=0;
	for(size_t snp = 0 ; snp < n_snps; ++snp){



		status = pio_next_row(input_file,buffer);
		locus = bim_get_locus(&input_file->bim_file, snp);

		
		const char * name =  (const char *)locus->name;
		int is_header = 0;
		for(int i = 0 ; i < n_headers; i++){
			if (!StringCmp(header_list[i], name, 0)){
				is_header = 1;
				break;
			}
		}
				

		if(is_header == 0){
			continue;
		}
		if(per_chromosome == 1 && locus->chromosome != chromosome) {
			
			for(size_t row = 0 ; row < n_subjects ; ++row) \
					for(size_t col = row;  col < n_subjects; ++col) GRM[row*n_subjects + col]/= (n_snps_matrix[row*n_subjects + col] + snp_count);
					
			char str_buffer[200];
			sprintf(str_buffer, "%s.chr%i.csv", basename, chromosome);
			write_matrix_to_file(str_buffer, input_file, GRM, n_subjects);					
			
			memset(n_snps_matrix, 0, sizeof(int)*n_subjects*n_subjects); 
			memset(GRM, 0, sizeof(double)*n_subjects*n_subjects); 
			snp_count = 0;
			
		}				
		

		if(status != PIO_OK){
			continue;
		}

		
		if(StringCmp(locus->allele1, "0",  0) || StringCmp(locus->allele2, "0", 0)){
			continue;
		}
		
		snp_count++;

		freq = 0.0;


		subject_count = n_subjects;
		homo_count = 0;
		heter_count = 0;







		snp_t *it_buffer = &buffer[0];

		for(size_t subject = 0 ; subject < n_subjects ; ++subject){

			switch(*it_buffer++){
			case 1:
				heter_count++;
				break;
			case 2:
				homo_count++;
				break;
			case 3:
				subject_count--;

				for(size_t col = subject; col < n_subjects ; ++col) n_snps_matrix[subject*n_subjects + col]--;

				for(size_t row = 0; row < subject; ++row) n_snps_matrix[row*n_subjects + subject]--;




				break;

			}




		}

		freq = (double)(homo_count*2.0 + heter_count)/(double)(2.0*subject_count);



		snp_t * row_value = &buffer[0];
		for(size_t row = 0; row < n_subjects; ++row, ++row_value){
			if(*row_value == 3) continue;
			GRM[row*n_subjects + row] += SQUARED((*row_value - 2.0*freq))/(freq*2.0*(1.0 - freq));
			snp_t * col_value = &buffer[n_subjects - 1];
			double factor = (*row_value - 2.0*freq)/(freq*(1.0 - freq)*2.0);
			for(size_t col = n_subjects - 1 ; col > row ; --col, --col_value) GRM[row*n_subjects + col] += factor*((*col_value != 3 ? *col_value : 2.0*freq)  - 2.0*freq);
		}



	}
	free(buffer);
	
	for(size_t row = 0 ; row < n_subjects ; ++row) \
	for(size_t col = row;  col < n_subjects; ++col) GRM[row*n_subjects + col]/= (n_snps_matrix[row*n_subjects + col] + snp_count);
	if(per_chromosome == 0){
		char str_buffer[200];
		sprintf(str_buffer, "%s.csv", basename);
		write_matrix_to_file(str_buffer, input_file, GRM, n_subjects);
	}else{
		char str_buffer[200];
		sprintf(str_buffer, "%s.chr%i.csv", basename, chromosome);
		write_matrix_to_file(str_buffer, input_file, GRM, n_subjects);
	}

	
	free(n_snps_matrix);

}





#undef SQUARED

