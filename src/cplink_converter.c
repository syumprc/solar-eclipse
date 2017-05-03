#include <stdio.h>
#include <math.h>
#include "plinkio.h"
#include "safelib.h"

static int count_chromosomes(struct pio_file_t * input, size_t n_snps){
	int n_chromo = 1;
	struct pio_locus_t * locus = bim_get_locus(&input->bim_file, 0);
	unsigned char chromosome = locus->chromosome;
	for(size_t snp = 1; snp < n_snps; snp++){
		locus =  bim_get_locus(&input->bim_file, snp);
		if(locus->chromosome != chromosome) {
			printf("here\n");
			n_chromo++;
			chromosome = locus->chromosome;
		}
	}
	return n_chromo;
}

void read_and_write_to_csv(struct pio_file_t  * input, const char* output_basename,
			  size_t col_size, size_t max_per_file, int bin_format, int per_chromo){

  unsigned int row_size = input->bed_file.header.num_samples;
  int num_files;

  if(per_chromo){
	  num_files = count_chromosomes(input, col_size);
	  max_per_file = col_size;
  }else if (max_per_file == 0){
      max_per_file = col_size;
      num_files = 1;
  }else{
      num_files = ceil((double)col_size/(double)max_per_file);
  }



  snp_t buffer[row_size];
  struct pio_locus_t * locus;
  FILE *  oss[num_files];

  FILE * oss_header[num_files];
  printf("files %i \n", num_files);
  if(per_chromo){
	 int chromo = 0;  
      for(unsigned int i = 0 ; i < num_files; i++){
		char output_filename[200];
      // = output_basename + ".csv";
		sprintf(output_filename, "%s.chr%i.csv",output_basename, chromo);
		oss[i] = fopen(output_filename, "w");
          if(bin_format){
			  char output_filename_header[200];// = output_basename + ".header.csv";
			  sprintf(output_filename_header, "%s.chr%i.header.csv",output_basename, chromo);              
              oss_header[i] = fopen(output_filename_header, "w");
          }
		
		chromo++;
      }	  

  }	 else if (num_files == 1){
      char output_filename[200];
      // = output_basename + ".csv";
      sprintf(output_filename, "%s.csv",output_basename);
      oss[0] = fopen(output_filename, "w");//.open(output_filename.c_str(), ios::out);
      if(bin_format){
          char output_filename_header[200];// = output_basename + ".header.csv";
          sprintf(output_filename_header, "%s.header.csv",output_basename);
          oss_header[0] = fopen(output_filename_header, "w");//.open(output_filename_header.c_str(), ios::out);
      }
  }else{

      for(unsigned int i = 0 ; i < num_files; i++){
		char output_filename[200];
      // = output_basename + ".csv";
		sprintf(output_filename, "%s_%i.csv",output_basename, i);
		oss[i] = fopen(output_filename, "w");
          if(bin_format){
			  char output_filename_header[200];// = output_basename + ".header.csv";
			  sprintf(output_filename_header, "%s_%i.header.csv",output_basename, i);              
              oss_header[i] = fopen(output_filename_header, "w");
          }
		

      }


  }




  fprintf(oss[0],"id,fid,sex,phenotype,");//oss[0] << "id,fid,sex,phenotype,";
  if(bin_format){
	fprintf(oss_header[0],"SNP,0,1,2\n");
  }
  int file = 0;
  short bits;
  locus = bim_get_locus(&input->bim_file, 0);
  unsigned char current_chromo = locus->chromosome;
  unsigned char new_chromo = current_chromo;
  int changed_chromo = 0;

  for(unsigned int i = 0 ; i < col_size ; i++){
	  
      if((max_per_file*(file + 1)) == i || (new_chromo != current_chromo && per_chromo)){
		
	  file++;
	  printf("%i %u \n", file, i);
	  current_chromo = new_chromo;
	  fprintf(oss[file],"id,fid,sex,phenotype,");
          if(bin_format) {
			fprintf(oss_header[file], "SNP,0,1,2\n");
          }
      }	  
	  if(i != col_size -1 && per_chromo){
		locus = bim_get_locus(&input->bim_file, i + 1);		  
		   new_chromo = locus->chromosome;
	   }	  
		locus = bim_get_locus(&input->bim_file, i);		
		
	  	  fprintf(oss[file], "snp_%s", locus->name);
          if(bin_format)
			fprintf(oss_header[file],"snp_%s,%s%s,%s%s,%s%s\n", locus->name , locus->allele1 , locus->allele1 \
              , locus->allele1, locus->allele2 ,  locus->allele2 , locus->allele2);		  

			if(i == ((file + 1)*(max_per_file) - 1) || (i == col_size - 1)|| (new_chromo != current_chromo && per_chromo)){
    	       fprintf(oss[file],"\n");
    	   }else{
    	       fprintf(oss[file],",");
    	   }      
  }

  struct pio_sample_t * sample;
  for(size_t row = 0 ; row < row_size ; row++){


      sample = fam_get_sample(&input->fam_file, row);
      file = 0;
      fprintf(oss[file], "%s,%s,%i,%f,", sample->iid, sample->fid, sample->sex, sample->phenotype);
	  locus = bim_get_locus(&input->bim_file, 0);     
      current_chromo = locus->chromosome;
      new_chromo = current_chromo;
      for(size_t col = 0 ; col < col_size ; col++){
		   pio_status_t status = pio_next_row(input, buffer);

			locus = bim_get_locus(&input->bim_file, col); 
			
		   if((max_per_file*(file + 1)) == col || (new_chromo != current_chromo && per_chromo)){
			file++;
			sample = fam_get_sample(&input->fam_file, row);
			current_chromo = new_chromo;
			fprintf(oss[file], "%s,%s,%i,%f,", sample->iid, sample->fid, sample->sex, sample->phenotype);
		   }
		   bits = buffer[row];
	
	   
		//   locus = bim_get_locus(&input->bim_file, col);
    	   if((!StringCmp(locus->allele1, "0", 0)
    	       || !StringCmp(locus->allele2, "0", 0))  && bits == 1){
    	       bits = 3;
    	   }else if((!StringCmp(locus->allele1, "0", 0) )
    		  && bits == 0){
    		  bits = 3;
    	   }else if((!StringCmp(locus->allele2, "0", 0))
 	     && bits == 2){
    	       bits = 3;
 	   }


    	   if(bits == 0){

    	       if(bin_format){
    		   fprintf(oss[file], "%i", bits);
    	       }else{
				  fprintf(oss[file], "%s%s", locus->allele1,locus->allele1);
    	       }
    	   }else if (bits == 1){
    	       if(bin_format){
    		   fprintf(oss[file], "%i", bits);
    	       }else{
			   fprintf(oss[file], "%s%s", locus->allele1,locus->allele2);
  
    	       }
    	   }else if (bits == 2){
    	       if(bin_format){
    		   fprintf(oss[file], "%i", bits);
    	       }else{
				 fprintf(oss[file], "%s%s", locus->allele2,locus->allele2);	   
    	       }

			}
		
	  if(col != col_size -1 && per_chromo){
		locus = bim_get_locus(&input->bim_file, col + 1);		  
		   new_chromo = locus->chromosome;
	   }				
		
    	   if(col == ((file + 1)*(max_per_file) - 1) || (col == col_size - 1)|| (new_chromo != current_chromo && per_chromo)){
    	       fprintf(oss[file],"\n");
    	   }else{
    	       fprintf(oss[file],",");
    	   }
      }
      pio_reset_row(input);
  }

  for(int f = 0 ; f < num_files ; f++){
	  fclose(oss[f]);
      if(bin_format)
          fclose(oss_header[f]);
  }

  pio_close(input);

  

}
