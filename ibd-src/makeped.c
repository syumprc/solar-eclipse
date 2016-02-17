
/*****************************************************************************/
/*                                                                           */
/*  Program:  MAKEPED.C                                                      */
/*                                                                           */
/*                                                                           */
/*  Used to convert pedigree files missing sib pointers to pedigree files    */
/*  with sib pointers.                                                       */
/*                                                                           */
/*  Revision history:                                                        */
/*                                                                           */
/*  6-23-88  The ind_lookup routine stopped looking for an id too soon.      */
/*           This showed up when strings where used as id's and the pedigree */
/*           was entered in the reverse order from the way you would draw    */
/*           a pedigree, that is, an idividual drawn at bottom of a pedigree */
/*           was listed as the first individual in the data file.            */
/*                                                                           */
/* 10-05-88  Character spacing of output when loops created an new individual*/
/*           has been corrected.                                             */
/*                                                                           */
/* 10-09-88  Loops have been corrected so that offspring refer to the newly  */
/*           created individual as a parent instead of the original parent.  */
/*                                                                           */
/* 24-12-88  Loops have been corrected so that nextoffspring points are      */
/*           kept for the original individual instead of new individual.     */
/*                                                                           */
/*****************************************************************************/

#define mversion 2.1

#include <stdio.h>
#include <ctype.h>
#ifdef MICROSOFT_C
  #include <malloc.h>
#elif TURBO_C
  #include <alloc.h>
#else
  #include <stdlib.h>
  #include <string.h>
#endif


#define	FALSE			0
#define	TRUE			1

typedef unsigned char   u_byte;
typedef          char   s_byte;
typedef unsigned short  u_word;
typedef          short  s_word;
typedef unsigned long   uns_long;
typedef          long   s_long;
typedef unsigned int    u_intg;
typedef          int    s_intg;
typedef          float  s_real;
typedef          double d_real;

        /* User defined constants */

#define maxlocus     1        /* maximum number of loci                      */
#define maxall       500      /* maximum number of alleles at one locus      */
#define maxallchars  200      /* maximum chars used in phenotypic data       */
#define maxind       50001    /* maximum number of individuals               */
#define maxped       5000     /* maximum number of pedigrees                 */
#define setindex     31       /* (setindex+1)*8 = num of elements in a set   */
#define maxname      20       /* maximum number of chars in ped or person id */
#define maxperped    5000     /* maximum number of persons per pedigree      */

        /* Quantitative trait */

#define maxtrait     1        /* maximum number of quantitative factors      */
#define missval      0        /* missing value for quantitative trait        */
#define affall       2        /* allele giving affval for quantitative       */
                              /* traits in males, sexlink                    */

        /* Affection status */

#define missaff      0        /* missing value for affection status          */
#define affval       2        /* code for affected individual                */
#define maxliab      2        /* maximum number of liability classes         */ 

        /* Binary factors */

#define maxphen      36       /* maximum number of binary codes at one locus */

        /* File pointers */

FILE *pedfile;                /* pointer to pedigree information             */
FILE *pedout;                 /* pointer to new pedigree output file         */
FILE *datafile;               /* pointer to locus definition file            */

        /* Others */

#define max_filespec   80     /* max number of chars in a filename           */


typedef u_byte binset[ setindex ];  /* type used to simulate sets */
typedef binset phenarray[maxall];   /* type used to hold phenotypic data */
typedef s_real means[maxall][maxall][maxtrait];
typedef s_real covmatrix[maxtrait][maxtrait];

struct locus_values{
  s_intg    nallele;
  s_real    freq[ maxall ];
  struct    locus_values *privlocus;
  enum {affection,quantitiative,binary,null} which;
};

struct phenotype {
  u_byte phen_chars[maxallchars];
};

struct ind{
  s_byte oldped_s[maxname]; /* original pedigree read from pedfile if string */      
  s_intg oldped;          /* original pedigree read from pedfile if integer  */      
  s_byte oldid_s[maxname];/* original person id read from pedfile if string  */ 
  s_intg oldid;           /* original person id read from pedfile if integer */
  s_intg ped;             /* new pedigree id                                 */
  s_intg id;              /* new person id                                   */
  s_intg paid;            /* new paternal id                                 */
  s_intg maid;            /* new maternal id                                 */
  s_intg offid;           /* new first offspring id                          */
  s_intg npaid;           /* new next paternal offspring id                  */
  s_intg nmaid;           /* new next maternal offspring id                  */
  s_intg sex;             /* sex of person                                   */
  s_intg proband;         /* proband                                         */
  struct ind *pa;         /* pointer to fathers ind structure                */
  struct ind *ma;         /* pointer to mothers ind structure                */
  struct ind *foff;       /* pointer to first offsprings ind structure       */
  struct ind *nextpa;     /* pointer to next paternal offspring ind struct   */
  struct ind *nextma;     /* pointer to next maternal offspring ind struct   */
  s_intg generations;     /* number of generations below this individual     */
  struct phenotype *phen; /* pointer to persons phenotype structure          */
  s_intg male;            /* boolean flag for sex of individual              */
  u_intg is_parent;       /* flag used to find people with no family         */
} ;

        /* arrays of individuals */

struct ind *person[ maxind ];  /* array of pointers to ind structures        */
struct ind *proband[ maxped ]; /* array of pointers to pedigree probands     */

        /* others */

s_intg nsystem;           /* number of systems                               */
s_intg risksys;           /* risk system                                     */
s_intg mutsys;            /* mutation system                                 */
s_intg autosomal;         /* autosomal data                                  */
s_intg mutlocus;          /* mutation locus                                  */
s_real mutrate;           /* mutation rate                                   */
s_intg hapfreq;
s_intg nupriv;
s_intg lastpriv;
s_intg nneeded;
s_intg whichlocus;
s_intg whichtype;         /* type of a locus at read time, binary, etc.      */
s_intg nuped;
s_intg interfer;
s_intg disequi     = FALSE;
s_intg sexlink     = FALSE;
s_intg mutup       = FALSE;
s_intg risk        = FALSE;
s_intg sexdif      = FALSE;
s_intg allnsystems = FALSE;
s_intg ped_integers= FALSE; /* Pedigree identifiers are integers if true.   */
s_intg ind_integers= FALSE; /* Individual identifiers are integers if true. */
s_intg next_id;             /* Person id's used if originals were strings.  */
s_intg family_count;        /* Number of people in family being processed.  */ 
s_intg family_offset;       /* Number of people in previous families.       */
s_intg nuperson;            /* Number of people read so far in a pedigree.  */
s_intg totperson;           /* Total num persons in pedigree file.          */
s_intg probands[maxped];    /* Indexes of people assigned probands.         */
s_intg loops[maxped];       /* Indexes of people assigned loops.            */
u_intg found_error;         /* Error flag for sex and phenotype check.      */
u_intg force_strings;       /* Force individual id's to be strings.         */
u_byte pifile[max_filespec];/* Pedigree input file.                         */
u_byte pofile[max_filespec];/* Revised pedigree output file.                */
s_intg biggest_p_id;        /* Largest pedigree id found or created.        */
s_intg biggest_i_id;        /* Largest individual id found or created.      */

 /* The following are set to true if the pedigree has allready been cleared */
 /* of probands set by auto_probands.  You only want to clear probands from */
 /* a pedigree once so that a user can set more than one proband/pedigree   */
 /* if he wants to.                                                         */

u_byte cleared[maxped] = {FALSE};


/****************************************************************************/
/*                                                                          */
/*			   read_pedigree   				    */
/*                                                                          */
/* Prompts for the identifer of a pedigree and verifies that it's an integer*/
/* if id's are supposed to be integers and a string if they are not.   If   */
/* the format is incorrect then the prompt is repeated.                     */
/*                                                                          */
/****************************************************************************/

void read_pedigree(pedigree_i, pedigree_s)
 s_intg *pedigree_i;    /* integer id's */
 s_byte pedigree_s[];   /* string id's  */
{

  s_intg good_format;
  s_intg i;

  do {
  good_format = TRUE;
  fprintf(stdout,"\n\tPedigree   -> ");
  fscanf(stdin,"%s",pedigree_s);

  if (pedigree_s[0] != '0')  /* Equals 0 when done entering pedigrees. */

    if (ped_integers) {
      i = 0;
      while((pedigree_s[i] != '\0') && (good_format))
        if (isalpha(pedigree_s[i++])) {
          good_format = FALSE;
          fprintf(stdout,"\tInteger identifier expected...\n");
        }
    }
    else {
      if ((!ped_integers) && (isdigit(pedigree_s[0]))) {
        fprintf(stdout,"\tString identifier expected...\n");
        good_format = FALSE;
      }
    }
 } while(!good_format);

 /* Convert to integer if needed. */

 if (ped_integers) sscanf(pedigree_s,"%d",pedigree_i);
}





/****************************************************************************/
/*                                                                          */
/*			   read_person   				    */
/*                                                                          */
/* Prompts for the identifer of a person and verifies that it is integer    */
/* if id's are supposed to be integers and a string if they are not.   If   */
/* the format is incorrect then the prompt is repeated.                     */
/*                                                                          */
/****************************************************************************/

void read_person(person_i, person_s)
 s_intg *person_i;   /* integer id's */
 s_byte person_s[];   /* string id's  */
{

  s_intg good_format;
  s_intg i;

  do {
  good_format = TRUE;
  fprintf(stdout,"\tPerson     -> ");
  fscanf(stdin,"%s",person_s);
  if (ind_integers) {
    i = 0;
    while((person_s[i] != '\0') && (good_format))
      if (isalpha(person_s[i++])) {
        good_format = FALSE;
        fprintf(stdout,"\tInteger identifier expected...\n");
      }
  }
  else {
    if ((!ind_integers) && (!force_strings) && (isdigit(person_s[0]))) {
      fprintf(stdout,"\tString identifier expected...\n");
      good_format = FALSE;
    }
  }  
 } while(!good_format);

  /* Convert to integer if needed. */

 if (ind_integers) sscanf(person_s,"%d",person_i);
}




/****************************************************************************/
/*                                                                          */
/*                            ind_lookup                                    */
/*                                                                          */
/****************************************************************************/

s_intg ind_lookup(name,sequence)
  s_byte name[];
  s_intg sequence;  /* the number of people in privious pedigrees */
{

/* Search through all people in this pedigree that have been read in so far */
/* and find the one whose original id matches that in the name parameter.   */
/* If one is found then return his id else return 0.                        */

s_intg i;

i = 1;
while (i <= nuperson) {
  if (person[sequence+i] == NULL) return(0);
  else
   if (!strcmp_i(person[sequence + i]->oldid_s,name) ) /* if they match */
    return(person[sequence + i]->id);
  else
    i++;
}
return(0);
}


/*****************************************************************************/
/*                                                                           */
/*                             getind                                        */
/*                                                                           */
/*****************************************************************************/

void getind(id,sequence,newped_i,newped_s,nuped)
     s_intg *id;
     s_intg sequence;
     s_intg newped_i;
     s_byte newped_s[];
     s_intg nuped;

{
  s_intg temp_i;               /* Holds id if id's are integers. */
  s_byte temp_s[maxname];      /* Holds id if id's are strings.  */
  s_intg found_id;

  /* Read the persons id and convert it to an initial index. */

  if (ind_integers) {
    fscanf(pedfile,"%d",id);

    temp_i = *id;

  }
  else {
    fscanf(pedfile,"%s",temp_s);
    if (temp_s[0]=='0') *id = 0;
    else {
      found_id = ind_lookup(temp_s,sequence);
      if (found_id) *id = found_id;
      else *id = next_id;
    }
  }

    if ((ind_integers) && (temp_i > maxind)) {
     fprintf(stderr,"ERROR: Ped: %d  Per: %d - ",newped_i,temp_i);
     fprintf(stderr,"maximum id of %d exceeded.\n",maxind - 1);
     exit(1);
    }


  /* If the initial index in not zero then compute the final index, */
  /* and allocate memory for their record if not already done, and  */
  /* initialize their record fields.                                */

  if  (*id != 0) {
      *id += sequence;
      if (person[*id] == NULL) {
	person[*id] = (struct ind *) calloc ( 1, sizeof( struct ind ) );

	if ( person[*id] == NULL ) {
	  fprintf(stderr,"ERROR: Cannot allocate memory.\n");
          exit(1);
	}

      if (ped_integers)
        person[*id]->oldped = newped_i;
      else
        strcpy(person[*id]->oldped_s,newped_s);

      if (ind_integers)
	person[*id]->oldid      = temp_i;
      else
        strcpy(person[*id]->oldid_s,temp_s);

      if (ind_integers){
	person[*id]->id       = temp_i;
        if (temp_i > biggest_i_id) biggest_i_id = temp_i;
      }
      else {
	person[*id]->id       = next_id;
        if (next_id > biggest_i_id) biggest_i_id = next_id;
        next_id++;
      }


      if (ped_integers){
        if (newped_i > biggest_p_id) biggest_p_id = newped_i;
      }
      else
        if (nuped > biggest_p_id) biggest_p_id = nuped;

      person[*id]->ped        = nuped;
      person[*id]->paid       = 0;
      person[*id]->maid       = 0;
      person[*id]->offid      = 0;
      person[*id]->npaid      = 0;
      person[*id]->nmaid      = 0;
      person[*id]->pa         = NULL;
      person[*id]->ma         = NULL;
      person[*id]->foff       = NULL;
      person[*id]->nextpa     = NULL;
      person[*id]->nextma     = NULL;
      nuperson++;

      }
    }
}


/*****************************************************************************/
/*                                                                           */
/*                             getphenotype                                  */
/*                                                                           */
/* This version of getphenotype simply reads in the phenotypic data as a     */
/* string and assigns it to the persons phen->phen_chars field.  This should */
/* be used when no interpretation of phenotypic data is intended.            */
/*                                                                           */
/*****************************************************************************/

void getphenotype(id)
     s_intg *id;

{
  s_intg i;
  s_byte c;

  person[*id]->phen =
    (struct phenotype *) calloc ( 1, sizeof( struct phenotype ) );
  if ( person[*id]->phen == NULL ) {
      fprintf(stderr,"ERROR: Cannot allocate memory.\n");
     exit(1);
  }

  i = 0;
  c = getc(pedfile);
  while(( c != '\n') && (!feof(pedfile))) {
    person[*id]->phen->phen_chars[i++] = c;
    c = getc(pedfile);
  }

  person[*id]->phen->phen_chars[i] = '\0';

}



/*****************************************************************************/
/*                                                                           */
/*                             readped                                       */
/*                                                                           */
/*****************************************************************************/

void readped()
{
  s_intg i;
  s_intg j;
  s_intg newid;
  s_intg sex;
  s_intg profield;
  s_intg sequence;                /* number of people read in pedigree */
  s_intg thisone;
  s_intg thisped_i;               /* integer pedigree id */
  s_byte thisped_s[maxname];      /* string pedigree id  */
  s_intg newped_i;                /* integer pedigree id */
  s_byte newped_s[maxname];       /* string pedigree id  */
  s_byte c;

  /* initialize some stuff */

  totperson      = 0;
  sequence       = 0;
  nuperson       = 0;
  nuped          = 1;
  next_id        = 1;
  proband[nuped] = NULL;
  for( i=0; i<maxind; i++) person[i] = NULL;

    /* Read first pedigree number*/

  fscanf(pedfile,"%s",newped_s);

    /* Test first character of first pedigree id to  */
    /* see if strings are being used as indentfiers. */

   if (isalpha(newped_s[0])) {
     ped_integers = FALSE;
     strcpy(thisped_s,newped_s);
   }  
   else {
     ped_integers = TRUE;
     sscanf(newped_s,"%d",&newped_i);  /* convert to integer */
     thisped_i = newped_i;
   }

  /* If force_strings is true then force individual id's to be strings. */
  /* Otherwise test the first character and determin wheather they are  */
  /* strings or integers.                                               */

  if (force_strings) ind_integers = FALSE;
  else {
    c = getc(pedfile);
    while  ( c == ' ' ||  c == '\t') c = getc(pedfile);
    if (isalpha(c)) ind_integers = FALSE;
    else ind_integers = TRUE;
  }

  /* Now the type is known for pedigree and individual id's, rewind */
  /* file, read past the first pedigree id, and leave the file pointer */
  /* pointing to the first individual id. */

  rewind(pedfile);
  fscanf(pedfile,"%s",newped_s);

  while (! feof(pedfile)) {

    /* Get person. */

    getind(&thisone,sequence,thisped_i,thisped_s,nuped); 

    /* Get persons father. */

    getind(&newid,sequence,thisped_i,thisped_s,nuped);
    if ((person[thisone]->pa = person[newid]) != NULL)
      person[thisone]->paid   = person[newid]->id;

    /* Get persons mother. */

    getind(&newid,sequence,thisped_i,thisped_s,nuped);
    if ((person[thisone]->ma = person[newid]) != NULL)
      person[thisone]->maid   = person[newid]->id;

    fscanf(pedfile,"%d",&sex);
    person[thisone]->sex = sex;
    if (sex == 1) person[thisone]->male = TRUE;
    else person[thisone]->male = FALSE;

    getphenotype(&thisone);

    /* Read in the next the pedigree id, if it's the start of the next */
    /* pedigree then update a few things and keep going.               */

    if (!feof(pedfile)) {
      if (ped_integers) {
        fscanf(pedfile,"%d",&newped_i);
        if (thisped_i != newped_i) {
          sequence += nuperson;
          nuperson = 0;
          thisped_i = newped_i;
          nuped++;
          next_id = 1; /* next id used if person id's were strings */
        }
      }
      else {
        fscanf(pedfile,"%s",newped_s);
        if (strcmp_i(thisped_s,newped_s) != 0) {
        sequence += nuperson;
        nuperson = 0;
        strcpy(thisped_s,newped_s);
        nuped++;
        next_id = 1;  /* next id used if person id's were strings */
        }
      }
    }
  }
  totperson = nuperson + sequence;
}


void pointers()
{
  struct ind *q;
  s_intg      i;
  s_intg      count;      /* number of people in current pedigree       */
  s_intg      ped_count;  /* number of people in all previous pedigrees */
  s_intg      ped_id;     /* current pedigree number                    */

     /* Note: these variable are easy to confuse... */
     /*                                             */
     /*           id's    pointer's                 */
     /*           ------  ---------                 */
     /*           paid    pa                        */
     /*           maid    ma                        */
     /*           offid   foff                      */
     /*           npaid   nextpa                    */
     /*           nmaid   nextma                    */

  count      = 0;
  ped_count  = 0;
  ped_id     = 0;

  for(i=1; i<=totperson; i++){
   if (person[i] != NULL) {

     if (person[i]->ped != ped_id) {  /* if begining of next pedigree... */
      ped_id    = person[i]->ped;
      ped_count += count;
      count     = 0;
    }
     ++count;

     if ( person[i]->paid != 0 ) {
       q = person[ person[i]->paid + ped_count ]; /* get pointer to father */

     /* If the fathers first offspring is not set then set it to person i */

       if (q->offid == 0) {
	 q->offid = i - ped_count;
	 q->foff  = person[i];
       }
       else {

     /* If the fathers first offspring is set then work your way down the */
     /* next paternal chain, starting from the first offspring, until it  */
     /* runs out. Then set the next paternal of the last person in the    */
     /* chain to the original person i that you started out with.         */

	 q = person[q->offid + ped_count];
	 while (q->nextpa != NULL) q = person[q->npaid + ped_count];
         q->npaid  =  i - ped_count;
	 q->nextpa = person[i];
       }
     }

     /* Repeat the above procedure for the mothers side. */

     if ( person[i]->maid != 0 ) {
       q = person[ person[i]->maid + ped_count ]; /* get pointer to mother */
       if (q->offid == 0) {
	 q->offid =  i - ped_count;
	 q->foff  = person[i];
       }
       else {
	 q = person[q->offid + ped_count];
	 while (q->nextma != NULL) q = person[q->nmaid + ped_count];
         q->nmaid  =  i - ped_count;
	 q->nextma = person[i];

       }
     }
   }
 }
}


/*****************************************************************************/
/*                                                                           */
/*                          save_loops                                       */
/*                                                                           */
/*****************************************************************************/

void save_loops(count)
    s_intg count;
{
  s_intg i;
  s_byte response;
  s_byte loop_file[max_filespec];
  FILE   *loopf;


  fprintf(stdout,"\n\nDo you want these selections saved ");
  fprintf(stdout,"for later use?  (y/n) -> ");
  fscanf(stdin,"%1s",&response);

  if ((response == 'y') || (response == 'Y')) {
    fprintf(stdout,"\nEnter filename -> ");
    fscanf(stdin,"%s",loop_file);
    if ( (loopf = fopen(loop_file,"w")) == NULL) {
      fprintf(stderr,"ERROR: Cannot open file %s\n",loop_file);
    }
    else {
      for(i=0; i<count; i++) {

       if (ped_integers)
         fprintf(loopf,"%d ",person[ loops[i] ]->oldped);
       else
         fprintf(loopf,"%s ",person[ loops[i] ]->oldped_s);

       if (ind_integers)
         fprintf(loopf,"%d\n",person[ loops[i] ]->oldid);
       else
         fprintf(loopf,"%s\n",person[ loops[i] ]->oldid_s);

      }
      fclose(loopf);
    }
  }
 }



/****************************************************************************/
/*                                                                          */
/*                            largest_id                                    */
/*                                                                          */
/* Given the index of a person, this routine returns the largest intger id  */
/* found in that persons pedigree.                                          */
/*                                                                          */
/****************************************************************************/

s_intg largest_id(person_index)
   s_intg person_index;
{

   s_intg pedigree_number;
   s_intg largest;
   s_intg i;

   largest = person[person_index]->id;

   i = person_index -1;
   while (( i >= 1 ) && (person[i]->ped == person[person_index]->ped)) {
     if (person[i]->id > largest )
       largest = person[i]->id;
     --i;
   }

   i = person_index +1;
   while((person[i] != NULL) && (person[i]->ped == person[person_index]->ped)){
     if (person[i]->id > largest )
       largest = person[i]->id;
     ++i;
   }

   return(largest);
 }










/*****************************************************************************/
/*                                                                           */
/*                            add_loop                                       */
/*                                                                           */
/*****************************************************************************/

void add_loop(start_of_ped,old)
  s_intg start_of_ped;
  s_intg old;
{
  s_intg i;
  s_intg new;
  s_intg next_possible_id;
  s_intg pedigree;

  /* Get next possible id for this pedigree */

  next_possible_id = largest_id(old) + 1;
  if (next_possible_id > biggest_i_id) biggest_i_id = next_possible_id;

  /* Open a slot in the person array to insert a new person. */

  i = totperson;
  while(i > old) {
   person[i+1] = person[i];
   i--;
  }
  new = i + 1;
  totperson++;

  person[new] = (struct ind *) calloc ( 1, sizeof( struct ind ) );
  if ( person[new] == NULL ) {
    fprintf(stderr,"ERROR: Cannot allocate memory.\n");
    exit(1);
  }

  /* Copy the original record to the new record  */
  /* but 0 out the parents of the new record and */
  /* and the children of the old record.         */

  if (ped_integers)
    person[new]->oldped = person[old]->oldped;
  else
    strcpy(person[new]->oldped_s,person[old]->oldped_s);

  if (ind_integers)
    person[new]->oldid = person[old]->oldid;
  else
    strcpy(person[new]->oldid_s,person[old]->oldid_s);

  person[new]->id         = next_possible_id;
  person[new]->ped        = person[old]->ped;
  person[new]->paid       = 0;
  person[new]->maid       = 0;
  person[new]->pa         = NULL;
  person[new]->ma         = NULL;
  person[new]->offid      = person[old]->offid;
  person[new]->foff       = person[old]->foff;
  person[new]->nextpa     = NULL;
  person[new]->nextma     = NULL;
  person[new]->sex        = person[old]->sex;
  person[new]->proband    = 2;
  person[new]->phen       = person[old]->phen;

  /* zero out the children of the original record */
  /* and assign the proband.                      */

  person[old]->offid      = 0;
  person[old]->npaid      = 0;
  person[old]->nmaid      = 0;
  person[old]->foff       = NULL;
  person[old]->proband    = 2;

  /* Scan through entire pedigree and look for people who had a parent */
  /* that was the original and switch them to the newly created parent */

       i = start_of_ped;
       pedigree = person[start_of_ped]->ped;

       while( (i<=totperson) &&
            (pedigree == person[i]->ped)){

           if (person[i]->paid == person[old]->id) {
              person[i]->paid = person[new]->id;
              person[i]->pa = person[new];
            }
	   if (person[i]->maid == person[old]->id) {
              person[i]->maid = person[new]->id;
              person[i]->ma = person[new];
           }
         i++;
        } /* end of while */
}




/*****************************************************************************/
/*                                                                           */
/*                             file_loops                                    */
/*                                                                           */
/*****************************************************************************/

void file_loops()
{
  s_byte response;
  s_intg pedigree;
  s_intg pedigree_i;
  u_byte pedigree_s[maxname];
  s_intg person_i;
  s_byte person_s[maxname];
  s_intg i;
  s_intg from_file;
  s_intg found_start_of_ped;
  s_intg start_of_ped;
  s_intg found_person;
  s_intg good_format;
  s_byte loop_file[max_filespec];
  s_byte c;
  FILE   *loopf;
 
  fprintf(stdout,"\nEnter filename -> ");
  fscanf(stdin,"%s",loop_file);
  if ( (loopf = fopen(loop_file,"r")) == NULL) {
   fprintf(stderr,"ERROR: Cannot open file %s\n",loop_file);
   exit(1);
  }

  while(!feof(loopf)) {

   fscanf(loopf,"%s",pedigree_s);
   fscanf(loopf,"%s",person_s);

   if (!feof(loopf)) {

  /* convert to integer if needed */

  if (ped_integers) sscanf((s_byte*)pedigree_s,"%d",&pedigree_i);
  if (ind_integers) sscanf(person_s,"%d",&person_i);

  found_person = FALSE;
  found_start_of_ped = FALSE;
  i = 1;
  while(( i <= totperson) && (!found_person)) {

    if ((!found_start_of_ped) && (person[i]->oldped == pedigree_i)) {
      start_of_ped = i;
      found_start_of_ped = TRUE;
    }

    if ( (ped_integers) &&
         (ind_integers) &&
         (pedigree_i == person[i]->oldped) &&
         (person_i == person[i]->oldid)) {
     found_person = TRUE; 
     add_loop(start_of_ped,i);
    }
    else
    if ( (!ped_integers) &&
         (ind_integers) &&
         (strcmp_i(pedigree_s,person[i]->oldped_s) == 0) &&
         (person_i == person[i]->oldid)) {
       found_person = TRUE; 
       add_loop(start_of_ped,i);
      }
      else
      if ( (ped_integers) &&
             (!ind_integers) &&
             (pedigree_i == person[i]->oldped)  &&
             (strcmp_i(person_s,person[i]->oldid_s ) == 0)) {
         found_person = TRUE; 
         add_loop(start_of_ped,i);
        }
        else
        if ( (!ped_integers) &&
             (!ind_integers) &&
             (strcmp_i(pedigree_s,person[i]->oldped_s) == 0) &&
             (strcmp_i(person_s,person[i]->oldid_s ) == 0)) {
           found_person = TRUE; 
           add_loop(start_of_ped,i);
        }

        i++;
	if (( i>totperson) && (!found_person)) {
  fprintf(stderr,"ERROR: Ped: %s Per: %s - Not found, check loop file.\n",
		  pedigree_s,
		  person_s);
	  exit(1);
	}
      }
     }
    }
     fclose(loopf);
}







/*****************************************************************************/
/*                                                                           */
/*                             some_loops                                    */
/*                                                                           */
/*****************************************************************************/

void some_loops()
{
  s_byte response;
  s_intg person_i;
  s_byte person_s[maxname];
  s_intg pedigree_i;
  s_byte pedigree_s[maxname];
  s_intg pedigree;
  s_intg i,j,k;
  s_intg start_of_ped;
  s_intg found_per;
  s_intg found_ped;
  s_intg good_format;
  s_intg count;             /* Number of people assigned loops.  */
  s_byte done;

  count = 0;

  fprintf(stdout,"\n\n\tEnter identifiers for ");
  fprintf(stdout,"each pedigree and person...\n");
  fprintf(stdout,"\tenter pedigree 0 when finished.\n");

  done = FALSE;
  while(!done) {

    read_pedigree(&pedigree_i,pedigree_s);

    if (pedigree_s[0] == '0') done = TRUE;
     else {

       found_ped = FALSE;
       i = 1;
	 while ((!found_ped) && (i<=totperson)) {

           if (ped_integers) {
             if (person[i]->oldped == pedigree_i) {
               found_ped = TRUE;
               start_of_ped = i;
               pedigree = person[i]->ped;
             }
           }

           else {  /* pedigree id's are strings */
             if (( strcmp_i(pedigree_s,person[i]->oldped_s) == 0)) {
	       found_ped = TRUE;
	       start_of_ped = i;
	       pedigree = person[i]->ped;
	     }
           }
	   i++;
         }

       if ((i>= totperson) && (!found_ped))
	 fprintf(stdout,"\tPedigree not found...\n");
       }

       if((!done) && (found_ped)) {
        found_per = FALSE;
        i = start_of_ped;
        while(!found_per) {

         read_person(&person_i,person_s);

	 while( (i<=totperson) &&
	       (pedigree == person[i]->ped) &&
	       (!found_per)){
		 if (ind_integers) {
		   if (person[i]->oldid == person_i) {
                     loops[count++] = i;
		     found_per = TRUE;
                     add_loop(start_of_ped,i);
		   }
		 }
		 else
		   if (!strcmp_i(person[i]->oldid_s,person_s)) {
                     loops[count++] = i;
		     found_per = TRUE;
                     add_loop(start_of_ped,i);
		   }

		 if (!found_per) i++;
	       } /* end of while */

	 /* The whole pedigree has been searched.  If no match  */
	 /* was found for the identifier then set i back to the */
	 /* begining of the pedigree and re-prompt the user for  */
	 /* for a new id.                                       */

	 if (!found_per) {
	   fprintf(stdout,"\tPerson not found...\n");
	   i = start_of_ped;
	 }
       }  /* end of while(!found_per) */
       } /* end of if(!done ) */
  }   /* end of while(!done) */

  save_loops(count);
}




/*****************************************************************************/
/*                                                                           */
/*                          get_loops                                        */
/*                                                                           */
/*****************************************************************************/

void get_loops()
{
  s_byte response;

  fprintf(stdout,"\n");
  fprintf(stdout,"Does your pedigree file contain any loops?    (y/n) -> ");
  fscanf(stdin,"%1s",&response);
  if ((response == 'y') || (response == 'Y')) {
    fprintf(stdout,"Do you have a file of loop assignments?       (y/n) -> ");
    fscanf(stdin,"%1s",&response);
    if ((response == 'y') || (response == 'Y')) file_loops();
    else some_loops();
  }
}



/****************************************************************************/
/*                                                                          */
/*			count_generations				    */
/*                                                                          */
/****************************************************************************/

s_intg count_generations(person_x)
  s_intg person_x;
{
  struct ind *decendent;
  s_intg count;

  count = 0;
  if ( person[person_x]->foff != NULL ) {
    count++;
    decendent = person[person_x]->foff;
    while (decendent->foff != NULL) {
      count++;
      decendent = decendent->foff;
    }
  }
  return(count);
}




/****************************************************************************/
/*                                                                          */
/*			   clear_proband   				    */
/*                                                                          */
/* Given the index of some person, this function clears al proband fields   */
/* in the same pedigree.                                                    */
/*                                                                          */
/****************************************************************************/

void clear_proband(person_index)
     s_intg person_index;
{

  s_intg pedigree;
  s_intg i;
  s_intg found_ped;
  s_intg end_of_ped;

  pedigree = person[person_index]->ped;

  found_ped = FALSE;
  i=1;
  while((!found_ped) && (i<=totperson))
   if (pedigree == person[i++]->ped) found_ped = TRUE;

  if (!cleared[pedigree]){
    i--;
    end_of_ped = FALSE;
    while ((!end_of_ped) && (i<=totperson)) {
      if (pedigree == person[i]->ped) {
	if (person[i]->proband != 2)  /* dont clear a loop person */
	  person[i]->proband = 0;
      }
      else end_of_ped = TRUE;
      i++;
    }
    cleared[pedigree] = TRUE;
  }
}





/*****************************************************************************/
/*                                                                           */
/*                          save_probands                                    */
/*                                                                           */
/*****************************************************************************/

void save_probands(count)
    s_intg count;
{
  s_intg i;
  s_byte response;
  s_byte proband_file[max_filespec];
  FILE   *prof;


  fprintf(stdout,"\n\nDo you want these selections saved ");
  fprintf(stdout,"for later use?  (y/n) -> ");
  fscanf(stdin,"%1s",&response);

  if ((response == 'y') || (response == 'Y')) {
    fprintf(stdout,"\nEnter filename -> ");
    fscanf(stdin,"%s",proband_file);
    if ( (prof = fopen(proband_file,"w")) == NULL) {
      fprintf(stderr,"ERROR: Cannot open file %s\n",proband_file);
    }
    else {
      for(i=0; i<count; i++) {

       if (ped_integers)
         fprintf(prof,"%d ",person[ probands[i] ]->oldped);
       else
         fprintf(prof,"%s ",person[ probands[i] ]->oldped_s);

       if (ind_integers)
         fprintf(prof,"%d\n",person[ probands[i] ]->oldid);
       else
         fprintf(prof,"%s\n",person[ probands[i] ]->oldid_s);

      }
      fclose(prof);
    }
  }
 }




/****************************************************************************/
/*                                                                          */
/*			   auto_probands   				    */
/*                                                                          */
/****************************************************************************/

void auto_probands()
{
  s_intg i;
  s_intg ped_num;		/* current pedigree */
  s_intg found;                 /* index of best choice so far */
  s_intg max_level;             /* max num generations found so far */

/* For each male with unknown parents, count the number */
/* of generations which are below him.                  */

  for(i=1; i<=totperson; i++)
   if ((person[i]->paid == 0) &&
       (person[i]->maid == 0) &&
       (person[i]->sex  == 1)) person[i]->generations = count_generations(i);

/* For each pedigree find the first male that has unknown */
/* parents and more generations below him than any other  */
/* male with unknown parents.                             */

i = 1;
while(i<=totperson) {
  found = i;
  ped_num = person[i]->ped;
  max_level = 0;  
  while((i<=totperson) && (ped_num == person[i]->ped)) {
    if ( (person[i]->paid == 0) &&
         (person[i]->maid == 0) &&
         (person[i]->sex  == 1) &&
         (person[i]->generations > max_level)) {
      found = i;
      max_level = person[i]->generations;
    }
    i++;
  }
  person[found]->proband = 1;
 }
}





/*****************************************************************************/
/*                                                                           */
/*                             file_probands                                 */
/*                                                                           */
/*****************************************************************************/

void file_probands()
{
  s_byte response;
  s_intg pedigree;
  s_intg pedigree_i;
  s_byte pedigree_s[maxname];
  s_intg person_i;
  s_byte person_s[maxname];
  s_intg i,j,k;
  s_intg from_file;
  s_intg start_of_ped;
  s_intg found;
  s_intg good_format;
  s_byte proband_file[max_filespec];
  s_byte c;
  FILE   *prof;
 
      fprintf(stdout,"\nEnter filename -> ");
      fscanf(stdin,"%s",proband_file);
      if ( (prof = fopen(proband_file,"r")) == NULL) {
       fprintf(stderr,"ERROR: Cannot open file %s\n",proband_file);
       exit(1);
      }

      auto_probands();  /* makes sure a proband is set for each pedigree */
                        /* even if the input file does not set it.       */


      while(!feof(prof)) {

       fscanf(prof,"%s",pedigree_s);
       fscanf(prof,"%s",person_s);

       if (!feof(prof)) {

       /* convert to integer if needed */

       if (ped_integers) sscanf(pedigree_s,"%d",&pedigree_i);
       if (ind_integers) sscanf(person_s,"%d",&person_i);

       found = FALSE;
       j = 1;
       while(( j<=totperson) && (!found)) {

        if ( (ped_integers) &&
             (ind_integers) &&
             (pedigree_i == person[j]->oldped) &&
             (person_i == person[j]->oldid)) {
          clear_proband(j);
          person[j]->proband = 1;
          found = TRUE;
        }
        else
        if ( (!ped_integers) &&
             (ind_integers) &&
             (strcmp_i(pedigree_s,person[j]->oldped_s) == 0) &&
             (person_i == person[j]->oldid)) {
          clear_proband(j);
          person[j]->proband = 1;
          found = TRUE;
        }
        else
        if ( (ped_integers) &&
             (!ind_integers) &&
             (pedigree_i == person[j]->oldped)  &&
             (strcmp_i(person_s,person[j]->oldid_s ) == 0)) {
          clear_proband(j);
          person[j]->proband = 1;
          found = TRUE;
        }
        else
        if ( (!ped_integers) &&
             (!ind_integers) &&
             (strcmp_i(pedigree_s,person[j]->oldped_s) == 0) &&
             (strcmp_i(person_s,person[j]->oldid_s ) == 0)) {
          clear_proband(j);
          person[j]->proband = 1;
          found = TRUE;
        }

        j++;
	if (( j>totperson) && (!found)) {
  fprintf(stderr,"ERROR: Ped: %s Per: %s - Not found, check proband file.\n",
		  pedigree_s,
		  person_s);
	  exit(1);
	}
       }
      }
     }
     fclose(prof);
}



/*****************************************************************************/
/*                                                                           */
/*                              all_probands                                 */
/*                                                                           */
/*****************************************************************************/

void all_probands()
{
  s_byte response;
  s_intg pedigree;
  s_byte person_s[maxname];
  s_intg person_i;
  s_intg i,j,k;
  s_intg start_of_ped;
  s_intg found;
  s_intg good_format;
  s_intg count;               /* number of people assigned probands.  */

  count = 0;

      fprintf(stdout,"\n\nEnter the identifier of the ");
      fprintf(stdout,"person who is to be the proband for...\n\n");
      pedigree = 0;
      for(i=1; i<=totperson; i++) {
        if (pedigree != person[i]->ped) {
          start_of_ped = i;
          pedigree = person[i]->ped;

          fprintf(stdout,"\n\n\tPedigree   -> ");
          if (ped_integers)
            fprintf(stdout,"%d\n",person[i]->oldped);
          else
            fprintf(stdout,"%s\n",person[i]->oldped_s);

          read_person(&person_i,person_s);

          found = FALSE;
          while( (i<=totperson) &&
                 (pedigree == person[i]->ped) &&
                 (!found)){

            if (ind_integers) {
                if (person[i]->oldid == person_i) {
                  person[i]->proband = 1;
                  probands[count++] = i;
                  found = TRUE;
                }
            }
            else
             if (!strcmp_i(person[i]->oldid_s,person_s)) {
               person[i]->proband = 1;
               probands[count++] = i;
               found = TRUE;
             }

             if (!found) i++;
           } /* end of while */

             /* The whole pedigree has been searched.  If no match  */
             /* was found for the identifier then set i back to the */
             /* begining of the pedigree and re-prompt the user for  */
             /* for a new id.                                       */

             if (!found) {
               fprintf(stdout,"\tPerson not found...\n");
               i = start_of_ped -1;  /* -1 because the for loop will +1 */
               pedigree = 0;
             }
        }
      }

  save_probands(count);
}

/*****************************************************************************/
/*                                                                           */
/*                              some_probands                                 */
/*                                                                           */
/*****************************************************************************/

void some_probands()
{
  s_byte response;
  s_intg person_i;
  s_byte person_s[maxname];
  s_intg pedigree_i;
  s_byte pedigree_s[maxname];
  s_intg pedigree;
  s_intg i,j,k;
  s_intg start_of_ped;
  s_intg found_per;
  s_intg found_ped;
  s_intg good_format;
  s_intg count;             /* number of people assigned probands.  */
  s_byte done;

  auto_probands();      /* Set all probands before the user selects a few. */
  count = 0;            /* Count number of successful entries. */

  fprintf(stdout,"\n\n\tEnter identifiers for ");
  fprintf(stdout,"each pedigree and person...\n");
  fprintf(stdout,"\tenter pedigree 0 when finished.\n");

  done = FALSE;
  while(!done) {

    read_pedigree(&pedigree_i,pedigree_s);

    if (pedigree_s[0] == '0') done = TRUE;
     else {

       found_ped = FALSE;
       i = 1;
	 while ((!found_ped) && (i<=totperson)) {

           if (ped_integers) {
             if (person[i]->oldped == pedigree_i) {
               found_ped = TRUE;
               start_of_ped = i;
               pedigree = person[i]->ped;
             }
           }

           else {  /* pedigree id's are strings */
             if (( strcmp_i(pedigree_s,person[i]->oldped_s) == 0)) {
	       found_ped = TRUE;
	       start_of_ped = i;
	       pedigree = person[i]->ped;
	     }
           }
	   i++;
         }

       if ((i>= totperson) && (!found_ped))
	 fprintf(stdout,"\tPedigree not found...\n");
       }

       if((!done) && (found_ped)) {
        found_per = FALSE;
        i = start_of_ped;
        while(!found_per) {

          read_person(&person_i,person_s);

	 while( (i<=totperson) &&
	       (pedigree == person[i]->ped) &&
	       (!found_per)){
		 if (ind_integers) {
		   if (person[i]->oldid == person_i) {
		     clear_proband(i);
		     person[i]->proband = 1;
                     probands[count++] = i;
		     found_per = TRUE;
		   }
		 }
		 else
		   if (!strcmp_i(person[i]->oldid_s,person_s)) {
		     clear_proband(i);
		     person[i]->proband = 1;
                     probands[count++] = i;
		     found_per = TRUE;
		   }

		 if (!found_per) i++;
	       } /* end of while */

	 /* The whole pedigree has been searched.  If no match  */
	 /* was found for the identifier then set i back to the */
	 /* begining of the pedigree and re-prompt the user for  */
	 /* for a new id.                                       */

	 if (!found_per) {
	   fprintf(stdout,"\tPerson not found...\n");
	   i = start_of_ped;
	 }
       }  /* end of while(!found_per) */
       } /* end of if(!done ) */
  }   /* end of while(!done) */

  save_probands(count);
}






/*****************************************************************************/
/*                                                                           */
/*                          get_probands                                     */
/*                                                                           */
/*****************************************************************************/

void get_probands()
{
  s_byte response;

  fprintf(stdout,"\n");
  fprintf(stdout,"Do you want probands selected automaticaly?   (y/n) -> ");
  fscanf(stdin,"%1s",&response);
  if ((response == 'y') || (response == 'Y')) auto_probands();
  else {

  fprintf(stdout,"Do you have a file of proband assignments?    (y/n) -> ");
  fscanf(stdin,"%1s",&response);
  if ((response == 'y') || (response == 'Y')) file_probands();
  else {
 
  fprintf(stdout,"Do you want to select all probands?           (y/n) -> ");
  fscanf(stdin,"%1s",&response);
  if ((response == 'y') || (response == 'Y')) all_probands();
  else some_probands();
}
}
}


/****************************************************************************/
/*                                                                          */
/*                           writeped                                       */
/*                                                                          */
/****************************************************************************/

void writeped()
{
  s_intg i;
  s_byte *ped_format;
  s_byte *ind_format;
  s_byte *format_1  = "%1d";
  s_byte *format_2  = "%2d";
  s_byte *format_3  = "%3d";
  s_byte *format_4  = "%4d";
  s_byte *format_5  = "%5d";
  s_byte *format_6  = "%6d";
  s_byte *format_7  = "%7d";
  s_byte *format_8  = "%8d";
  s_byte *format_9  = "%9d";
  s_byte *format_10 = "%10d";

  /* setup pedigree printing format */

  if (biggest_p_id >= 10e10)
    ped_format = format_10;
  else
  if (biggest_p_id >= 10e09)
    ped_format = format_9;
  else
 if (biggest_p_id >= 10e08)
    ped_format = format_8;
  else
  if (biggest_p_id >= 10e07)
    ped_format = format_7;
  else;
  if (biggest_p_id >= 10e06)
    ped_format = format_6;
  else
  if (biggest_p_id >= 10e05)
    ped_format = format_5;
  else
  if (biggest_p_id >= 10e04)
    ped_format = format_4;
  else
  if (biggest_p_id >= 10e03)
    ped_format = format_3;
  else
  if (biggest_p_id >= 10e02)
    ped_format = format_2;
  else
    ped_format = format_1;

  /* setup individual printing format */

  if (biggest_i_id > 9999)
    ind_format = format_6;
  else
  if (biggest_i_id > 999)
    ind_format = format_5;
  else
  if (biggest_i_id > 99)
    ind_format = format_4;
  else
  if (biggest_i_id > 9)
    ind_format = format_3;
  else
    ind_format = format_2;



  for( i=1; i<=totperson; i++) {
    if (ped_integers )
      fprintf(pedout,ped_format,person[i]->oldped);
    else
      fprintf(pedout,ped_format,  person[i]->ped);

    fprintf(pedout,ind_format,  person[i]->id);

    if (person[i]->pa != NULL)
      fprintf(pedout,ind_format,  person[i]->pa->id);
    else
      fprintf(pedout,ind_format,  0);

    if (person[i]->ma != NULL)
      fprintf(pedout,ind_format,  person[i]->ma->id);
    else
      fprintf(pedout,ind_format,  0);

    if (person[i]->foff != NULL)
      fprintf(pedout,ind_format,  person[i]->foff->id);
    else
      fprintf(pedout,ind_format,  0);

    if (person[i]->nextpa != NULL)
      fprintf(pedout,ind_format,  person[i]->nextpa->id);
    else
      fprintf(pedout,ind_format,  0);

    if (person[i]->nextma != NULL)
      fprintf(pedout,ind_format,  person[i]->nextma->id);
    else
      fprintf(pedout,ind_format,  0);

    fprintf(pedout,"%2d", person[i]->sex);
    fprintf(pedout,"%2d", person[i]->proband);
    fprintf(pedout,"%s", person[i]->phen->phen_chars);

    if (ped_integers)
      fprintf(pedout,"  Ped: %d",person[i]->oldped);
    else
      fprintf(pedout,"  Ped: %s",person[i]->oldped_s);

    if (ind_integers)
      fprintf(pedout,"  Per: %d\n",person[i]->oldid);
    else
      fprintf(pedout,"  Per: %s\n",person[i]->oldid_s);
  }

}





/*****************************************************************************/
/*                                                                           */
/*                             strcmp_i                                      */
/*                                                                           */
/* A case insensitive string compare. Returns 0 if equal, 1 if not.          */
/*                                                                           */
/*****************************************************************************/

s_intg strcmp_i(s1,s2)
     register u_byte *s1;
     register u_byte *s2;
{
  while(((*s1 >= 'a' && *s1 <= 'z')?
	 (*s1 & 0xdf):
	 (*s1))==
	((*s2 >= 'a' && *s2 <= 'z')?
	 (*s2++ & 0xdf):
	 (*s2++)))
    if (*s1++ == '\0')
      return(0);
  return(1);
}



/*****************************************************************************/
/*                                                                           */
/*                               check_sequential                            */
/*                                                                           */
/* If person identifiers are integers you must verify that there are no      */
/* missing sequential id's.  Exit the program if any are found because       */
/* the other data-checking routines will blow up if there is a missing id.   */
/* If id's are strings then substitute id's are assigned automaticaly and    */
/* there shouldn't be any problem.                                           */
/*                                                                           */
/* This function was changed to return an error code so that the program can */
/* reset itself and go with string id's for individuals.  (When string id's  */
/* are used sequential id's are assigned automatically.                      */
/*                                                                           */
/*****************************************************************************/

/* void check_sequential() */
s_intg check_sequential()
{

  s_intg i;
  u_intg bailout;
  s_intg missing;

  bailout = FALSE;
  missing = 0;
  for( i=1; i<=totperson+missing; i++) {
    if (person[i]==NULL) {
      if (ped_integers){
        while((person[i]==NULL) && (i<=totperson+missing)) {
	  i++;
	  missing++;
	}
        fprintf(stderr,"WARNING: Ped: %d  Per: %d - Out of sequence.\n",
	          person[i]->oldped,person[i]->oldid);
	bailout = TRUE;
      }
      else {
        while((person[i]==NULL) && (i<=totperson+missing)) {
	  i++;
	  missing++;
	}
	fprintf(stderr,"ERROR: Ped: %s  Per: %d - Out of sequence.\n",
		person[i]->oldped_s,person[i]->oldid);
	bailout = TRUE;
      }
    }
  }
  if (bailout) return(1);
  else return(0);
}

/*****************************************************************************/
/*                                                                           */
/*                               check_sex                                   */
/*                                                                           */
/*****************************************************************************/

void check_sex()
{
  s_intg i;

  for(i=1; i<=totperson; i++) {

    /* Verify that each person has either 0 or 2 parents. */

    if(((person[i]->pa != NULL) && (person[i]->ma == NULL)) ||
       ((person[i]->pa == NULL) && (person[i]->ma != NULL))){
     if ((ped_integers) && (ind_integers)){
       fprintf(stderr,"ERROR: Ped: %d  Per: %d - Only one parent.\n",
	      person[i]->oldped,person[i]->oldid);
       found_error = TRUE;
     }
     else if ((ped_integers) && (!ind_integers)) {
       fprintf(stderr,"ERROR: Ped: %d  Per: %s - Only one parent.\n",
	      person[i]->oldped,person[i]->oldid_s);
       found_error = TRUE;
     }
     else if ((!ped_integers) && (ind_integers)) {
       fprintf(stderr,"ERROR: Ped: %s  Per: %d - Only one parent.\n",
	      person[i]->oldped_s,person[i]->oldid);
       found_error = TRUE;
     }
     else if ((!ped_integers) && (!ind_integers)) {
       fprintf(stderr,"ERROR: Ped: %s  Per: %s - Only one parent.\n",
	      person[i]->oldped_s,person[i]->oldid_s);
       found_error = TRUE;
     }

    }

    /* Verify that father is male. */

   if((person[i]->pa != NULL) && (person[i]->pa->sex != 1)) {
     if ((ped_integers) && (ind_integers)) {
       fprintf(stderr,"ERROR: Ped: %d  Per: %d - Sex of father.\n",
	      person[i]->oldped,person[i]->oldid);
       found_error = TRUE;
     }
     else if ((ped_integers) && (!ind_integers)) {
       fprintf(stderr,"ERROR: Ped: %d  Per: %s - Sex of father.\n",
	      person[i]->oldped,person[i]->oldid_s);
       found_error = TRUE;
     }
     else if ((!ped_integers) && (ind_integers)) {
       fprintf(stderr,"ERROR: Ped: %s  Per: %d - Sex of father.\n",
	      person[i]->oldped_s,person[i]->oldid);
       found_error = TRUE;
     }
     else if ((!ped_integers) && (!ind_integers)) {
       fprintf(stderr,"ERROR: Ped: %s  Per: %s - Sex of father.\n",
	      person[i]->oldped_s,person[i]->oldid_s);
       found_error = TRUE;
     }
    }

    /* Verify that mother is female. */

   if((person[i]->ma != NULL) && (person[i]->ma->sex != 2)) {
     if ((ped_integers) && (ind_integers)) {
       fprintf(stderr,"ERROR: Ped: %d  Per: %d - Sex of mother.\n",
	      person[i]->oldped,person[i]->oldid);
       found_error = TRUE;
     }
     else if ((ped_integers) && (!ind_integers)) {
       fprintf(stderr,"ERROR: Ped: %d  Per: %s - Sex of mother.\n",
	      person[i]->oldped,person[i]->oldid_s);
       found_error = TRUE;
     }
     else if ((!ped_integers) && (ind_integers)) {
       fprintf(stderr,"ERROR: Ped: %s  Per: %d - Sex of mother.\n",
	      person[i]->oldped_s,person[i]->oldid);
       found_error = TRUE;
     }
     else if ((!ped_integers) && (!ind_integers)) {
       fprintf(stderr,"ERROR: Ped: %s  Per: %s - Sex of mother.\n",
	      person[i]->oldped_s,person[i]->oldid_s);
       found_error = TRUE;
     }
    }
  }
}





/*****************************************************************************/
/*                                                                           */
/*                           check_no_phen                                   */
/*                                                                           */
/* If a person has no phenotypic data then they were referenced as a parent  */
/* in the datafile but did not have a record of their own.                   */
/*                                                                           */
/*****************************************************************************/

void check_no_phen()
{

  s_intg i;

  for(i=1; i<=totperson; i++) {

   if (person[i]->phen == NULL){
     if ((ped_integers) && (ind_integers)) {
       fprintf(stderr,"ERROR: Ped: %d  Per: %d - No data.\n",
	      person[i]->oldped,person[i]->oldid);
       found_error = TRUE;
     }
     else if ((ped_integers) && (!ind_integers)) {
       fprintf(stderr,"ERROR: Ped: %d  Per: %s - No data.\n",
	      person[i]->oldped,person[i]->oldid_s);
       found_error = TRUE;
     }
     else if ((!ped_integers) && (ind_integers)) {
       fprintf(stderr,"ERROR: Ped: %s  Per: %d - No data.\n",
	      person[i]->oldped_s,person[i]->oldid);
       found_error = TRUE;
     }
     else if ((!ped_integers) && (!ind_integers)) {
       fprintf(stderr,"ERROR: Ped: %s  Per: %s - No data.\n",
	      person[i]->oldped_s,person[i]->oldid_s);
       found_error = TRUE;
     }
   }
 }
}




/*****************************************************************************/
/*                                                                           */
/*                           check_no_family                                 */
/*                                                                           */
/* If a person has no parents and no children then report it.                */
/*                                                                           */
/*****************************************************************************/

void check_no_family()
{
  u_intg i;

  /* Mark each person that is referenced as a parent of another person. */

  for( i=1; i<=totperson; i++) {
    if (person[i]->pa != NULL) person[i]->pa->is_parent = TRUE;
    if (person[i]->ma != NULL) person[i]->ma->is_parent= TRUE;
  }

  /* Report all persons who do not have their is_parent  */
  /* flag set and do not have any parents of their own.  */

  for( i=1; i<=totperson; i++ ) {
    if ((!person[i]->is_parent) &&
        (person[i]->pa == NULL) &&
        (person[i]->ma == NULL))
      if ((ped_integers) && (ind_integers)) {
	fprintf(stderr,"ERROR: Ped: %d  Per: %d - No family.\n",
		person[i]->oldped,person[i]->oldid);
	found_error = TRUE;
      }
      else if ((ped_integers) && (!ind_integers)) {
	fprintf(stderr,"ERROR: Ped: %d  Per: %s - No family.\n",
		person[i]->oldped,person[i]->oldid_s);
	found_error = TRUE;
      }
      else if ((!ped_integers) && (ind_integers)) {
	fprintf(stderr,"ERROR: Ped: %s  Per: %d - No family.\n",
		person[i]->oldped_s,person[i]->oldid);
	found_error = TRUE;
      }
      else if ((!ped_integers) && (!ind_integers)) {
	fprintf(stderr,"ERROR: Ped: %s  Per: %s - No family.\n",
		person[i]->oldped_s,person[i]->oldid_s);
	found_error = TRUE;
      }
  }
}


/*****************************************************************************/
/*                                                                           */
/*                           reset_for_strings                               */
/*                                                                           */
/*****************************************************************************/

void reset_for_strings()
{
  u_intg i;

  force_strings = TRUE;
  fclose(pedfile);
  fclose(pedout);

  for (i=1; i<=totperson; i++) {
    if (person[i] != NULL) {
      free (person[i]->phen);
      free (person[i]);
    }
  }
}

/*****************************************************************************/
/*                                                                           */
/*                                 main                                      */
/*                                                                           */
/*****************************************************************************/

main(argc,argv)
     s_intg argc;
     char *argv[];
{

  u_byte response;

  fprintf(stdout,"\n           MAKEPED Version %3.1f\n\n",mversion);

  if (argc > 3) {
    fprintf(stderr,"ERROR: Two many command line arguments");
    exit(1);
  }

  pifile[0] = '\0';
  pofile[0] = '\0';

  if( argc > 1) {                        /* GET FILESPEC IF ON COMMAND LINE */
    strcpy((char*) pifile, argv[1]);
  }

  else{                                 /* FILES ARE NOT ON COMMAND LINE */
    while ( pifile[0] == '\0' ) {
      fprintf(stdout,"Pedigree file -> ");
      gets  ((s_byte*)pifile);
    }
    while ( pofile[0] == '\0') {
      fprintf(stdout,"Output file   -> ");
      gets  ((s_byte*)pofile);
    }
  }


  if( argc > 2) {                       /* GET OUTPUT FILE IF ON COM LINE */
    strcpy((char*) pofile, argv[2]);
  }
  else {
    while ( pofile[0] == '\0') {
      fprintf(stdout,"Output file   -> ");
      gets  ((s_byte*)pofile);
    }
  }

  force_strings = FALSE;
  found_error = FALSE;

 try_again:
  if ((pedfile = fopen((s_byte*)pifile, "r")) == NULL){
   fprintf(stderr,"ERROR: Cannot open %s\n",pifile);
   exit(1);
 }

  if ((pedout = fopen((s_byte*)pofile, "w")) == NULL){
   fprintf(stderr,"ERROR: Cannot open %s\n",pofile);
   exit(1);
 }

  readped();
	
  /* As an after-thought, it was decided that if a sequential-id error was */
  /* found that the program should reset itself and run with string id's,  */
  /* goto try_again takes care of this.                                    */

  if ((ind_integers) && (check_sequential())) {
      reset_for_strings();
      goto try_again;
    }

  check_sex();
  check_no_phen();
/*  check_no_family();*/
  if(found_error) exit(1);
  pointers();
  get_loops(); 
  get_probands();
  writeped();

  fclose(pedfile);
  fclose(pedout);
  exit (0);
}
