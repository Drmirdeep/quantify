#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

void complement(char s[]){
  int i;
  for (i = 0;i < strlen(s); i++) {
	if(s[i] == 'a'){
	  s[i] = 't';
	}else if(s[i] == 'c'){
	  s[i] = 'g';
	}else if(s[i] == 'g'){
	  s[i] = 'c';
	}else if(s[i] == 't'){
	  s[i] = 'a';
	}else if(s[i] == 'u'){
	  s[i] = 'a';
	}else if(s[i] == 'A'){
	  s[i] = 'T';
	}else if(s[i] == 'C'){
	  s[i] = 'G';
	}else if(s[i] == 'G'){
	  s[i] = 'C';
	}else if(s[i] == 'T'){
	  s[i] = 'A';
	}else if(s[i] == 'U'){
	  s[i] = 'A';
	}else if(s[i] == '.'){
		s[i] = 'N';
	}else{
	  
	  ;
	}
  }
}

void lower(char s[]){
  int i;
  for (i = 0;i < strlen(s); i++) {
	if(s[i] < 97){
	  s[i]+=32;
	}
  }
}

void aclean(char s[]){
  int j=0;
  while(s[j] != '\0'){
	s[j] = '\0';
	j++;
  }
}

void reverse(char s[])
{
      int c, i, j;
      for (i = 0, j = strlen(s)-1; i < j; i++, j--) {
         c = s[i];
         s[i] = s[j];
         s[j] = c;
      }
}





void rcomplement(char s[]){
  reverse(s);
  complement(s);
}





int main(int argc, char *argv[]){
  if(argc < 4){
	printf("\nUsage:\nprepro_deep_seq_reads  max_read_len  read_file outfile [options]\nadd_offset_of_read_id");
	return 0;
  }
  FILE *file;




  if(fopen(argv[2],"r") == NULL){
	printf("could not open file %s for reading\n",argv[2]);
	return 0;
  }else{
	file=fopen(argv[2],"r");
  }

  int RLEN=atoi(argv[1]);  
  /*   create a file pointer */
  /*   now open file for reading to get the current index */
  /*   variable declarations */
  char id[RLEN+1]; /* create id array of length 200 */
  char strand;
  char chr[RLEN+1];
  int start;
  char seq[RLEN+1];
  char gseq[RLEN+1];
  char qual[RLEN+1];
  char changes[RLEN+1];
  char edit[RLEN+1];
  int freq;
  int mm=0;
  int len;
  /*   make outfile pointer */
  FILE *outfile;
  /*   either append or create if not exists yet */
  if(atoi(argv[3]) == 1){
	  outfile=stdout;
  }else{
	  outfile=fopen(argv[3],"w");
  }
  
  /*   now parse and print as fastq */
  int count=0;
  char tmp[RLEN+1];
  //  while(fscanf(file,"%s %c %s %d %s %s %d %s",id,&strand,chr,&start,seq,qual,&freq,changes)!=EOF){
  // printf("Die max. Groesse des Puffers: %d\n",BUFSIZ); 
  //static char puffer[BUFSIZ];
  int offset=0;
  if(argc > 4){ offset=atoi(argv[4]);}

  

 setvbuf(file,NULL,_IOFBF,0); 
// printf("Die max. Groesse des Puffers: %d\n",BUFSIZ);

  while(fscanf(file,"%s",tmp)!=EOF){
	count++;
	if(count==1){
	  strcpy(id,tmp);
	}else if(count == (2+offset)){
	  strand=tmp[0];
	  
	  //	  printf("%c\n",strand);return 0;
	}else if(count == (3+offset)){
	  strcpy(chr,tmp);
	}else if(count == (4+offset)){
	  start=atoi(tmp);
	}else if(count == (5+offset)){
	  strcpy(seq,tmp);
	  if(strand == '-'){
		rcomplement(seq);
	  }
	  lower(seq);
	  strcpy(gseq,seq);
	  
	  // set edit string
	  int i=0;
	  while(seq[i] != '\0'){
		edit[i]='m';
		i++;
	  }
	  len=i;
	  edit[i]='\0';
	  
	}else if(count == (6+offset)){
	  strcpy(qual,tmp);

	}else if(count == (7+offset)){
	  freq=atoi(tmp);

	}else if(count == (8+offset) && (tmp[0] > 47 && tmp[0] < 58)){
	  count=0;
	  strcpy(changes,tmp);
	  if(strand == '-') complement(changes);
	// do this stuff here
		  
	  int i=0;
	  int pos=0;
	  char startp[20];
//	  	printf("%s\n",changes);
	  //	return 0;
	  
	  while(changes[i] != '\0'){
		if(changes[i] > 47 && changes[i] < 58){
		  startp[pos]=changes[i];
		  pos++;
		  i++;
		}else if(changes[i] == ':'){
		  i++;
		  int val=atoi(startp);
		 // 		  printf("%d\n",val);
		  //       return 0;
		  //printf("%d\t%d\n",val,i-1);
		  //printf("%s\n",changes);
		  //if(strand == '-') complement(changes);
		  //printf("%s\n",changes);return 0;
		  gseq[val]=changes[i]+32;
		  edit[val]='M';
		  mm++;
		  // if we reach the next entry now
		}else if(changes[i] == ','){
		  int j=0;
		  pos=0;
		  // clean array now
		  while(startp[j] != '\0'){
			startp[j]='\0';
		  j++;
		  }
		  i++;
		}else{i++;}
	  }
	  fprintf(outfile,"%s\t%d\t1\t%d\t%s\t%s\t%d\t%d\t%d\t%s\t%c\t%d\t%s\n",id,len,len,seq,chr,len,start+1,start+len,gseq,strand,mm,edit);

	  aclean(changes);
	  aclean(edit);
	  aclean(startp);
	
	  mm=0;


	  
	  // no mismatches occured so we are here and have to print out now the previous line
	}else if(count == (offset+8)){// this is the new id because there were no mismatches
	  fprintf(outfile,"%s\t%d\t1\t%d\t%s\t%s\t%d\t%d\t%d\t%s\t%c\t%d\t%s\n",id,len,len,seq,chr,len,start+1,start+len,gseq,strand,mm,edit);
	  
	  count=1;
	  strcpy(id,tmp);
	}


	
	//	printf("%s\t%d\t1\t%d\t%s\t%s\t%d\t%d\t%d\t%s\t%c\t%d\t%s\n",id,len,len,seq,chr,len,start,start+len,gseq,strand,mm,edit);
	




	//cind++;

  }
  // only print if last line in file does not have any edits otherwise it was printed already
  
  if(count == (7+offset))
	fprintf(outfile,"%s\t%d\t1\t%d\t%s\t%s\t%d\t%d\t%d\t%s\t%c\t%d\t%s\n",id,len,len,seq,chr,len,start+1,start+len,gseq,strand,mm,edit);
  fclose(file);
  fclose(outfile);
  //  cf=fopen(argv[3],"w");
  //fprintf(cf,"%d\n",cind);
  //fclose(cf);
  return 0;
}
