/*

The MIT License

Copyright (c) 2007 Guillermo L. Taboada (taboada@udc.es)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*/

/*

Implementacion en C del algoritmo Needleman-Wunsch
Version 2. 26/11/2007.

*/

/*


*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <fcntl.h>
#include <sys/stat.h> 


#define FILENAME1 "file1.adn"
#define FILENAME2 "file2.adn"
#define FILENAME3 "file3.adn"

/*

Matriz de distancias entre nucleótidos:
- 	 A 	 G 	 C 	 T
A 	10 	-1 	-3 	-4
G 	-1 	 7 	-5 	-3
C 	-3 	-5 	 9 	 0
T 	-4 	-3 	 0 	 8
*/
static int s[4][4]={{10,-1,-3,-4},{-1,7,-5,-3},{-3,-5,9,0},{-4,-3,0,8}};

/*

Matriz de distancias entre nucleótidos:
- 	 A 	 G 	 C 	 T
A 	 1 	-1 	-1 	-1
G 	-1 	 1	-1 	-1
C 	-1 	-1 	 1 	-1
T 	-1 	-1 	-1 	 1
*/
//static int s[4][4]={{1,-1,-1,-1},{-1,1,-1,-1},{-1,-1,1,-1},{-1,-1,-1,1}};


static penalty = -5;//-2

//int m=8,n=11;//4,3

static int *a;//m+1
static int *b;//n+1

int DEBUG = 0;
int PRINT_RESULT = 0;

/*

  Pseudoce of the computeF function:

  for i=0 to length(A)-1
    F(i,0) <- penalty*i
  for j=0 to length(B)-1
    F(0,j) <- penalty*j
  for i=1 to length(A)
    for j = 1 to length(B)
    {
      Choice1 <- F(i-1,j-1) + S(A(i), B(j))
      Choice2 <- F(i-1, j) + penalty
      Choice3 <- F(i, j-1) + penalty
      F(i,j) <- max(Choice1, Choice2, Choice3)
    }

*/

int *computeF(int m, int n){
  int i = 0, j = 0;
  int choice1,choice2,choice3,max;
  int *f;

  if (DEBUG) printf("\n computeF with mxn dimensions: %dx%d",m,n);

  f = (int *)malloc(m*n*sizeof(int));

  if (f == NULL) {
    printf("\n ERROR computeF(): memoria insuficiente\n");
    exit(EXIT_FAILURE);
  }

  for(i=0;i<m;i++)
	*(f+i*n) = penalty*i;

  for(i=0;i<n;i++)
	*(f+i) = penalty*i; 

  for(i=1;i<m;i++) {
    for(j=1;j<n;j++) {
	choice1 = *(f+(i-1)*n+j-1)+s[a[i-1]][b[j-1]];
	choice2 = *(f+(i-1)*n+j)+penalty;
	choice3 = *(f+i*n+j-1)+penalty;

	if (choice1>choice2) 
	   max=choice1;
	else 
	   max=choice2;
	if (choice3>max) 
	   max=choice3;
	
	*(f+i*n+j) = max;
    } 
  }

  return f;
}


/*

   Pseudocode of the generation os sequence:

  AlignmentA <- ""
  AlignmentB <- ""
  i <- length(A) - 1
  j <- length(B) - 1
  while (i > 0 AND j > 0)
  {
    Score <- F(i,j)
    ScoreDiag <- F(i - 1, j - 1)
    ScoreUp <- F(i, j - 1)
    ScoreLeft <- F(i - 1, j)
    if (Score == ScoreDiag + S(A(i), B(j)))
    {
      AlignmentA <- A(i-1) + AlignmentA
      AlignmentB <- B(j-1) + AlignmentB
      i <- i - 1
      j <- j - 1
    }
    else if (Score == ScoreLeft + d)
    {
      AlignmentA <- A(i-1) + AlignmentA
      AlignmentB <- "-" + AlignmentB
      i <- i - 1
    }
    otherwise (Score == ScoreUp + d)
    {
      AlignmentA <- "-" + AlignmentA
      AlignmentB <- B(j-1) + AlignmentB
      j <- j - 1
    }
  }
  while (i > 0)
  {
    AlignmentA <- A(i-1) + AlignmentA
    AlignmentB <- "-" + AlignmentB
    i <- i - 1
  }
  while (j > 0)
  {
    AlignmentA <- "-" + AlignmentA
    AlignmentB <- B(j-1) + AlignmentB
    j <- j - 1
  }

*/


char getNucleotido(int i);

void genAlignment(int *f, int m, int n){
 int i,j,k;
 int Score,ScoreDiag,ScoreUp,ScoreLeft;


 char AlignmentA[m+n];
 char AlignmentB[n+m];
 int indexA=m+n-1, indexB=m+n-1; 

 i=m-1;
 j=n-1;

 while ((i>0)&&(j>0)) {
  Score = *(f+i*n+j);
  ScoreDiag = *(f+(i-1)*n+(j-1));
  ScoreUp = *(f+(i-1)*n+j);
  ScoreLeft = *(f+i*n+j-1);
  
  if (DEBUG) printf("\n Values of Score=%d ScoreDiag=%d ScoreUp=%d ScoreLeft=%d ",Score,ScoreDiag,ScoreUp,ScoreLeft);

  if (Score == (ScoreDiag + s[a[i-1]][b[j-1]]) ) {
    AlignmentA[indexA--] = getNucleotido(a[i-1]);
    AlignmentB[indexB--] = getNucleotido(b[j-1]);
    if (DEBUG) printf("\n %c_%c",getNucleotido(a[i-1]),getNucleotido(b[j-1]));
    i--;
    j--;
  }
  else if (Score == (ScoreLeft + penalty))
  {
    AlignmentA[indexA--] = getNucleotido(b[j-1]);
    AlignmentB[indexB--] = '-';   
    if (DEBUG) printf("\n %c_-",getNucleotido(b[j-1]));
    j--;
  }
  else if (Score == (ScoreUp + penalty)) 
  {
    AlignmentA[indexA--] = '-';
    AlignmentB[indexB--] = getNucleotido(a[i-1]);  
    if (DEBUG) printf("\n -_%c",getNucleotido(a[i-1]));   
    i--;
  }
 }//end while

  while (i > 0)
  {
    AlignmentA[indexA--] = '-';
    AlignmentB[indexB--] = getNucleotido(a[i-1]);
    if (DEBUG) printf("\n %c_' '",getNucleotido(a[i-1])); 
    i--;
  }
  while (j > 0)
  {
    AlignmentA[indexA--] = getNucleotido(b[j-1]);
    AlignmentB[indexB--] = '-';   
    if (DEBUG) printf("\n ' '-%c",getNucleotido(b[j-1]));
    j--;
  }


  int stepPrintf=70;
  if (DEBUG) {
  while (i<(m+n)){

  printf("\n A=");
  for(i=indexA+1;(i<(m+n))&&(i<(indexA+1+stepPrintf));i++){
   printf("%c",AlignmentA[i]);
  }

  printf("\n   ");
  for(i=indexA+1;(i<(m+n))&&(i<(indexA+1+stepPrintf));i++){
   if ((AlignmentA[i]=='-')||( AlignmentB[i]=='-'))
      printf(" ");
   else if (AlignmentA[i]==AlignmentB[i]) 
      printf("|");
   else
      printf(".");
  }


  printf("\n B=");
  for(i=indexB+1;(i<(m+n))&&(i<(indexB+1+stepPrintf));i++){
   printf("%c",AlignmentB[i]);
  }

  if (i<(m+n)) {
   indexA+=stepPrintf;
   indexB+=stepPrintf;
  }
  }//end-while
  }//end DEBUG


  if (PRINT_RESULT) 
	printf("\n\n\n Length=%d",(m+n-(indexA+1)));

  k=0;
  for(i=indexB+1;i<(m+n);i++){
     if (AlignmentB[i]==AlignmentA[i]) k++;
  }
  if (PRINT_RESULT)
  	printf("\n Similarity=%.2f%%  (%d/%d) \n",(double)(100*k)/(m+n-(indexA+1)),k,m+n-(indexA+1));


}



int convNucleotido(char c){
  switch (c) {
  case 'A': return 0;break;
  case 'G': return 1;break;
  case 'C': return 2;break;
  case 'T': return 3;break;
//  default: return -1;
  }

}

char getNucleotido(int i){
  switch (i) {
   case 0: return 'A';break;
   case 1: return 'G';break;
   case 2: return 'C';break;
   case 3: return 'T';break;
   default: return 'X';
  }
}



int sizeFile(char *fileName)
{
  FILE *fp;

  if ((fp = fopen(fileName, "r")) == NULL)
    printf("Error: File could not be opened\n");
  else
  {
    fseek(fp,0,SEEK_END);
    int file_size = ftell(fp);
    if (DEBUG) printf("\n %s: file size=%d", fileName, file_size);
    return file_size;
  }

}


char *readFile(char *fileName)
{
  FILE *fp;
  int DEBUG2 = 0;

  if ((fp = fopen(fileName, "r")) == NULL)
    printf("Error: File could not be opened\n");
  else
  {

    int file_size = sizeFile(fileName);;
    if (DEBUG2) printf(" file size=%d", file_size);

    char *cadena = (char *)malloc(file_size);

    int read_chars = 1;
    int read_lines = 0;
    while(read_chars<(file_size-read_lines)) {
        fgets(&(cadena[read_chars-1]),file_size,fp);
        if (cadena[read_chars-1]!='>') {
 	       read_chars=strlen(cadena);
	}
        read_lines++;
        if (DEBUG2) printf(" read_chars=%d", read_chars);
        if (DEBUG2) printf(" read_lines=%d", read_lines);
    }

    if (DEBUG2) {
        int i;
		char *cadena2 = cadena;
		for(i=1;i<=read_chars;i++){
		   printf("[%d,%c],",i,*cadena2++);
		}
    }

    return cadena;
   }

}


printF(int *f, int m, int n){
	int i,j,value;

	printf("\n\n Matriz F \n-------------------------\n        ");
	for(j=0;j<n;j++) {
		printf("%c,  ",getNucleotido(b[j]));
	}
	printf("\n");

	for(i=0;i<=m;i++){
		if (i>0) printf("%c ",getNucleotido(a[i-1]));
		else printf("  ");
		for(j=0;j<=n;j++) {
			value = *(f+i*(n+1)+j);
			printf("%3d,",value);
		}
		printf("\n");
	}
	printf("-------------------------\n");
}


int main(int argc, char **argv){

	int i,j,value;
	char *cad1, *cad2;
	int m,n;

	int myrank,
		numprocs,
		itask;
	double mytime,   /*variables used for gathering timing statistics*/
		   maxtime,
		   mintime,
		   avgtime;

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	 
	int sendcounts[numprocs];
	int displs[numprocs];

	PRINT_RESULT = 1;
	DEBUG=1;

//  For debugging:
//char cadena1[] = "CGAGACGT";   //2,1,0,1,0,2,1,3,3,
//char cadena2[] = "AGACTAGTTAC";//0,1,0,2,3,0,1,3,3,0,2,2,0

  
	if (argc<3) {
		printf("\n Usage: ./cadenas-adn filey.adn filex.adn \n\n");
		exit(EXIT_FAILURE);
	}

	if (!myrank) {
	  /* Realizar toda la entrada/salida de fichero desde el procesador 0 */
		printf("\n ---------------- Begin of Needleman-Wunsch implementation -------------------------\n");  
	
	

	//  char *cad1=readFile(FILENAME1); 
	//  char *cad2=readFile(FILENAME2);

		cad1=readFile(argv[1]);
		cad2=readFile(argv[2]);

		m = strlen(cad1)-1;  
		n = strlen(cad2)-1;

//TODO: Lo de obtener las matrices a enteros no lo podría hacer cada proceso?
	
		a = malloc(sizeof(int) * m);
		b = malloc(sizeof(int) * n);

		DEBUG=0;

		if (DEBUG) printf("\n a[]=");
		for(i=0;i<m;i++){
			a[i] = convNucleotido(cad1[i]);
			if (DEBUG) printf("%d,",a[i]);
		}
		if (DEBUG) printf("\n");


		if (DEBUG) printf("\n b[]=");
		for(j=0;j<n;j++){
			b[j] = convNucleotido(cad2[j]);
			if (DEBUG) printf("%d,",b[j]);
		}
		if (DEBUG) printf("\n");

	
	}


	int errorcodeM = MPI_Bcast (&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (errorcodeM != MPI_SUCCESS)
		MPI_Abort(MPI_COMM_WORLD,errorcodeM);
	
	int errorcodeN = MPI_Bcast (&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (errorcodeN != MPI_SUCCESS)
		MPI_Abort(MPI_COMM_WORLD,errorcodeN);

	if(myrank)
		a = malloc(sizeof(int) * m);
	int errorcode = MPI_Bcast (a, m, MPI_INT, 0, MPI_COMM_WORLD);
	if (errorcode != MPI_SUCCESS)
		MPI_Abort(MPI_COMM_WORLD,errorcode);


	
	displs[0] = 0;
	sendcounts[0] = (n / numprocs) + ((0 < (n%numprocs))?1:0);
	for(itask=1; itask<numprocs; itask++){
		sendcounts[itask] = (n / numprocs) + ((itask < (n%numprocs))?1:0);
		displs[itask]=displs[itask-1] + sendcounts[itask-1];
	}

	if(!myrank){
		printf("\n");
		for(i=0;i<numprocs;i++){
			printf("Task %d: Vector sent to %d: sum=%d size= %d\n",
		    myrank,i,displs[i],sendcounts[i]);//fflush(stdout);
		}
	}

	int nlocal = (n / numprocs) + ((myrank < (n%numprocs))?1:0);

	/**TODO
	DUDA (el proc 0 tambien hace trabajo ¿?)
	n = 10
	numprocs = 4
	-------------
	nlocal(0) = 2.5 + (0<2) = 3.5 -> 3
	nlocal(1) = 2.5 + (1<2) = 3.5 -> 3
	nlocal(2) = 2.5 + (2<2) = 2.5 -> 2
	nlocal(3) = 2.5 + (3<2) = 2.5 -> 2
	*/
	
	// Create a buffer that will hold a subset of the random numbers
	int *sub_b = malloc(sizeof(int) * nlocal);
	// Scatter b to all processes
	MPI_Scatterv(b, sendcounts, displs, MPI_INT, sub_b, nlocal, MPI_INT, 0, MPI_COMM_WORLD);

//DEBUG=1;
	mytime = MPI_Wtime();  /*get time just before work section */
	int *f = computeF(m+1,nlocal+1);
 
	if (DEBUG) printF(f,m,nlocal);

	genAlignment(f,m+1,nlocal+1);

	mytime = MPI_Wtime() - mytime;  /*get time just after work section*/

/* Reduce and print results - Insert here those lines */

	MPI_Finalize();
	
	if (!myrank) {
	/* Realizar toda la entrada/salida de fichero desde el procesador 0 */
		printf("\n ---------------- End of Needleman-Wunsch implementation -------------------------\n");
	}
	
	exit(EXIT_SUCCESS);
}
