/*
Program Name	:- SAGPAR.c
Program Topic	:- StructurAl Grammar-based automated PAthway Reconstruction algorithm for automated reconstruction of any metabolic pathway from a given set of randomly chosen metabolites.
*/

//Importing header files
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

//Global variable declaration
char react[40][80];
int rmass[40];
float ss[40][40],mscore[40][40],fscore[40][40];
int NO;

//Calculating best possible scores
void largest()
{
	float large=0,seclarge=0;
	int i,j,m,n;
	
	printf("*****************************************************************************");
	printf("\n Some best possible scores among metabolites\n"); 
	printf("Note: i & j values represent the metabolite-pairs, that are formed from metabolites present in the input file\n"); 
	for(i=0;i<NO;i++)
	{
		m=0,n=0;
		large=fscore[i][0];
		seclarge=fscore[i][0];
		
		for(j=0;j<NO;j++)
		{
			if(fscore[i][j]>large)
			{
				large=fscore[i][j];
				m=j;
			}
		}
		for(j=0;j<NO;j++)
		{
			if((fscore[i][j]>seclarge) && (j!=m))
			{
				seclarge=fscore[i][j];
				n=j;
			}
		}
		printf("Row[%d] :- (i=%d,j=%d)BEST = %f\t (i=%d,j=%d)NEXTBEST = %f\n",i,i,m,large,i,n,seclarge);
	}
}

//Function finalscore() for calculating the final score for every metabolite pair
void finalscore()
{
	//Declaring local variables
	int i,j;

	//Loop for calculating final score
	for(i=0;i<NO;i++)
	{
		for(j=0;j<NO;j++)
		{
			//Score calculation and storing them in an array
			fscore[i][j]=ss[i][j]*mscore[i][j];
		}
	}
	
	//Printing the final score
	for(i=0;i<NO;i++)
	{
		for(j=0;j<NO;j++)
		{
			printf("fscore[%d][%d] = %f\t",i,j,fscore[i][j]);
		}
		printf("\n");
	}
}//Function finalscore() ends

//Function simscore() for calculation of similarity scores
void simscore(int no)
{
	//File handles
	FILE *fp,*fr;

	//Local variables declaration
	int i,r,a=0,b=0,ca,d,n,k,dtemp,blen,plen,cloop;

	char breact[245000]={'\0'},bprod[245000]={'\0'};
	float sscore;

	//Opening files in read mode
	if((fr=fopen("Rbit.txt","r"))==NULL)
	{
		printf("Cannot open file\n");
		exit(0);
	}
	for (r=0;r<no;r++)
	{
		fscanf(fr,"%s",breact);

	//Parsing the file and extracting the required output
	//Calculating the length of the string array
		blen=strlen(breact);
	//Calculating 'a'
		for(k=0;k<blen;k++)
		{
			if(breact[k]=='1') a++;
		}
	
		if((fp=fopen("Rbit.txt","r"))==NULL)
		{
			printf("Cannot open file\n");
			exit(0);
		}


	//Loop for comparison among two metabolites 
		for(i=0;i<no;i++)
		{
		//Local variable declaration
			b=0;
			ca=0;
			dtemp=0;
		
		//while loop begins
			fscanf(fp,"%s",bprod);
		//Calculating the string length
			plen=strlen(bprod);

		//Calculation of 'b'
			for(k=0;k<plen;k++)
			{
				if(bprod[k]=='1') b++;
			}
		
		//Calculation of 'n'
			n=blen<plen?plen:blen;
			cloop=blen>plen?plen:blen;
		
		//Calculation of 'ca'
			for(k=0;k<cloop;k++)
			{
				if((breact[k]=='1') && (bprod[k]=='1'))
				{
					dtemp++;
					ca++;
				}
			}
		
		//Calculation of 'd'
			for(k=cloop;k<n;k++)
			{
				if(blen>plen)
				{
					if(breact[k]=='1') dtemp++;
				}
				else if(bprod[k]=='1') 	dtemp++;
			}
			d=n-dtemp;

		//Calculation of final sscore
			sscore=sqrt((ca+d)*1./n);

		//Storing sscore in an array
			ss[r][i]=sscore;
		}
		fclose(fp);
	}	
	fclose(fr);
	//Closing file

}//Function simscore() ends

//Function massscore() for calculation of masscore begins
void masscore(int no)
{
	//Declaration of local variables
	int i,k,j,rlen;
	float mc,md,me,mf,car=12,oxy=16,hyd=1,pho=30.9,nit=14,sul=32.1,resmass=0;

	//Loop starts
	for(i=0;i<no;i++)
	{
		for(k=0;k<no;k++)
		{
			mc=md=me=mf=0;

			//String comparision....For side reactuct determination
			if(strcmp(react[i],"C(C1C(C(C(C(O1)O)O)O)O)O")==0 && strcmp(react[k],"C(C1C(C(C(C(O1)O)O)O)O)OP(=O)(O)O")==0)
			{
				mc=7*car+12*hyd+4*nit+5*oxy+1*pho;
				md=7*car+10*hyd+4*nit+2*oxy;
			}
			else if(strcmp(react[i],"C(C1C(C(C(O1)(CO)O)O)O)OP(=O)(O)O")==0 && strcmp(react[k],"C(C1C(C(C(O1)(COP(=O)(O)O)O)O)O)OP(=O)(O)O")==0)
			{
				mc=10*car+16*hyd+5*nit+13*oxy+3*pho;
				md=10*car+15*hyd+5*nit+10*oxy+2*pho;
			}
			else if(strcmp(react[i],"C(C1C(C(C(O1)(COP(=O)(O)O)O)O)O)OP(=O)(O)O")==0 && strcmp(react[k],"C(C(=O)COP(=O)(O)O)O")==0)
			{      
				mc=3*car+7*hyd+6*oxy+pho;
			}
			else if(strcmp(react[i],"C(C(C=O)O)OP(=O)(O)O")==0 && strcmp(react[k],"C(OP(O)(O)=O)(=O)C(O)COP(O)(O)=O")==0)
			{
				mc=21*car+28*hyd+7*nit+14*oxy+2*pho;
				md=3*hyd+4*oxy+pho;
				me=21*car+29*hyd+7*nit+14*oxy+2*pho;
				mf=hyd;
			}
			else if(strcmp(react[i],"C(OP(O)(O)=O)(=O)C(O)COP(O)(O)=O")==0 && strcmp(react[k],"C(C(C(=O)O)O)OP(=O)(O)O")==0)
			{
				mc=10*car+16*hyd+5*nit+13*oxy+3*pho;
				md=10*car+15*hyd+5*nit+10*oxy+2*pho;
			}
			else if(strcmp(react[i],"C(C(C(=O)O)OP(=O)(O)O)O")==0 && strcmp(react[k],"C=C(C(=O)O)OP(=O)(O)O")==0)
			{
				mc=2*hyd+oxy;
			}
			else if(strcmp(react[i],"C=C(C(=O)O)OP(=O)(O)O")==0 && strcmp(react[k],"CC(=O)C(=O)[O-]")==0)
			{
				mc=10*car+16*hyd+5*nit+13*oxy+3*pho;
				md=10*car+15*hyd+5*nit+10*oxy+2*pho;
			}

			//String the mass scores in an array
			mscore[i][k]=(rmass[i]+rmass[k])/(rmass[i]+rmass[k]+mc+md+me+mf);
		}
	}
}//Function masscore() ends

//Function mass() for calculation of metabolite mass begins
void mass(int i)
{
	//Local variable declaration
	int k,rlen;
	float resmass=0;

	//String length calculation
	rlen=strlen(react[i]);
	resmass=0;

	//Calculation of mass
	for(k=0;k<rlen;k++)
	{
		if(react[i][k]=='C') resmass+=12;
		else if(react[i][k]=='O') resmass+=16;
		else if(react[i][k]=='H') resmass+=1;
		else if(react[i][k]=='P') resmass+=30.9;
		else if(react[i][k]=='N') resmass+=14;
		else if(react[i][k]=='S') resmass+=32.1;
		else if(react[i][k]=='(' || react[i][k]==')' || react[i][k]=='1' || react[i][k]=='=') resmass+=0;
	}

	//Storing the mass in an array
	rmass[i]=resmass;
}//Function mass() ends

//Function searchsmile() for storing the smiles string in arrays begins
int searchsmile()
{
	//Local variable declaration
	char name1[60],name2[60],inpfile[100];
	int i=0;

	//File Handle
	FILE *f1;
	
	printf("Enter the input file: (Note: File name should be less than 100 characters)...");
	scanf("%s",&inpfile);
	//Opening file
	//if((f1=fopen("newsmilesPent.txt","r"))==NULL)
	if((f1=fopen(inpfile,"r"))==NULL)
	{
		printf("Cannot open file\n");
		exit(0);
	}

	//Parsing the file and string the data in arrays
	while(fscanf(f1,"%s",react[i])==1)
	{
		i++;
	}

	//Closing the file
	fclose(f1);
	return i;
}//Function searchsmile() ends 

//Function pattern() for generating a list of patterns
void pattern(int no)
{
	//Declaring the file handles
	FILE *fp,*fp1,*fp2;
	
	//Declaring the local variables

	char rpat[11050][90];
	int i,j,c,i1,j1,c1,n=strlen(react[no]);

	int rptlen,oldlen,k,k1,oldlen1; 
	char rptstr[6]; 
	
	c=0;

	fp2=fopen("Pat.txt","a");
	//Loop begins for generating patterns 
	for(i=0;i<n;i++)
	{
		for(j=0;j<n-i;j++)
		{
			strncpy(rpat[c],react[no]+j,i+1);
			rpat[c][i+1]='\0';
			fprintf(fp2,"%s\n",rpat[c]);
			c++;
		}
	}
	
	//Opening a file in append mode
	fp=fopen("Rbit.txt","a");

	oldlen=0;

	//Replacing the patterns by 1,0
	for(i=0;i<c;i++)
	{
		rptlen=strlen(rpat[i]);
		for(k=0;k<rptlen;k++)
		{
			strcpy(rptstr,"");
			if(rpat[i][k]=='C') strcpy(rptstr,"111111");
			else if(rpat[i][k]=='O') strcpy(rptstr,"111110");
			else if(rpat[i][k]=='H') strcpy(rptstr,"111100");
			else if(rpat[i][k]=='P') strcpy(rptstr,"111000");
			else if(rpat[i][k]=='N') strcpy(rptstr,"110000");
			else if(rpat[i][k]=='S') strcpy(rptstr,"100000");
			else if(rpat[i][k]=='(' || rpat[i][k]==')' || rpat[i][k]=='1' || rpat[i][k]=='=') strcpy(rptstr,"");
			fprintf(fp,"%s",rptstr);
		}
	}
	fprintf(fp,"\n");

	//Closing the file handle
	fclose(fp);

	fflush(stdin);
}//Function pattern() ends

//main() begins
main()
{
	//Declaring local variables
	int j;
	char initial[60],target[60],smi1[90],smi2[90];

	//File handles
	FILE *fp,*fp1,*fpr;
//	clrscr();
	
	//Function call
	NO=searchsmile();
	
	fp=fopen("Rbit.txt","w"); 
	fclose(fp);
	for(j=0;j<NO;j++)
		pattern(j);
	
	for(j=0;j<NO;j++)
		mass(j);
	masscore(NO);
//	getch();
	simscore(NO);

	finalscore();
	largest();
	exit(1);
}
//Program Ends
