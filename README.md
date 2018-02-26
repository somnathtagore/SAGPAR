# SAGPAR

*************************************************************************************************************************************************************
                                       SAGPAR -> StructurAl Grammar-based automated PAthway Reconstruction
                                                HELP FILE AND MANUAL FOR RUNNING THE PROGRAM
*************************************************************************************************************************************************************

Contents
--------

1. Introduction
2. Platform Requirements
3. Hardware Requirements
4. Software Requirements
5. Usage
6. Input file format
7. Output
8. Sample file
9. Rejoinder files 
10. Interpretation
11. Default messages
----------------------------------------------------------------------------------------------------------------------------------------------
######################################################################
Introduction
######################################################################

SAGPAR is an automated pathway reconstruction algorithm, capable of 
reconstructing a metabolic pathway from a given set of randomly chosen 
metabolites. It finds links among metabolites based upon some predefined 
criteria, conditions and thresholds. It can even predict multi-links 
among metabolites. 

######################################################################
Platform Requirements
######################################################################

The program can run on any Unix and Windows platform. 

######################################################################
Hardware Requirements
######################################################################

An Intel Pentium 4, 2.56 GHz processor, 256 MB RAM and sufficient virtual
memory is the minimum requirement.

######################################################################
Software Requirements
######################################################################

The program can run on any Unix and Windows platform. For running in Unix
platform, a 'gcc' package should be installed, whereas, in Windows platform,
any typical C/C++ compiler can be installed. For example, Turbo C++, Borland
C++ etc.

######################################################################
Usage
######################################################################

    ==== How to Compile ====

****** Unix ******
gcc -ld SAGPAR.c

****** Windows ******
Save > Compile

    ==== How to Run ====

****** Unix ******
./a.out

Enter the input file:(Note: File name should be less than 100 characters)...

****** Windows ******

> Run

Enter the input file:(Note: File name should be less than 100 characters)...

Example)
Enter the input file:(Note: File name should be less than 100 characters)...
newsmilesPent.txt

######################################################################
Input File Format
######################################################################

Note: The input file shoul be a "text file", containing SMILES string format of the
input metabolites. Each metabolite's SMILES should be in a new line.

Example:
C(C1C(C(C(C(O)O)O)O)O)O
C(C(C(O)O)O
C(C1C(C(C(O)O)O)O)O

######################################################################
Output
######################################################################

****** Unix ******

Note: There are 2 ways in which output can be viewed. 
1. Using the terminal itself.
2. Using a separate output file.
In 2nd case, the file name should be provide. See example below.

Example:
>>> 1. Using terminal
./a.out
Enter the input file:(Note: File name should be less than 100 characters)...
newsmilesPent.txt

>>> 2. Using output file
./a.out > OutputFile.txt
Enter the input file:(Note: File name should be less than 100 characters)...
newsmilesPent.txt

****** Windows ******

Note: There are 2 ways in which output can be viewed. 
1. Using the command prompt itself.
2. Using a separate output file.
In 2nd case, the file name should be provide. See example below.

Example:
>>> 1. Using command prompt
>Run
Enter the input file:(Note: File name should be less than 100 characters)...
newsmilesPent.txt

>>> 2. Using output file
>Run OutputFile.txt
Enter the input file:(Note: File name should be less than 100 characters)...
newsmilesPent.txt

######################################################################
Sample file
######################################################################

A sample file named "newsmilesPent.txt" is provided for the sake of
convenience.

######################################################################
Rejoinder files
######################################################################

1. A rejoinder file called "Pat.txt" would be created upon running the
program. Users can check the "patterns" generated from the SMILES in it.
2. A rejoinder file called "Rbit.txt" would be created upon running the program.
Users can check the "binary strings" produced from the patterns in it for 
cross-verification.

######################################################################
Interpretation
######################################################################

The final output consists of the best possible scores among the metabolite
pairs. The better the score, the more probable the pairing is. The links are
made by identifying these best scores among the metabolites (identified in 
terms of i & j values) produced by the program. These i & j values correspond
to the pairing of ith & jth metabolites present in the input file.

######################################################################
Default messages
######################################################################

Note: Some default messages can be thrown after the program is executed. These are
not errors and can be easily avoided. These messages are,
1. "Segmentation fault" - To debug this message, users need to increase the
size of arrays "rpat", "breact" and "bprod" according to their need of storing 
the values.
2. "File not found" - TO debug this message, re-confirm the input file name 
provided.

######################################################################
END OF MANUAL
######################################################################
