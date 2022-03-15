# BE_Tiling_Library

Generate sgRNA tiling library for cytosine base editor and adenine base editor.
Input is a single gene symbol or list of gene symbol.


required file
1. make directory named hg19 and hg38 in EssentialData
2. then download fasta file of GRCh37 and GRCh38 divided by chromosome. 
3. rename each fasta file (eg hg19_1.fasta)


example
1. you can add interest gene symbol in Input/input.txt, don't delete first row 'gene'
2. then Execute BE_SGE.py, if you want find NG PAM library, you can change pamtype NGG into NG
