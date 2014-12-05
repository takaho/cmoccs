cmoccs
======

MOCCS implementation in C++

INSTALL
==

tar xfvz moccs-X.Y.Z.tar.gz
./configure && make && make install

How to use
==

moccs [OPTIONS]

OPTIONS
==
	-i <filename>
	   input filename, multiple fasta files containing sequences around peaks

	-o <filename>
	   output filename

	-m <motif unit>
	   motif size

	-e <expected>
	   cutoff value of occurrence, motifs observed more than X times than expected will be displayed

	-n <max_number>
	   maximum number of displayed motifs

	-w <window size>
	   length of sequences around peak
	   
	--strand
	    strand-specific search
	
	--verbose
	    verbose mode
