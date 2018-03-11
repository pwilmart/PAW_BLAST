# PAW_BLAST
A utility for blasting one protein FASTA file (queries) against another (targets) to find orthologs.

**What is this and why would I use it:**
This compares protein sequences from one FASTA protein database against another. It determines the reciprocal best matches (a basic ortholog definition). It uses a downloaded local installation of the NCBI BLAST program. A tab-delimited text file of the matches is written to disk.

There are many uses of such a tool. The concept is to take a list of identified proteins from a proteomics experiment and create a subset FASTA file with just those sequences. The second FASTA database can be several things: identifications from another proteomics experiment (example: compare mouse eye lens proteins to human eye lens proteins), a better annotated related proteome (example: Leshmainia donovani IDs versus Leshmania infantum reference proteome), a reference proteome for a high quality model system (example: rat identifications versus a mouse Swiss-Prot reference proteome), a set of sequences from another publication (assuming that protein sequences (not accessions) can be obtained).

The text summary file includes both the accessions and descriptions from the target database, so this is an easy way to fetch protein descriptions from another source for poorly annotated protein databases. Some model systems have impressive annotations and others do not. Rat is poor compared to mouse. Monkey is poor compared to human. Rat and monkey should have considerable sequence similarity to mouse and human, respectively, so annotations from orthologs should be very informative.

There is some computational overhead to running BLAST, hence the idea of using just a smaller list of identified proteins as the query FASTA database. Generally, a few thousand query sequences against 20K references sequences seems to work fine. There is a Word file with more details on running the program and interpreting the output.

Creating the query FASTA database is up to the user. The formats of lists of protein identifications vary by search engine and there is the usual issue of protein groups. Generally, proteins that get grouped in most protein inference algorithms will have sequence similarity and picking one protein sequence to represent the group will be fine. Groups can also be expanded, but it is a little safer to use one representative. This is a reciprocal best match so lower protein redundancy tends to work better (both for the query database and the target database).
