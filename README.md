# cycle_finder README


## AUTHOR
Yoshiki Tanaka wrote the original code.


## VERSION
1.0.0


## DESCRIPTION
cycle_finder is a tool for detecting tandem repeat and interspersed repeat from short reads. The execusion procedures are as follows:
1. Extracting high frequency k-mer from short reads.
2. Finding cycle from de Bruijn graph and detect tandem repeats.
3. Clustering and copy estimating from detected tandem repeats.
4. Finding path from de Bruijn graph that removed cycle and detect interspersed repeats.
5. Clustering and copy estimating from detected interspersed repeats.


## INSTALATION
### Using Bioconda (Linux)
```bash
conda install -c bioconda -c conda-forge cycle_finder
```
### From source
```bash
tar zxfv cycle_finder_<version>.tar.gz
cd cycle_finder_<version>
make
cp cycle_finder <installation_path>
```


## SYNOPSIS
### single mode
```bash
cycle_finder all -f <SHORT_READS>.fastq
```

### compare mode
```bash
cycle_finder all -f1 <SHORT_READS1.fastq> -f2 <SHORT_READS2.fastq>
```


## TEST
```bash
# Get the test dataset, which includes simulated reads from repeats-inserted E. coli genomes.  
wget https://github.com/rkajitani/cycle_finder/releases/download/v1.0.0/test_data.tar.gz
tar xzfv test_data.tar.gz

# Test for tandem repeats.
# The reads, genome, tandem repeat unit (253 bp; mutation rate, 2%), and correct result are
# reads.fq, genome.fa, rep_unit.fa, and result/*, respectively.
cd test_data/tandem/
bash cmd.sh
# Output: out_T.fa out_T.tsv

cd ../..

# Test for interspersed repeats.
# The reads, genome, interspersed repeat unit (253 bp; mutation rate, 2%), and correct result are
# reads.fq, genome.fa, rep_unit.fa, and result/*, respectively.
cd test_data/interspersed/
bash cmd.sh
# Output: out_I.fa out_I.tsv
```


## DEPENDENCY
- GCC  

- OpenMP  

- Jellyfish (>= 2.2.6)  
G. Marçais and C. Kingsford, “A fast, lock-free approach for efficien
parallel counting of occurrences of k-mers,” Bioinformatics, vol. 27
no. 6, pp. 764–770, 2011.
https://academic.oup.com/bioinformatics/article/27/6/764/234905

- TRF (Tandem Repeat Finder; >= 4.07b)  
G. Benson, “Tandem repeats finder : a program to analyze DNA sequence,”
vol. 27, no. 2, pp. 573–580, 1999.

- BLAST+ (>= 2.2.31+)  
S. F. Altschul, T. L. Madden, A. A. Schäffer, J. Zhang, Z. Zhang, W. Miller,
and D. J. Lipman, “Gapped BLAST and PSI-BLAST: A new generation of proten
database search programs,” Nucleic Acids Res., vol. 25, no. 17
pp. 3389–3402, 1997.

- CD-HIT (>= 4.6)  
W. Li, L. Fu, B. Niu, S. Wu, and J. Wooley, “Ultrafast clustering algorithm
for metagenomic sequence analysis,” Brief. Bioinform., vol. 13, no. 6, pp
656–668, 2012.


## USAGE

### COMMON OPTIONS
```
-t INT    : Number of threads (default 1)

-o STR    : Prefix of output files (default out)
```


### cycle_finder extract [OPTIONS]
Extracting high frequency k-mer from short reads.

#### INPUT OPTIONS
    
##### single mode
```
-f FILE1 [FILE2 ...]: Reads file (fasta or fastq format)
```

##### compare mode
```
-f1 FILE1 [FILE2 ...]: Reads file (fasta or fastq format)
-f2 FILE1 [FILE2 ...]: Reads file (fasta or fastq format)
```

    
#### OTHER OPTIONS
```
-c1 INT  : copy difference for detecting repeat

-c2 INT  : copy difference for copy estimation
```

##### single mode
```
-jf FILE : if you have jf file that output by jellyfish,
            you can skip jellyfish count.

-dm FILE : if you have dump file that output by jellyfish,
            you can skip jellyfish dump.

-hs FILE : if you have histo file that output by jellyfish,
            you can skip jellyfish histo.
```

##### compare mode
```
-jf1 FILE : if you have jf file that output by jellyfish, 
            you can skip jellyfish count(corresponds to FILE1).

-jf2 FILE : if you have jf file that output by jellyfish,
            you can skip jellyfish count(corresponds to FILE2).

-dm1 FILE : if you have dump file that output by jellyfish, 
            you can skip jellyfish dump(corresponds to FILE1).

-dm2 FILE : if you have dump file that output by jellyfish, 
            you can skip jellyfish dump(corresponds to FILE2).

-hs1 FILE : if you have histo file that output by jellyfish,
            you can skip jellyfish histo(corresponds to FILE1).

-hs2 FILE : if you have histo file that output by jellyfish,
            you can skip jellyfish histo(corresponds to FILE2).
```


#### OUTPUT FILES
```
PREFIX_for_detect${c1}.fa           : high frequency k-mer for detecting repeat

PREFIX_for_estimate${c2}.fa         : high frequency k-mer for copy estimating

PREFIX_kmer_peak                    : k-mer peak

PREFIX.jf /PREFIX.jf1 / PREFIX.jf2  : jf file output by jellyfish

PREFIX.dm / PREFIX.dm1 / PREFIX.dm2 : dump file output by jellyfish

PREFIX.hs / PREFIX.hs1 / PREFIX.hs2 : histogram of k-mer frequency
```


### cycle_finder cycle [OPTIONS]
Finding cycle from de Bruijn graph and detect tandem repeats.

#### INPUT OPTIONS
```
-f FILE   : high frequency k-mer for detecting repeat

-r FASTQ  : short read fastq file(only one file)

-c INT    : single mode->1 compare mode->2

-l INT    : threshold of detecting repeat length

-n INT    : threshold of searching nodes

-d INT    : threshold of searching depth

-p INT    : k-mer peak (single mode)

-p1 INT   : k-mer peak (compare mode corresponds to FILE1)

-p2 INT   : k-mer peak (compare mode corresponds to FILE2)

-rc FLOAT : Down sampling coverage (default 0.5)
```

#### OUTPUT FILES
```
PREFIX_T_repeat_num     : temporal copy estimation

PREFIX_T_repeat_num_min : temporal minimum copy estimation
```


## cycle_finder cluster [OPTIONS]
Clustering and copy estimating from detected repeat.

### INPUT OPTIONS
```
-f1 FILE : high frequency k-mer for copy estimating

-f2 FILE : temporal copy estimation

-f3 FILE : temporal minimum copy estimation

-c  INT  : single mode->1 compare mode->2

-p INT   : k-mer peak (single mode)

-p1 INT  : k-mer peak (compare mode corresponds to FILE1)

-p2 INT  : k-mer peak (compare mode corresponds to FILE2)

-m  INT  : threshold mismatch of kmer alignment for copy estimation (default 0)
```

### OUTPUT FILES
#### tandem repeat
```
PREFIX_T_blst.blastn : k-mer alignment file

PREFIX_T.fa          : detected tandem repeat

PREFIX_T.tsv         : copy esitimation

PREFIX_T_clst.fa     : fasta file of detected tandem repeats
```

#### interspersed repeat
```
PREFIX_I.fa          : detected interspersed repeat

PREFIX_I.tsv         : copy estimation

PREFIX_I_clst.fa     : fasta file of detected tandem repeats
```


## cycle_finder intersperse [OPTIONS]
Finding path from de Bruijn graph that removed cycle and detect interspersed 
repeats.

### INPUT OPTIONS
```
-f FILE   : high frequency k-mer for detecting repeat

-b FILE   : alignment file (output from cycle command)

-r FASTQ  : short read fastq file(only one file)

-c INT    : single mode->1 compare mode->2

-L INT    : threshold of detecting repeat length

-N INT    : threshold of searching nodes

-D INT    : threshold of searching depth

-p INT    : k-mer peak (single mode)

-p1 INT   : k-mer peak (compare mode corresponds to FILE1)

-p2 INT   : k-mer peak (compare mode corresponds to FILE2)

-rc FLOAT : Down sampling coverage (default 0.5)
```


### OUTPUT FILES
```
PREFIX_for_detect${c1}.fa_no_cycle : high frequency k-mer for detect without cycle

PREFIX_I_repeat_num                          : temporal copy estimation

PREFIX_I_repeat_num_min                      : temporal minimum copy estimation
```


## cycle_finder all [OPTIONS]
Run whole pipeline:  
extract -> cycle -> cluster -> intersperse > cluster

### INPUT OPTIONS
#### single mode
```
-f FILE1 [FILE2 ...]: Reads file (fasta or fastq format)
```

#### compare mode
```
-f1 FILE1 [FILE2 ...]: Reads file (fasta or fastq format)
-f2 FILE1 [FILE2 ...]: Reads file (fasta or fastq format)
```

### OTHER OPTIONS
```
-c1 INT  : copy difference for detecting repeat

-c2 INT  : copy difference for copy estimation
```

#### single mode
```
-jf FILE : if you have jf file that output by jellyfish,
            you can skip jellyfish count.

-dm FILE : if you have dump file that output by jellyfish,
            you can skip jellyfish dump.

-hs FILE : if you have histo file that output by jellyfish,
            you can skip jellyfish histo.
```

#### compare mode
```
-jf1 FILE : if you have jf file that output by jellyfish, 
            you can skip jellyfish count(corresponds to FILE1).

-jf2 FILE : if you have jf file that output by jellyfish,
            you can skip jellyfish count(corresponds to FILE2).

-dm1 FILE : if you have dump file that output by jellyfish, 
            you can skip jellyfish dump(corresponds to FILE1).

-dm2 FILE : if you have dump file that output by jellyfish, 
            you can skip jellyfish dump(corresponds to FILE2).

-hs1 FILE : if you have histo file that output by jellyfish, 
            you can skip jellyfish histo(corresponds to FILE1).

-hs2 FILE : if you have histo file that output by jellyfish, 
            you can skip jellyfish histo(corresponds to FILE2).

-l INT    : threshold of detecting repeat length

-n INT    : threshold of searching nodes

-d INT    : threshold of searching depth

-L INT    : threshold of detecting repeat length

-N INT    : threshold of searching nodes

-D INT    : threshold of searching depth

-rc FLOAT : Down sampling coverage (default 0.5)

-m  INT   : threshold mismatch of kmer alignment for copy estimation (default 0)
```

### OUTPUT FILES
```
PREFIX_for_detect${c1}.fa          : high frequency k-mer for detecting repeat

PREFIX_for_estimate${c2}.fa          : high frequency k-mer for copy estimating

PREFIX_kmer_peak                          : k-mer peak

PREFIX.jf(PREFIX.jf1, PREFIX.jf2)   : jf file output by jellyfish

PREFIX.dm / PREFIX.dm1 / PREFIX.dm2 : dump file output by jellyfish

PREFIX.hs / PREFIX.hs1 / PREFIX.hs2 : histogram of k-mer frequency

[tandem repeat]			

PREFIX_T_repeat_num     : temporal copy estimation

PREFIX_T_repeat_num_min : temporal minimum copy estimation

PREFIX_T.fa          : detected tandem repeat

PREFIX_T.tsv         : copy esitimation

[interspersed repeat]

PREFIX_I_repeat_num     : temporal copy estimation

PREFIX_I_repeat_num_min : temporal minimum copy estimation

PREFIX_I.fa          : detected interspersed repeat

PREFIX_I.tsv         : copy estimation
```


## FILE FORMAT
### PREFIX_T_repeat_num / PREFIX_T_repeat_num_min / PREFIX_I_repeat_num / PREFIX_I_repeat_num_min

- column1  : sequence

- column2  : sequence length

- column3  : explanation

- column4  : difference between column5 and column4

- column5  : the copy number (FILE or FILE1)

- column6  : the copy number (FILE2 ; in single mode, this column shows 0)

- column7  : explanation

- column8  : difference between column9 and column10

- column9  : the number of bases (FILE or FILE1)

- column10 : the number of bases (FILE2 ; in single mode, this column shows 0)


### PREFIX_T.tsv / PREFIX_I.tsv

- column1 : "Family" name (the number next to '_' shows "Family number", length,
            sequence number from left)

- column2 : difference between column3 and column4

- column3 : the number of bases consists of the "Family" (FILE or FILE1)

- column4 : the number of bases consists of the "Family" (FILE2 ; in single mode,
            this column shows 0)

- column5 : difference between normalized column3 and normalized column4
    	      (normalized means devided by k-mer peak, so it shows difference 
	      of the copy number)


### PREFIX_T_clst.fa / PREFIX_I_clst.fa
fasta file of detected repeats. The sequence which has "*" is a representative
sequence of "family".


### PREFIX_T_blst.blastn
- column1 : sequence name

- column2 : k-mer aligned to the repeat

- column3 : the number of match
    
### PREFIX_T.element / PREFIX_I.element
after ">" shows Family, and other lines shows the "family number".
