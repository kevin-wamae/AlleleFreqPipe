# R-Toolkit_FreqSNP
An R pipeline for extracting single nucleotide polymorphisms (SNPs) from DNA 
sequence alignments and computing their frequencies.

The script takes as input a DNA sequence alignment in Phylip format and generates 
a frequency table of all SNPs and as well as a list of the polymorphic loci.

I am still trying to work out how to present the output "example_snpFreq.csv" (shown below),
so as to have one row for bi-allelic SNPs, two rows for tri-allelic SNPs and so on,
contrary to what is seen e.g. with the SNP at locus 250 that takes up two rows.

| locus | Wildtype | Wildtype_Freq | Mutant | Mutant_Freq |
|-------|----------|---------------|--------|-------------|
| 211   | A        | 2.8 [3]       | G      | 97.2 [104]  |
| 250   |          |               | G      | 3.8 [4]     |
| 250   |          |               | T      | 96.2 [102]  |
| 1116  |          |               | C      | 98.9 [88]   |
| 1116  |          |               | T      | 1.1 [1]     |


Below is the expected output (I moved the entries manually).

| locus | Wildtype | Wildtype_Freq | Mutant | Mutant_Freq |
|-------|----------|---------------|--------|-------------|
| 211   | A        | 2.8 [3]       | G      | 97.2 [104]  |
| 250   | G        | 3.8 [4]       | T      | 96.2 [102]  |
| 1116  | C        | 98.9 [88]     | T      | 1.1 [1]     |

NB: Also, the script does not incorporate information on wildtype vs. mutant alleles, 
hence, "Adenine" is assumed to always be the wild-type (see SNP at locus 211), followed 
by any other nucleotide in alphabetical order
