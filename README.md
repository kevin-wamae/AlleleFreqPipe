# R-Toolkit - FreqAllele
An R script for generating allele frequencies from DNA/amino acid sequence alignments.

**NB**: _It's worth noting that the script (inside the `01_R_scripts` directory) ignores indels/gaps in the alignment_

* Input is an alignment in Phylip format which can be generated in many programs such as [Aliview](https://ormbunkar.se/aliview/).

* Incase you choose to use Aliview, load your alignmnent and click on:
  * `File -> Save as Phylip (full names & padded)`

* Point the script to your input file in the `03_data_input` directory (line 25).

* Run the entire script and you'll find your output in the `04_data_output` directory, which will look like the example below.

* The dominant/major allele is grouped under `major`, while minor variants are grouped under `minor_1`, `minor_2`, e.t.c.

| locus | major         | minor_1       | minor_2      |
|-------|---------------|---------------|--------------|
| 64    | A / 98 [49]   | G / 2 [1]     |              |
| 79    | A / 97.9 [47] | G / 2.1 [1]   |              |
| 241   | A / 98 [50]   | T / 2 [1]     |              |
| 448   | T / 96 [48]   | A / 4 [2]     |              |
| 467   | G / 98 [49]   | A / 2 [1]     |              |
| 470   | A / 82 [41]   | G / 18 [9]    |              |
| 502   | G / 98.1 [53] | C / 1.9 [1]   |              |
| 515   | A / 98.2 [56] | G / 1.8 [1]   |              |
| 592   | A / 98.6 [69] | T / 1.4 [1]   |              |
| 593   | C / 98.6 [69] | A / 1.4 [1]   |              |
| 676   | A / 93.4 [71] | G / 3.9 [3]   | R / 2.6 [2]  |
| 820   | G / 56.2 [45] | A / 40 [32]   | R / 3.8 [3]  |
| 835   | G / 65 [52]   | A / 33.8 [27] | R / 1.2 [1]  |
