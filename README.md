# HIV consensus calling
## Molecular epidemiology of HIV-1 in Oryol Oblast, Russia

This repository contains iterative consensus calling pipeline used in 10.1101/2021.10.26.21265513.
The pipeline consists of the two parts:
1. <b>blastn</b> sequencing reads of a sample against a set of curated reference HIV genomes to identify the closest reference sequence
2. iterative mapping and variant calling of sequencing reads on the selected reference sequence
