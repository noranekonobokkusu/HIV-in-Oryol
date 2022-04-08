## HIV-in-Oryol
# Molecular epidemiology of HIV-1 in Oryol Oblast, Russia

This repository contains iterative consensus calling pipeline used in 10.1101/2021.10.26.21265513.
The pipeline consists of the two parts:
1. identifying the closest reference sequence from a set of given sequences based on the <b>blastn</b> results
2. iterative mapping and variant calling of sequencing reads on the selected reference sequence
