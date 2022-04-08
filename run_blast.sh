### this script blasts trimmed paired Illumina reads against a set of reference HIV genomes lanl_refs.fasta (can be downloaded from https://www.hiv.lanl.gov/content/sequence/NEWALIGN/align.html)

# create an indexed database for blastn
/home/tools/blast-2.2.31/bin/makeblastdb -in lanl_refs.fasta -out lanl -dbtype 'nucl'

# ${sample} is the name of the sample passed via qsub -v sample=sample_name..
# convert .fastq read files to fasta
seqtk fq2fa ${sample}.R1.fq > ${sample}.fasta
seqtk fq2fa ${sample}.R2.fq >> ${sample}.fasta

# run blastn
/home/tools/blast/blast-2.2.31/bin/blastn -db lanl -query ${sample}.fasta -outfmt "6 qseqid sseqid pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore" -out ${sample}.out -evalue 1e-10 -num_threads 1

# sort blastn results to include only one best hit per read, based on bitscore and evalue 
sort -k1,1 -k13,13gr -k12,12g ${sample}.out | sort -u -k1,1 --merge > ${sample}.sorted.out

