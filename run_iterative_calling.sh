#### This scripts identifies the most closely related reference sequence based on the results produced by run_blast.sh and uses it to iteratively call consensus sequence

# software requirements:
# seqtk
# samtools-1.2
# smalt
# ivar
# bedtools-2.29.2
# lofreq_star-2.1.5
# htslib
# bcftools-1.10.2
# bedops

# ${sample} is the name of the sample passed via qsub -v sample=sample_name..

# identify reference LANL sequence as the sequence supported by the majority of sequencing reads among sequences collected in Russia (RU)
awk '{print $2}' ${sample}.sorted.out | sort | uniq -c | sort -k1,1gr | grep -m 1 RU | awk '{print $2}' > ${sample}.ref
refname=`cat ${sample}.ref`
cat ${sample}.ref | sed 's/^/>/' > ${sample}.ref.fasta
seqtk subseq ref_lanl.fasta ${sample}.ref | grep -v ">" | tr -d '\n' | sed 's/-//g' >> ${sample}.ref.fasta
ref=${sample}.ref.fasta
/home/tools/samtools/samtools-1.2/samtools faidx $ref


# trim primers from reads; for this purpose, i'll map reads onto the reference genome (HXB2; K03455). Smalt mapping parameters should be softer!
reads1=${sample}.R1.fq
reads2=${sample}.R2.fq

# bed/ is a directory that contains primer coordinates in BED format for each reference sequence; example is provided in the repository
bed=bed/$refname.bed

smalt index -k 10 -s 3 ref $ref
smalt map -o raw.bam -n 1 -i 1000 -x ref $reads1 $reads2
/home/tools/samtools/samtools-1.2/samtools view -f 3 -b -o f3.raw.bam raw.bam
/home/tools/samtools/samtools-1.2/samtools sort -m 10000000000 f3.raw.bam f3.raw.sorted
/home/tools/samtools/samtools-1.2/samtools index f3.raw.sorted.bam

# ivar softclips primers in the provided .bam file based on the primer coordinates from $bed
/home/noraneko/Software/bin/ivar trim -i f3.raw.sorted.bam -b $bed -e -p ivar.bam
/home/tools/samtools/samtools-1.2/samtools sort -m 10000000000 ivar.bam ivar.sorted
/home/tools/samtools/samtools-1.2/samtools index ivar.sorted.bam

# python script removes softclipped regions; https://github.com/ngsutils/ngsutils/blob/master/ngsutils/bam/removeclipping.py
python removeclipping.py -f ivar.sorted.bam ivar.sorted.clipped.bam

/home/tools/samtools/samtools-1.2/samtools index ivar.sorted.clipped.bam
# samtools1.2 sort -n doesn't work for some reason
/home/tools/samtools/samtools-1.10/bin/samtools sort -n -o ivar.namesorted.clipped.bam ivar.sorted.clipped.bam
/home/tools/bedtools/bedtools-2.29.2/bin/bedtools bamtofastq -i ivar.namesorted.clipped.bam -fq r1.clipped.fastq -fq2 r2.clipped.fastq

# reads are now primer-free!
reads1=r1.clipped.fastq
reads2=r2.clipped.fastq

# the reference sequence is a full-length genome; let's crop it to the region covered by reads
smalt index -k 10 -s 3 ref $ref
smalt map -o temp.bam -n 1 -i 1000 -x ref $reads1 $reads2
/home/tools/samtools/samtools-1.2/samtools sort -m 10000000000 temp.bam temp.sorted
/home/tools/samtools/samtools-1.2/samtools index temp.sorted.bam
bam=temp.sorted.bam

/home/tools/bedtools/bedtools-2.29.2/bin/bamToBed -i $bam > ${sample}.bed
/home/tools/bedtools/bedtools-2.29.2/bin/mergeBed -i ${sample}.bed -d 100 | awk '{if ($3-$2>1000) print $0}' > ${sample}.reg
cat ${sample}.reg | awk '{$2=$2-1; print $1"\t"$2"\t"$3}' > ${sample}.0based.reg
nregions=`cat ${sample}.0based.reg | wc -l`

# filter the required pol gene region, in case the regions outside the required fragment got non-zero coverage
if (( $nregions > 1 ))
then
   awk '{if ($3<3500) print $0}' ${sample}.0based.reg > temp
   mv temp ${sample}.0based.reg
 fi

/home/tools/bedtools/bedtools-2.29.2/bin/fastaFromBed -fi ${sample}.ref.fasta -bed ${sample}.0based.reg | sed 's/:.*//' > cropped.${sample}.fasta

# $ref is the starting reference sequence used in the iterative calling procedure below
ref=cropped.${sample}.fasta
/home/tools/samtools/samtools-1.2/samtools faidx $ref


# usually converges after 2 or 3 steps
for iter in 1 2 3 4 5
do

  smalt index -k 10 -s 3 ref $ref
  smalt map -o iter.bam -n 1 -i 1000 -x ref $reads1 $reads2
  /home/tools/samtools/samtools-1.2/samtools sort -m 10000000000 iter.bam iter.sorted
  /home/tools/samtools/samtools-1.2/samtools index iter.sorted.bam
  bam=iter.sorted.bam

  /home/noraneko/Software/lofreq_star-2.1.5/bin/lofreq viterbi --ref $ref $bam > viterbi.bam
  /home/tools/samtools/samtools-1.2/samtools sort -m 10000000000 viterbi.bam viterbi.sorted
  /home/tools/samtools/samtools-1.2/samtools index viterbi.sorted.bam
  viterbi=viterbi.sorted.bam
  rm iq.$viterbi
  /home/noraneko/Software/lofreq_star-2.1.5/bin/lofreq indelqual -f $ref -o iq.$viterbi --dindel $viterbi
  iq=iq.$viterbi
  /home/tools/samtools/samtools-1.2/samtools index $iq
  
  #call variants using  min depth 2
  vcf=$iter.lofreq.mindepth4.nofilt.vcf
  /home/noraneko/Software/lofreq_star-2.1.5/bin/lofreq call -f $ref -o $vcf -C 4 --call-indels --no-default-filter --force-overwrite --use-orphan $iq
  vars=$iter.lofreq.dp4.vars.vcf

  rm $vars
  /home/noraneko/Software/lofreq_star-2.1.5/bin/lofreq filter -Q 20 -K 20 --no-defaults -v 4 -V 0 -a 0.500001 -A 0 -i $vcf -o $vars
  new_vars=`grep -v "#" $vars`

  if [ "$new_vars" == "" ]
  then
    echo $sample, $iter
    break
  fi

  /home/tools/htslib/htslib-1.11/bin/bgzip -c $vars > $vars.gz
  /home/tools/htslib/htslib-1.11/bin/tabix -p vcf $vars.gz
  vars=$vars.gz

  /home/tools/bedtools/bedtools-2.29.2/bin/bedtools genomecov -bga -ibam $iq | awk '$4<4' | /home/tools/bedtools/bedtools-2.29.2/bin/bedtools merge > $iter.low_cov.bed
  zcat $vars | grep INDEL | grep -v "#" | /home/noraneko/Software/bedops/bin/vcf2bed --deletions | /home/noraneko/Software/bedops/bin/bedops --range 1:0 --everything - > $iter.deletions.bed
  /home/noraneko/Software/bedops/bin/bedops --difference $iter.low_cov.bed $iter.deletions.bed > $iter.mask.bed

  /home/noraneko/Software/bcftools-1.10.2/bcftools index -f $vars
  /home/noraneko/Software/bcftools-1.10.2/bcftools consensus -f $ref -m $iter.mask.bed $vars > $iter.consensus.fa

  ref=$iter.consensus.fa
  /home/tools/samtools/samtools-1.2/samtools faidx $ref

done

cat $ref | sed "s/>.*/>$sample/" > consensus.fa

### we may also want to call ambiguous nucleotides
mkdir call_ambiguous
cd call_ambiguous
cp ../consensus.fa .
ref=consensus.fa
/home/tools/samtools/samtools-1.2/samtools faidx $ref
# take the vcf file from the latest iteration
vcf=`ls -tr ../*.lofreq.mindepth4.nofilt.vcf | tail -n 1`
vars=lofreq.dp4.vars.vcf
rm $vars
# filter minor variants at frequency of at least 20%
/home/noraneko/Software/lofreq_star-2.1.5/bin/lofreq filter -Q 20 -K 20 --only-snvs --no-defaults -v 4 -V 0 -a 0.2 -A 0 -i $vcf -o $vars

# a custom python script that applies minor variants to produce a degenerate consensus
python ../split_calls.py -i $vars -f $ref -o deg_ref.fa
