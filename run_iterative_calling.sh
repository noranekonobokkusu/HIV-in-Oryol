# software requirements:
# seqtk
# samtools-1.2

# awk '{print $2}' /home/noraneko/projects/hiv/11_10_20_call_consensus/blast_reads/$sample/${sample}.sorted.out | sort | uniq -c | sort -k1,1gr | grep -m 1 RU | awk '{print $2}' > ${sample}.ref
# refname=`cat ${sample}.ref`
# cat ${sample}.ref | sed 's/^/>/' > ${sample}.ref.fasta
# seqtk subseq /home/noraneko/projects/hiv/11_10_20_call_consensus/blast_reads/HIV1_COM_2017_genome_DNA.fasta ${sample}.ref | grep -v ">" | tr -d '\n' | sed 's/-//g' >> ${sample}.ref.fasta
# ref=${sample}.ref.fasta
# /home/tools/samtools/samtools-1.2/samtools faidx $ref
#
# ### trim primers for reads; for this purpose, i'll map reads onto the reference genome (hxb2). Smalt mapping parameters should be softer!
#
# reads1=/home/noraneko/projects/hiv/11_10_20_reads_grouped_by_samples/${sample}.R1.fq
# reads2=/home/noraneko/projects/hiv/11_10_20_reads_grouped_by_samples/${sample}.R2.fq
# bed=/home/noraneko/projects/hiv/ref_alignments/bed/$refname.bed
#
# smalt index -k 10 -s 3 ref $ref
# smalt map -o raw.bam -n 1 -i 1000 -x ref $reads1 $reads2
# /home/tools/samtools/samtools-1.2/samtools view -f 3 -b -o f3.raw.bam raw.bam
# /home/tools/samtools/samtools-1.2/samtools sort -m 10000000000 f3.raw.bam f3.raw.sorted
# /home/tools/samtools/samtools-1.2/samtools index f3.raw.sorted.bam
# #
# /home/noraneko/Software/bin/ivar trim -i f3.raw.sorted.bam -b $bed -e -p ivar.bam
# /home/tools/samtools/samtools-1.2/samtools sort -m 10000000000 ivar.bam ivar.sorted
# /home/tools/samtools/samtools-1.2/samtools index ivar.sorted.bam
#
# python ../removeclipping.py -f ivar.sorted.bam ivar.sorted.clipped.bam
#
# /home/tools/samtools/samtools-1.2/samtools index ivar.sorted.clipped.bam
# ### samtools1.2 sort -n doesn't work..
# /home/tools/samtools/samtools-1.10/bin/samtools sort -n -o ivar.namesorted.clipped.bam ivar.sorted.clipped.bam
# /home/tools/bedtools/bedtools-2.29.2/bin/bedtools bamtofastq -i ivar.namesorted.clipped.bam -fq r1.clipped.fastq -fq2 r2.clipped.fastq
#
#
# reads1=r1.clipped.fastq
# reads2=r2.clipped.fastq
#
# smalt index -k 10 -s 3 ref $ref
# smalt map -o temp.bam -n 1 -i 1000 -x ref $reads1 $reads2
# /home/tools/samtools/samtools-1.2/samtools sort -m 10000000000 temp.bam temp.sorted
# /home/tools/samtools/samtools-1.2/samtools index temp.sorted.bam
# bam=temp.sorted.bam
#
# /home/tools/bedtools/bedtools-2.29.2/bin/bamToBed -i $bam > ${sample}.bed
# /home/tools/bedtools/bedtools-2.29.2/bin/mergeBed -i ${sample}.bed -d 100 | awk '{if ($3-$2>1000) print $0}' > ${sample}.reg
# cat ${sample}.reg | awk '{$2=$2-1; print $1"\t"$2"\t"$3}' > ${sample}.0based.reg
# nregions=`cat ${sample}.0based.reg | wc -l`
# if (( $nregions > 1 ))
# then
#    awk '{if ($3<3500) print $0}' ${sample}.0based.reg > temp
#    mv temp ${sample}.0based.reg
#  fi
#
# /home/tools/bedtools/bedtools-2.29.2/bin/fastaFromBed -fi ${sample}.ref.fasta -bed ${sample}.0based.reg | sed 's/:.*//' > cropped.${sample}.fasta
#
# ref=cropped.${sample}.fasta
# /home/tools/samtools/samtools-1.2/samtools faidx $ref
#
# for iter in 1 2 3 4 5
# do
#
#   smalt index -k 10 -s 3 ref $ref
#   smalt map -o iter.bam -n 1 -i 1000 -x ref $reads1 $reads2
#   /home/tools/samtools/samtools-1.2/samtools sort -m 10000000000 iter.bam iter.sorted
#   /home/tools/samtools/samtools-1.2/samtools index iter.sorted.bam
#   bam=iter.sorted.bam
#
#   /home/noraneko/Software/lofreq_star-2.1.5/bin/lofreq viterbi --ref $ref $bam > viterbi.bam
#   /home/tools/samtools/samtools-1.2/samtools sort -m 10000000000 viterbi.bam viterbi.sorted
#   /home/tools/samtools/samtools-1.2/samtools index viterbi.sorted.bam
#   viterbi=viterbi.sorted.bam
#   rm iq.$viterbi
#   /home/noraneko/Software/lofreq_star-2.1.5/bin/lofreq indelqual -f $ref -o iq.$viterbi --dindel $viterbi
#   iq=iq.$viterbi
#   /home/tools/samtools/samtools-1.2/samtools index $iq
#   #
#   #  #call variants using  min depth 2
#   vcf=$iter.lofreq.mindepth4.nofilt.vcf
#   /home/noraneko/Software/lofreq_star-2.1.5/bin/lofreq call -f $ref -o $vcf -C 4 --call-indels --no-default-filter --force-overwrite --use-orphan $iq
#   vars=$iter.lofreq.dp4.vars.vcf
#
#
#   rm $vars
#   /home/noraneko/Software/lofreq_star-2.1.5/bin/lofreq filter -Q 20 -K 20 --no-defaults -v 4 -V 0 -a 0.500001 -A 0 -i $vcf -o $vars
#   new_vars=`grep -v "#" $vars`
#
#   if [ "$new_vars" == "" ]
#   then
#     echo $sample, $iter
#     break
#   fi
#
#   /home/tools/htslib/htslib-1.11/bin/bgzip -c $vars > $vars.gz
#   /home/tools/htslib/htslib-1.11/bin/tabix -p vcf $vars.gz
#   vars=$vars.gz
#
#   /home/tools/bedtools/bedtools-2.29.2/bin/bedtools genomecov -bga -ibam $iq | awk '$4<4' | /home/tools/bedtools/bedtools-2.29.2/bin/bedtools merge > $iter.low_cov.bed
#   zcat $vars | grep INDEL | grep -v "#" | /home/noraneko/Software/bedops/bin/vcf2bed --deletions | /home/noraneko/Software/bedops/bin/bedops --range 1:0 --everything - > $iter.deletions.bed
#   /home/noraneko/Software/bedops/bin/bedops --difference $iter.low_cov.bed $iter.deletions.bed > $iter.mask.bed
#
#   /home/noraneko/Software/bcftools-1.10.2/bcftools index -f $vars
#   /home/noraneko/Software/bcftools-1.10.2/bcftools consensus -f $ref -m $iter.mask.bed $vars > $iter.consensus.fa
#
#   ref=$iter.consensus.fa
#   /home/tools/samtools/samtools-1.2/samtools faidx $ref
#
# done
#
# cat $ref | sed "s/>.*/>$sample/" > consensus.fa

cd ${sample}
mkdir call_ambiguous
cd call_ambiguous
cp ../consensus.fa .
ref=consensus.fa
/home/tools/samtools/samtools-1.2/samtools faidx $ref
vcf=`ls -tr ../*.lofreq.mindepth4.nofilt.vcf | tail -n 1`
vars=lofreq.dp4.vars.vcf
rm $vars
/home/noraneko/Software/lofreq_star-2.1.5/bin/lofreq filter -Q 20 -K 20 --only-snvs --no-defaults -v 4 -V 0 -a 0.2 -A 0 -i $vcf -o $vars
cat $vars | grep -v "^#" | cut -f 2 | sort | uniq -c | grep -v "1 " | awk '{print $2}' > multiple_alleles.txt
python ../../split_calls.py -i $vars -p multiple_alleles.txt -f $ref -o new_ref.fa
ref=new_ref.fa
/home/tools/samtools/samtools-1.2/samtools faidx $ref
smalt index -k 10 -s 3 ref $ref
smalt map -o iter.bam -n 1 -i 1000 -x ref $reads1 $reads2
/home/tools/samtools/samtools-1.2/samtools sort -m 10000000000 iter.bam iter.sorted
/home/tools/samtools/samtools-1.2/samtools index iter.sorted.bam
bam=iter.sorted.bam
