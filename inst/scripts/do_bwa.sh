#!/bin/bash
#./do_bwa.sh
#SBATCH -n6 -t 4-0 -p campus --mail-type=ALL --mail-user=cwon2@fhcrc.org -A tapscott_s
fq1=$1
bwaGenome=$2
bwaOutDir=$3
base=`basename ${fq1}`
sampleName=`echo ${base//.fastq}`
echo $sampleName

# use BWA
bwaVersion=/home/solexa/apps/bwa/bwa-0.7.10
## bwaGenome=/shared/solexa/solexa/Genomes/iGenome/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa

## need to define mainDir and sampleDir?
mkdir -p $bwaOutDir

# umask 0002
export PATH=/usr/kerberos/bin:/usr/local/bin:/bin:/usr/bin:/usr/X11R6/bin:/opt/moab/bin:$bwaVersion:/home/solexa/apps/samtools/samtools-0.1.19:/home/solexa/apps/FastQC

# Alignment
bwa aln -t 6 -R 1 -n 1 $bwaGenome $fq1 > $bwaOutDir/$sampleName.sai

bwa samse -n 1 $bwaGenome $bwaOutDir/$sampleName.sai $fq1 | samtools view -bS -> $bwaOutDir/$sampleName.bam

# clean up
rm $bwaOutDir/$sampleName.sai

# sort and index bam
cd $bwaOutDir
samtools sort $sampleName.bam $sampleName.bam.sorted
mv $sampleName.bam.sorted.bam $sampleName.bam
samtools index $sampleName.bam
touch $sampleName.bwaDone.txt
exit 0


