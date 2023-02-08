#UnclassifiedKrakenReadsEmmaJackson22/23

#Denovo genome assembly. This script was used to generate contigs through overlapping reads, obtained through shotgun sequencing. 

#!/bin/bash -e

#SBATCH -c 10 --mem=8Gb --account=<project_code> --time 00:30:00 -J 12566_SPAdes_EJ --array=0 --output=%A_%aemmaA.out --error=%A_%aemmaA.err

files=( $(ls *copyR1.fastq.gz) )
file=${files[$SLURM_ARRAY_TASK_ID]}

module purge
module load SPAdes/3.15.4-gimkl-2022a-Python-3.10.5

spades.py --isolate -1 $file -2 ${file/copyR1/copyR2} -o ./output/${file%%*}/

$module load BBMap/39.01-GCC-11.3.0
$stats.sh contigs.fasta

#Visualising reads mapped to the metagenome in IGV
#Sequences from the NCBI database were downloaded and uploaded to NeSI. 
#The genomes which were added to the metagenome was each top result from the initial BLAST search made from the contigs. This made up the main metagenome. 
#Those which were most commonly detected were added to another metagenome

#Main metagenome:
#Staphylococcus taiwanensis NTUH-S172
#Staphylococcus haemolyticus GDY8P80P
#Staphylococcus haemolyticus 4M 
#Staphylococcus haemolyticus BC5211
#Staphylococcus haemolyticus NY5
#Staphylococcus haemolyticus 12b
#Staphylococcus haemolyticus 7b
#Staphylococcus haemolyticus 1b
#Staphylococcus haemolyticus 20b
#Staphylococcus haemolyticus SCAID PHRX1-2019
#Staphylococcus haemolyticus SCAID URN1-2019
#Staphylococcus haemolyticus ATCC 29970
#Staphylococcus haemolyticus SGAir0252
#Staphylococcus haemolyticus SE2.14
#Staphylococcus haemolyticus SE3.8
#Staphylococcus haemolyticus SE3.9
#Staphylococcus haemolyticus VB5326
#Staphylococcus haemolyticus MSA_JNM56C1
#Staphylococcus aureus 280
#Staphylococcus aureus 281
#Staphylococcus aureus GHA4
#Staphylococcus aureus UP_1313
#Staphylococcus aureus NAS_AN_157
#Staphylococcus chromogenes 20B
#Staphylococcus chromogenes IVB6200
#Staphylococcus hominis S34-1
#Staphylococcus simiae NCTC13838
#Staphylococcus simulans NCTC7944
#Staphylococcus simulans IVB6207
#Staphylococcus equorum C2014
#Macrococcus epidermidis Epi3002-OL
#Streptococcus parauberis strain KSP10
#Streptococcus parauberis strain KSP14
#Macrococcus bohemicus strain 19Msa0383
#Macrococcus bohemicus strain 19Msa0936
#Macrococcus brunensis strain 18KM571
#Macrococcus brunensis strain 18KM1742
#Macrococcus equipercicus strain Epi0143-OL
#Macrococcus armenti JEK46
#Macrococcus armenti JEK12
#Macrococcus canis SD607
#Macrococcus canis Epi0100-OL
#Macrococcus sp IME1552
#Corynebacterium striatum strain KC-Na-01
#Corynebacterium resistens DSM 45100
#Corynebacterium kalinowskii strain 1959
#Arthrobacter koreensis strain DL
#Candida tropicalis strain MYA-3404
#Rothia nasimurium E1706032
#Rothia terrae KJZ-14
#Rothia amarae KJZ-9
#Aerococcus urinaeequi CCUG28094
#Pantoea agglomerans TH81
#Nocardia nova SH22a

#Most common metagenome:
#Staphylococcus chromogenes 20B
#Staphylococcus chromogenes IVB6200
#Staphylococcus haemolyticus 12b
#Staphylococcus haemolyticus 1b
#Staphylococcus haemolyticus 4M 
#Staphylococcus haemolyticus 7b
#Staphylococcus haemolyticus BC5211
#Staphylococcus haemolyticus GDY8P80P
#Staphylococcus haemolyticus SE2.14
#Staphylococcus haemolyticus SE3.8
#Staphylococcus haemolyticus SE3.9
#Staphylococcus taiwanensis NTUH-S172
#Macrococcus brunensis strain 18KM1742
#Macrococcus epidermidis Epi3002-OL

#Once the genomes were obtained, they were concatenated together using the cat command
cat <fasta file of all genomes> > ncbinewmetagenome.fa

cat <fasta file of most common genomes> > ncbimetagenome.fa

#The reads from each sequenced genome were mapped to the metagenome
#!/bin/bash -e
#SBATCH -c 10 --mem=8Gb --account=<project_code> --time 00:10:00 -J 12566_bwamem_EJ

module load BWA/0.7.17-GCC-9.2.0
module load SAMtools/1.13-GCC-9.2.0

REF=/nesi/project/nesi00187/emjac0/ncbimetagenome.fasta
bwa mem -t 6 -R"@RG\tID:<sample_code>\tPL:ILLUMINA\tSM:<sample_code>_March22" $REF /nesi/nobackup/nesi00187/Unclassified_Reads/sean/<sample_code>*R1.fastq.gz /nesi/nobackup/nesi00187/Unclassified_Reads/sean/<sample_code>*R2.fastq.gz | samtools view - -Sb | samtools sort - -@10 -o LIC353_Illuminareadsmappedtometagenome.cram

#The cram file is indexed using samtools
$module load SAMtools 
$samtools index <readmappingtosample.cram>

#The metagenome, cram file and cram.crai file and downloaded from NeSI
#The metagenome, cram file and cram.crai file are then transferred to IGV. Once IGV is opened, the tools button is pressed, then run igvtools is pressed. The input file is selected as the cram file, the output file is generated as a .cram.tdf file. This is inputted into IGV. 

#Making a phylogenetic tree 
#Use the metagenome as a reference to map all the sequenced files to it using the VariantCalling.sh script
#!/bin/bash -e
#SBATCH -c 6 --mem=16Gb --account=<project_code> --time 1-00 -J 12566_VariantCalling_EJ

module load BWA/0.7.17-GCC-7.4.0
module load SAMtools/1.9-GCC-7.4.0
module load BCFtools/1.9-GCC-7.4.0


mkdir -p METAbamemma_files_Sequenced_Files_April2022

DAT=`find ${1} -name "*R1*.fastq.gz"`
REF=/nesi/project/nesi00187/emjac0/metagenomes/<metagenome.fasta>
for i in ${DAT[@]}
do
        FQ1=${i}
        FQ2=`echo ${i} | sed s/R1/R2/`
        echo $FQ1
        echo $FQ2
        base=$(basename ${i} _R1.fastq.gz)
        echo "basename is ${base}"
        bwa mem -t 6 -R "@RG\tID:${base}\tPL:ILLUMINA\tSM:LIC_${base}" ${REF} ${FQ1} ${FQ2} | samtools sort -@ 6 -m 1G -O BAM -o METAbamemma_files_Sequenced_Files_April2022/LIC_${base}_Illuminaextractedreadsmappedtometa.bam

#Put all the paths to the bam files into a .txt file
$ find "$PWD" -name \*MSSA476.bam > locationofbamfilesGITHUB.txt

#Use mpileup to compare the sequence from each bam and the reference and combine these into one file
#!/bin/bash -e

#SBATCH -c 2 --mem=16Gb --account=<project_code> --time 72:00:00 -J 12566_bcftools_mpileup_EJ
module purge
module load GCCcore/11.3.0
module load BCFtools/1.15.1-GCC-11.3.0

bcftools mpileup -Ou -f /nesi/project/nesi00187/emjac0/metagenomes/ncbimetagenome.fasta -b /nesi/project/nesi00187/emjac0/Metagenome_Tree/locationofbamfiles.txt | bcftools call --ploidy 1 -mv -Ob -o /nesi/project/nesi00187/emjac0/meta_LIC_VariantCalls.bcf

$Bcftools view -Ov meta_LIC_VariantCalls.bcf > meta_LIC_VariantCalls.vcf

$module load VCF-kit/0.2.6-gimkl-2020a-Python-3.8.2
$vk phylo tree nj meta_LIC_VariantCalls.vcf > METAtreefeb.nwk


#Basecalling and further analysis of nanopore results
#Guppy - basecaller
#!/bin/bash -e

#SBATCH --account <project_code> -J 12566_guppy_gpu --gpus-per-node A100:1 --mem 6G --cpus-per-task 4 --time 10:00:00 --output slurmout.%j.out

module purge
module load ont-guppy-gpu/6.2.1

INPUT=/nesi/project/nesi00187/LIC_raw_reads/nanopore/EmmaNovel1/no_sample/20211210_1013_MC-110462_FAP17229_b3786ead
FOLDER=/nesi/project/nesi00187/LIC_raw_reads/nanopore/EmmaNovel1/Guppy_Super_Accurate_Model1

guppy_basecaller -i $INPUT -s $FOLDER -c /opt/nesi/CS400_centos7_bdw/ont-guppy-gpu/6.2.1/data/dna_r9.4.1_450bps_sup.cfg --device auto --recursive --detect_mid_strand_adapter --min_qscore 7 --barcode_kits SQK-NBD110-24

#Porechop - a tool for finding and removing adapters from nanopore reads
#!/bin/bash -e
#SBATCH -c 15 --mem=16Gb --account=<project_code> --time 00:120:00 -J 12566_porechop_EJ

module purge
module load Porechop/0.2.4-gimkl-2020a-Python-3.8.2

for input in /nesi/project/nesi00187/LIC_raw_reads/nanopore/EmmaNovel1/allbarcodesfastq/Emma_Barcode*_Guppy.fastq

do
base=$(basename ${input} _Guppy.fastq)
echo "working with file $input"
echo "basename is $base"

OUTPUT=/nesi/project/nesi00187/LIC_raw_reads/nanopore/EmmaNovel1/EmmaNovel_${base}_Guppy_Filtlong_Porechop.fastq
porechop -i $input -o $OUTPUT --discard_middle

#NanoStat - calculates various statistics
$Nanostat -fastq <filename.fastq> 

#NanoFilt - filters and trims long reads
$Nanofilt -q <mean read quality> | gzip > EmmaNovel_Emma_Barcode<barcode_number>_Guppy_Filtlong_Porechop.fastq.gz

#Canu - tool used to assemble long reads
#!/bin/bash -e
#SBATCH --nodes 1 --cpus-per-task 1 --ntasks 10 --mem=100Gb --account <project_code> --time 72:00:00 -J 12566_Canu_EJ --hint=nomultithread

module purge
module load Canu/2.1.1-GCC-9.2.0

#Need to change path of input to where the Nanostat/nanofilt/poerechop/guppy file is
#Change Directory to what name and path you want to store Canu in make this specific for each filename
#change name to be whatever you want to call the canu output files in the directory. i.e. saureus_barcodenumber

INPUT=/nesi/project/nesi00187/LIC_raw_reads/nanopore/EmmaNovel1/EmmaNovel_Emma_Barcode<barcode_number>_Guppy_Filtlong_Porechop.fastq.gz
DIRECTORY=/nesi/project/nesi00187/LIC_raw_reads/nanopore/EmmaNovel1/<barcode_file>
NAME=BAR<barcode_number>canuNovel


canu -d $DIRECTORY -p $NAME genomesize=2.8m useGrid=true gridOptions="--time=72:00:00 --account=<project_code> --hint=nomultithread" -nanopore-raw $INPUT

#Racon - polishing step which corrects any raw contigs

#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 10
#SBATCH --job-name 12566_ej_racon
#SBATCH --mem=20G
#SBATCH --time=00:40:00
#SBATCH --account=<project_code>
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

module load Racon/1.4.21-GCC-9.2.0-CUDA-11.2.0-hybrid
module load minimap2/2.24-GCC-9.2.0

READS=/nesi/project/nesi00187/LIC_raw_reads/nanopore/EmmaNovel1/EmmaNovel_Porechop_Nanofilt_Nanostat<barcode_number>.fastq.gz
CONTIG_FILE=/nesi/project/nesi00187/LIC_raw_reads/nanopore/EmmaNovel1/<barcode_file>/<barcode_number>canuNovel.contigs.fasta

#======> Correction 1 X Racon
minimap2 -t 10 $CONTIG_FILE $READS > NovelBAR<barcode_number>racon1.paf
racon -t 10 $READS NovelBAR<barcode_number>racon1.paf $CONTIG_FILE > NovelBAR<barcode_number>racon1.fasta

#======> Correction 2 X Racon
minimap2 -t 10 NovelBAR<barcode_number>racon1.fasta $READS > NovelBAR<barcode_number>racon2.paf
racon -t 10 $READS NovelBAR<barcode_number>racon2.paf NovelBAR<barcode_number>racon1.fasta > NovelBAR<barcode_number>racon2.fasta

#======> Correction 3 X Racon
minimap2 -t 10 NovelBAR<barcode_number>racon2.fasta $READS > NovelBAR<barcode_number>racon3.paf
racon -t 10 $READS NovelBAR<barcode_number>racon3.paf NovelBAR<barcode_number>racon2.fasta > NovelBAR<barcode_number>racon3.fasta

#======> Correction 4 X Racon
minimap2 -t 10 NovelBAR<barcode_number>racon3.fasta $READS > NovelBAR<barcode_number>racon4.paf
racon -t 10 $READS NovelBAR<barcode_number>racon4.paf NovelBAR<barcode_number>racon3.fasta > NovelBAR<barcode_number>racon4.fasta

#======> Correction 5 X Racon
#minimap2 -t 10 NovelBAR<barcode_number>racon4.fasta $READS > NovelBAR<barcode_number>racon5.paf
#racon -t 10 $READS NovelBAR<barcode_number>racon5.paf NovelBAR<barcode_number>racon4.fasta > NovelBAR<barcode_number>racon5.fasta

#BWA - used to map the sequences against a reference genome

$module load BWA/0.7.17-gimkl-2017a
$bwa index NovelBAR<barcode_number>racon5.fasta


#!/bin/bash -e
#SBATCH -c 10 --mem=8Gb --account=<project_code> --time 00:10:00 -J 12566_bwamem_EJ

module load  BWA/0.7.17-GCC-9.2.0
module load SAMtools/1.13-GCC-9.2.0

REF=/nesi/project/nesi00187/emjac0/NovelBAR<barcode_number>racon4.fasta
bwa mem -t 6 -R"@RG\tID:<sample_code>\tPL:ILLUMINA\tSM:<sample_code>_seqJan23" $REF /nesi/nobackup/nesi00187/Unclassified_Reads/LIC037_90%_Unclassified_Reads/<sample_code>*_R1.fastq.gz /nesi/nobackup/nesi00187/Unclassified_Reads/LIC037_90%_Unclassified_Reads/<sample_code>*_R2.fastq.gz | samtools view - -Sb | samtools sort - -@10 -o NovelBAR<barcode_number>_BWA.bam

#samtools - generate a BAM file
#!/bin/bash -e
#SBATCH -c 4 --mem=8Gb --account=<project_code> --time 00:10:00 -J 12566_samtoolssort_EJ

module load SAMtools

samtools sort --reference /nesi/project/nesi00187/emjac0/NovelBAR<barcode_number>racon4.fasta -o  NovelBAR<barcode_number>_<sample_code>_PolishedGenome.bam -O BAM /nesi/project/nesi00187/emjac0/NovelBAR<barcode_number>_BWA.bam

#Pilon - Polishes the assembled genome 
#!/bin/bash -e
#SBATCH -c 4 --mem=10Gb --account=<project_code> --time 00:10:00 -J 12566_PILON_EJ

module load Pilon/1.23-Java-1.8.0_144

mkdir -p /nesi/project/nesi00187/emjac0/pilon_genomeBAR<barcode_number>/

GENOMEPATH=/nesi/project/nesi00187/emjac0/pilon_genomeBAR<barcode_number>/

BAM=/nesi/project/nesi00187/emjac0/NovelBAR<barcode_number>_<sample_code>_PolishedGenome.bam

GENOME=/nesi/project/nesi00187/emjac0/NovelBAR<barcode_number>racon4.fasta

java -Xmx10G -jar $EBROOTPILON/pilon.jar --genome $GENOME --fix all --changes --frags $BAM --threads 10 --output $GENOMEPATH/pilon_round1 | tee $GENOMEPATH/round1.pilon

##First time you run this you are using what I have included in the brackets
##First time = genome is path to racon4 --frags is path to sortedbamfile from the samtools sort script --output is the path to pilon directory you make e.g. /nesi/project/nesi00187/spg12/pilon keep the pilonround1 same with the tee bit path to pilon directory /round1.pilon
##Second time you run this --genome pathtopilon1.fasta --frags pathtosortedpilon1.bam --output pilonpath/pilon_round2 tee pathtopilondirectory/round2.pilon
##Third time you run this --genome pathtopilon2.fasta --frags pathtosortedpilon2.bam --output pilonpath/pilon_round3 tee pathtopilondirectory/round3.pilon
##The final pilon polish the change file should be empty and the pilon3 is the hybrid genome then runKraken again

$bwa index pilon_round1.fasta 

#run bwamem_pilon.sh
#!/bin/bash -e
#SBATCH -c 5 --mem=10Gb --account=<project_code> --time 00:10:00 -J 12566_bwamempilon_EJ

module load GCC/9.2.0

module load BWA/0.7.17-GCC-9.2.0
module load SAMtools/1.13-GCC-9.2.0

#First step is your path to pilon_round1.fasta
##Second step is the paths to illumina reads 1 and 2
##Third step is the output path to sorted bam file i gave a example bam file

##Second time you run this first step pathtopilon_round2.fasta, the illumina read paths stay the same, the output changes to Pilon2_sorted.bam


bwa mem -t 5 /nesi/project/nesi00187/emjac0/pilon_genomeBAR<barcode_number>/pilon_round2.fasta /nesi/nobackup/nesi00187/Unclassified_Reads/LIC037_90%_Unclassified_Reads/<sample_code>*R1.fastq.gz /nesi/nobackup/nesi00187/Unclassified_Reads/LIC037_90%_Unclassified_Reads/<sample_code>*R1.fastq.gz | samtools view - -Sb | samtools sort - -@14 -o /nesi/project/nesi00187/emjac0/pilon2_BAR<barcode_number>_sorted.bam

$samtools index pilon1_sorted.bam

#run pilon.sh again

$bwa index pilon_round2.fasta

#run bwamem_pilon.sh again

$samtools index pilon2_sorted.bam

#run pilon.sh again

#Kraken2 
$sbatch -J 12566_Kraken2Nanopore --mem 75G --time 00:05:00 -c 12 --account <project_code> --wrap 'set -x; module load Kraken2; kraken2 --confidence 0.1 --db /nesi/project/nesi00187/anwal2/KrakenDBs/Dec21/ --report BAR<barcode_number>nanopore_NovelEmmaseqJan23.report --thread 24 /nesi/project/nesi00187/LIC_raw_reads/nanopore/EmmaNovel1/BAR<barcode_number>/BAR<barcode_number>canuNovel.contigs.fasta '
