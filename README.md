# unclassifiedreadsEJ2023-23

#Denovo genome assembly. This script was used to generate contigs through overlapping reads, obtained through shotgun sequencing. 

#!/bin/bash -e

#SBATCH -c 10 --mem=8Gb --account=<project_code> --time 00:30:00 -J 12566_SPAdes_EJ --array=0 --output=%A_%aemmaA.out --error=%A_%aemmaA.err

files=( $(ls *copyR1.fastq.gz) )
file=${files[$SLURM_ARRAY_TASK_ID]}

module purge
module load SPAdes/3.15.4-gimkl-2022a-Python-3.10.5

spades.py --isolate -1 $file -2 ${file/copyR1/copyR2} -o ./output/${file%%*}/

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
#Macrococcus armenti JEK29
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


