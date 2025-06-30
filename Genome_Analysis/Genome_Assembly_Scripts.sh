#10X Genomics data assembly
supernova run --id=EC_assembly_output --sample=EC1,ec2,Ec3 --fastqs=EC_input1,EC_input2,EC_input3 --maxreads=all --localcores=40 --localmem=980 --nodebugmem --accept-extreme-coverage
supernova mkoutput --style=pseudohap2 --asmdir=EC_assembly_output/outs/assembly --outprefix=pseudo2_output_EC --minsize=0 --index

longranger basic --id=longranger_out_EC --fastqs=/DATA/EC/EC_input1,/DATA/EC/EC_input2,/DATA/EC/EC_input3 --sample=EC1,ec2,Ec3 --localcores=36 --localmem=500
samtools faidx pseudo2_output_EC.fasta
bwa index pseudo2_output_EC.fasta
bwa mem -t46 -p -C pseudo2_output_EC.fasta /DATA/EC/longranger_out_EC/outs/barcoded.fastq.gz | samtools sort -@8 -tBX -o draft.reads.sortbx.bam
tigmint-molecule draft.reads.sortbx.bam | sort -k1,1 -k2,2n -k3,3n > draft.reads.molecule.bed
tigmint-cut -p46 -o draft.tigmint.fa pseudo2_output_EC.fasta draft.reads.molecule.bed

#Illumina assembly
spades.py --only-assembler -1 EC_illumina_R1.fastq -2 EC_illumina_R2.fastq -s EC_illumina_unpaired.fastq -k 107 -t 68 -m 970 -o SPADES_output_EC_k107

#Scaffolding both 10X Genomics and Illumina-based assemblies
#Using Illumina data
platanus_allee consensus -c EC_assembly.fasta -IP1 EC_illumina_R1.fastq EC_illumina_R2.fastq -t 46 -o EC_platanus 
#Using 10X Genomics data
arcs-make arcs draft=EC_platanus_consensusScaffold reads=barcoded_longranger 
LINKS -f EC_arcs_scaffolded.fasta -s empty.fof -b EC_scaffolding_ont -l 5 -t 2 -a 0.3 
#Using RNA-Seq data

#Merging of the assemblies
/home/Softwares/quickmerge-0.3-pl526he1b5a44_0/bin/merge_wrapper.py 10X_Genomics_scaffolded.fasta Illumina_scaffolded.fasta -l 38968 -ml 5000

#Gap-closing
GapCloser -b config.txt -a Quickmerge_merged.fasta -l 155 -t 46 -o EC_assembly_gapclosed.fasta
config.txt file format:
#maximal read length
max_rd_len=155
[LIB]
#average insert size
avg_ins=350
#if sequence needs to be reversed
reverse_seq=0
#in which part(s) the reads are used
asm_flags=4
#use only first 100 bps of each read
rd_len_cutoff=150
q1=EC_illumina_R1.fastq
q2=EC_illumina_R2.fastq
q=EC_illumina_unpaired.fastq

#Polishing 
bwa index EC_assembly_gapclosed.fasta
bwa mem -t 70 -M EC_assembly_gapclosed.fasta EC_illumina_R1.fastq EC_illumina_R2.fastq -o Pilon_EC.sam
samtools view -S Pilon_EC.sam -b -o Pilon_EC.bam -@ 70
samtools sort Pilon_EC.bam -o Pilon_EC_sorted.bam -@ 70 -m 10000000000
samtools index Pilon_EC_sorted.bam
pilon -Xmx950G --genome EC_assembly_gapclosed.fasta --bam Pilon_EC_sorted.bam --output Pilon_output_EC --outdir Pilon_outdir --threads 70 --changes

#Filtering to keep the contigs >=5 Kb
perl calc_length_separate.pl Pilon_output_EC.fasta 4999
mv Pilon_output_EC.fasta.annot Final_EC_assembly.fasta
