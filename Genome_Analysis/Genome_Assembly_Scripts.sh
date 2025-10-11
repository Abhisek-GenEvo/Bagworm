###################10X Genomics data assembly
supernova run --id=EC_assembly_output --sample=EC1,ec2,Ec3 --fastqs=EC_input1,EC_input2,EC_input3 --maxreads=all --localcores=40 --localmem=800 --nodebugmem --accept-extreme-coverage
supernova mkoutput --style=pseudohap2 --asmdir=EC_assembly_output/outs/assembly --outprefix=pseudo2_output_EC --minsize=0 --index
#Processing of the barcoded reads for downstream usage
longranger basic --id=longranger_out_EC --fastqs=EC_input1,EC_input2,EC_input3 --sample=EC1,ec2,Ec3 --localcores=40 --localmem=800
#Correction of the mis-assembly regions using Tigmint
samtools faidx pseudo2_output_EC.fasta
bwa index pseudo2_output_EC.fasta
bwa mem -t40 -p -C pseudo2_output_EC.fasta longranger_out_EC/outs/barcoded.fastq.gz | samtools sort -@8 -tBX -o draft.reads.sortbx.bam
tigmint-molecule draft.reads.sortbx.bam | sort -k1,1 -k2,2n -k3,3n > draft.reads.molecule.bed
tigmint-cut -p40 -o EC_assembly.fasta pseudo2_output_EC.fasta draft.reads.molecule.bed

###################Illumina assembly
#Data-filtration
java -jar /home/Softwares/Trimmomatic-0.39/trimmomatic-0.39.jar PE EC_Illumina_raw_R1.fastq EC_Illumina_raw_R2.fastq EC_illumina_R1.fastq EC_illumina_R1_unpaired.fastq EC_illumina_R2.fastq EC_illumina_R2_unpaired.fastq ILLUMINACLIP:/home/Softwares/Trimmomatic-0.39/adapters/all_trueseq_PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:60 -threads 15
#Assembly using SPAdes
spades.py --only-assembler -1 EC_illumina_R1.fastq -2 EC_illumina_R2.fastq -s EC_illumina_unpaired.fastq -k 107 -t 68 -m 970 -o SPADES_output_EC_k107

###################Scaffolding both 10X Genomics-based and Illumina-based assemblies (Here, the commands are given for only Illumina-based assembly)
#Using Illumina data
platanus_allee consensus -c EC_assembly.fasta -IP1 EC_illumina_R1.fastq EC_illumina_R2.fastq -t 40 -o EC_platanus 
#Using 10X Genomics data
arcs-make arcs draft=EC_platanus_consensusScaffold reads=barcoded_longranger 
#Using RNA-Seq data
/home/Softwares/augustus/bin/augustus --species=fly --gff3=on EC_scaffolded_10x.fa --outfile=augustus_EC_genome
hisat2-build EC_scaffolded_10x.fa Genome_Index
hisat2 -x Genome_Index -1 EC_RNAseq_R1_paired.fastq -2 EC_RNAseq_R2_paired.fastq -p 40 -S Illumina_Genome.sam
samtools view -S Illumina_Genome.sam -b -o Illumina_Genome.bam -@ 40
python2.7 /home/Softwares/AGOUTI/agouti.py scaffold -assembly EC_scaffolded_10x.fa -bam Illumina_Genome.bam -gff augustus_EC_genome.gff3 -outdir ./agouti_output
#The scaffolding steps were similarly performed for the Illumina data-based assembly

###################Merging of the assemblies
/home/Softwares/quickmerge-0.3-pl526he1b5a44_0/bin/merge_wrapper.py 10X_Genomics_assembly_AGOUTI_scaffolded.fasta Illumina_assembly_AGOUTI_scaffolded.fasta -l 38968 -ml 5000

###################Gap-closing
GapCloser -b config.txt -a Quickmerge_merged.fasta -l 155 -t 40 -o EC_assembly_gapclosed.fasta
##config.txt file format
-----
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
-----

###################Polishing 
bwa index EC_assembly_gapclosed.fasta
bwa mem -t 40 EC_assembly_gapclosed.fasta EC_illumina_R1.fastq EC_illumina_R2.fastq -o Pilon_EC.sam
samtools view -S Pilon_EC.sam -b -o Pilon_EC.bam -@ 40
samtools sort Pilon_EC.bam -o Pilon_EC_sorted.bam -@ 40 -m 10000000000
samtools index Pilon_EC_sorted.bam
pilon -Xmx600G --genome EC_assembly_gapclosed.fasta --bam Pilon_EC_sorted.bam --output Pilon_output_EC --outdir Pilon_outdir --threads 40 --changes

###################Filtering to keep the contigs >=5 Kb
perl calc_length_separate.pl Pilon_output_EC.fasta 4999
mv Pilon_output_EC.fasta.filtered Final_EC_assembly.fasta

###################Genome assembly quality evaluation
#Running QUAST
python /home/Softwares/quast-5.2.0/quast.py --min-contig 0 Final_EC_assembly.fasta -t 50
#Running BUSCO
/home/Softwares/busco-5.4.4/bin/busco --config config.ini -m genome -i Final_EC_assembly.fasta --augustus --augustus_species fly -o EC_BUSCO_output -l Insecta_odb10

###################Contamination removal from the assembly
#Mapping the reads onto the genome
bwa index Final_EC_assembly.fasta
bwa mem Final_EC_assembly.fasta EC_illumina_R1.fastq EC_illumina_R2.fastq -o mapping_output.sam -t 40 
/home/anaconda3/bin/samtools view -S mapping_output.sam -b -o mapping_output.bam -@ 40 
/home/anaconda3/bin/samtools sort mapping_output.bam -o mapping_output_sorted.bam -@ 40 -m 9000000000 
/home/anaconda3/bin/samtools index mapping_output_sorted.bam
#Mapping the genome against NCBI-nt database
/home/anaconda3/bin/blastn -db ./nt_database/nt -query Final_EC_assembly.fasta -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -max_target_seqs 10 -max_hsps 1 -evalue 1e-5 -num_threads 60 -out EC_blastn_for_blobtools
#Running BlobTools
conda activate blobtools2
blobtools taxify -f EC_blastn_for_blobtools -m nucl_gb.accession2taxid -s 1 -t 2 -o EC_blastn_blobtools_taxified.out
cut -f1,2,3 EC_blastn_blobtools_taxified.out > EC_blastn_blobtools_taxified_mod.out
blobtools create -i Final_EC_assembly.fasta -b mapping_output_sorted.bam -t EC_blastn_blobtools_taxified_mod.out -o BlobTools_EC --nodes data/nodes.dmp --names data/names.dmp
blobtools view -i BlobTools_EC.blobDB.json
blobtools plot -i BlobTools_EC.blobDB.json

###################Constructing a pseudochromosome-level genome assembly via genome-wide synteny alignment with Luffia ferchaultella genome assembly (GenBank - GCA_949709985.1)
export SATSUMA2_PATH=/home/vinit/SOFTWARES/satsuma2/bin/
Chromosemble -q BW_Final_CLEAN_genome.fasta -t LufFer_nuclear_genome.fasta -o BW_pseudochromosomes -n 60
