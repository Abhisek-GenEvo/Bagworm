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

#Scaffolding of both assemblies
platanus_allee consensus -c EC_assembly.fasta -IP1 BW_R1_paired.fastq BW_R2_paired.fastq -t 46 -o BW_platanus
arcs-make arcs draft=BW_platanus_consensusScaffold_single reads=barcoded
LINKS -f BW_arcs_scaffolded.fasta -s empty.fof -b BW_scaffolding_ont -l 5 -t 2 -a 0.3

#Merging of the assemblies
./quickmerge-0.3-pl526he1b5a44_0/bin/merge_wrapper.py agouti_10x.fasta agouti_Illumina_annot.fasta -l 38968 -ml 5000

#Gap-closing
GapCloser -b config.txt -a merged_out.fasta -l 155 -t 46 -o BW_merged_Illu_gapclosed.fasta

#Polishing 

#Filtering to keep the contigs >=5 Kb
perl calc_length_separate.pl Final_pilon3_output_bw_single.fasta 4999
