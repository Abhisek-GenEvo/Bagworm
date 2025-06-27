#10X Genomics data assembly
supernova run --id=Assembly_Output --sample=BW,bw,Bw --fastqs=input1,input2,input3 --maxreads=all --localcores=40 --localmem=980 --nodebugmem --accept-extreme-coverage
supernova mkoutput --style=pseudohap2 --asmdir=Assembly_Output/outs/assembly --outprefix=pseudo2_output --minsize=0 --index

longranger basic --id=longranger_out_BW --fastqs=/DATA/BW_abhisek/BW_input1,/DATA/BW_abhisek/BW_input2,/DATA/BW_abhisek/BW_input3 --sample=BW,bw,Bw --localcores=36 --localmem=500
samtools faidx pseudo2_BW_All_output_single1.fasta
bwa index pseudo2_BW_All_output_single1.fasta
bwa mem -t46 -p -C pseudo2_BW_All_output_single1.fasta /DATA/BW_abhisek/longranger_out_BW/outs/barcoded.fastq.gz | samtools sort -@8 -tBX -o draft.reads.sortbx.bam
tigmint-molecule draft.reads.sortbx.bam | sort -k1,1 -k2,2n -k3,3n > draft.reads.molecule.bed
tigmint-cut -p46 -o draft.tigmint.fa pseudo2_BW_All_output_single1.fasta draft.reads.molecule.bed

#Illumina assembly
spades.py --only-assembler -1 BW_illumina_10x_R1.fastq -2 BW_illumina_10x_R2.fastq -s BW_illu_10x_unpaired.fastq -k 107 -t 68 -m 970 -o SPADES_Output_final_BW_107

#Scaffolding of both assemblies
platanus_allee consensus -c assembly.fasta -IP1 BW_R1_paired.fastq BW_R2_paired.fastq -t 46 -o BW_platanus
arcs-make arcs draft=BW_platanus_consensusScaffold_single reads=barcoded
LINKS -f BW_arcs_scaffolded.fasta -s empty.fof -b BW_scaffolding_ont -l 5 -t 2 -a 0.3

#Merging of the assemblies
./quickmerge-0.3-pl526he1b5a44_0/bin/merge_wrapper.py agouti_10x.fasta agouti_Illumina_annot.fasta -l 38968 -ml 5000

#Gap-closing
GapCloser -b config.txt -a merged_out.fasta -l 155 -t 46 -o BW_merged_Illu_gapclosed.fasta

#Polishing 

#Filtering to keep the contigs >=5 Kb
perl calc_length_separate.pl Final_pilon3_output_bw_single.fasta 4999
