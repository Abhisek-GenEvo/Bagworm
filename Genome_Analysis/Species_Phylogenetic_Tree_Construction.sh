###################Orthogroups construction
orthofinder -f ./input_files_dir/ -S diamond -t 40 -a 40

###################Extraction of the fuzzy one-to-one orthogroups
python2.7 /home/Softwares/kinfin/src/kinfin.py -p SpeciesIDs.txt -g Orthogroups.txt -c config_working.txt -o kinfin_EC_phylogeny -s SequenceIDs.txt -a ./fasta_files
#Fuzzy one-to-one orthogroup IDs were extracted from the output file - all.all.cluster_1to1s.txt 
python2.7 /home/Softwares/kinfin/scripts/get_protein_ids_from_cluster.py -g Orthogroups.txt --cluster_ids 1to1_55speciesIDs.txt
python pick_longest_ortholog.py 
for i in ./138_Longest_Orthogroups/*; do cut -d"#" -f1,2 "$i" > "$i".fasta; done;

###################Orthogroups alignment
parallel -j 20 '/home/anaconda3/bin/mafft {} > {/.}.einsi.faa' ::: ./138_Longest_Orthogroups_mod/*.fasta

###################Alignment processing
perl /home/Softwares/BeforePhylo/BeforePhylo.pl -type=protein -conc=raxml -trim ./138_Longest_Orthogroups_alignment/*einsi.faa

###################Tree construction
/home/anaconda3/bin/raxmlHPC-PTHREADS-SSE3 -m PROTGAMMAAUTO -f a -x 12345 -# 100 -T 30 -p 12345 -s output.phy -n EC_55species_Tree
