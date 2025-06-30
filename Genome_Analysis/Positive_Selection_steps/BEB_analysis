###################Bayes Empirical Bayes (BEB) analysis to identify the genes with positively selected sites (PSS) using the CodeML output directories
cat list_genes | nohup parallel -j 150 'python extractsites_pss.py {}' &
ls *.pss | nohup parallel -j 150 "grep -v "Positive" {} > {}.positive" &
ls *.positive | nohup parallel -j 150 "sed -i '/^$/d' {}" &
ls *.positive | nohup parallel -j 150 "sed 's/^[ \t]*//' {} > {}.space" &
ls *.space | nohup parallel -j 150 "grep -v "The grid" {} > {}.grid" &
ls *.grid | nohup parallel -j 150 "sed 's/ /\t/g' {} > {}.mod" &
cat *.grid.mod > EC_genes_with_PSS.txt

###################Extracting the genes with PSS (>95% probability)
awk -F'\t' '{if ($3 > 0.95 ) {print $1"\t"$2"\t"$3}}' EC_genes_with_PSS.txt > EC_selected_genes_with_PSS.txt
