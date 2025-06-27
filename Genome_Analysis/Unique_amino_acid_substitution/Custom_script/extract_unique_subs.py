import re
import sys
d = {}
species = ["Bombyx#mori", "Danaus#plexippus", "Heliconius#melpomene", "Melitaea#cinxia","eumeta2#variegata","Eumeta1#crameri"]
dir1 = "/run/media/abhisek/HDD3/WB_data/Orthofinder_analysis/MSA/OrthoFinder/Results_Oct13/kinfin_analysis/MSA_orthogroups/crameri/mafft_aligned_orthogroups_single_crameri/"
for f1 in open(dir1 + sys.argv[1], 'r'):
  f1 = f1.rstrip()
  if (">" in f1):
    ann = f1.replace(">", "")
    ann_a = ann.split("#")
    ann1 = ann_a[0] + "#" + ann_a[1]
  else:
    d[ann1] = f1

d2 = {}
humanseq = list(d[species[0]])
for i in range(0, len(humanseq)):
  d2[i] = 1

for i in range(0, len(species)-1):
  seq = list(d[species[i]])
  for j in range(0, len(seq)):
    if (seq[j] == "-"):
      a = j-10
      if (a < 0):
        a = 0
      #sys.stdout.write(str(a) + " ")
      for k in range(a, (j+10)):
        #sys.stdout.write(str(k) + " ")
        d2[k] = 0
    if (d2[j] == 1):
      if (seq[j] != humanseq[j]):
        d2[j] = 0

i = len(species)-1
seq = list(d[species[i]])
for j in range(0, len(seq)):
  if (seq[j] == "-"):
    a = j-10
    if (a < 0):
      a = 0
    for k in range(a, j+10):
      d2[k] = 0
for j in range(0, len(seq)):
  if ((d2[j] == 1) and (seq[j] != humanseq[j]) and (seq[j] != 'X')):
    pos = j + 1
    for k in range(0, pos):
      if (humanseq[k] == "-"):
        pos = pos - 1
    #sys.stdout.write(humanseq[j] + str(pos) + seq[j] + ";")
    a = [item.replace("-", "") for item in humanseq]
    a2 = "".join(a)
    target1 = open('crameri_subs_10/' + sys.argv[1] + '.fasta', 'w')
    target1.write(">" + sys.argv[1] + "\n" + a2 + "\n")
    target1.close()
    target2 = open('crameri_subs_10/' + sys.argv[1] + '.subst', 'a')
    target2.write(humanseq[j] + str(pos) + seq[j] + "\n")
    target2.close()
