import sys
import re

mylines = []                                # Declare an empty list.

with open ("FOLDER_crameri/"+sys.argv[1]+"/"+sys.argv[1]+"_variable.txt", 'r') as myfile:    # Open lorem.txt for reading text.
    for myline in myfile:                   # For each line in the file,
        mylines.append(myline.rstrip('\n')) # strip newline and add to list.

#for element in mylines:                     # For each element in the list,
#    print(element)

o = open(sys.argv[1]+".pss", "w")

start = 0
end = 0
for i in range(len(mylines)):
    if mylines[i].startswith("Bayes Empirical Bayes (BEB) analysis"):
        start = i
    if mylines[i].startswith("The grid (see ternary graph for p0-p1)"):
        end = i

o.write("\n".join(mylines[start+1:end+1]))

o.close()
