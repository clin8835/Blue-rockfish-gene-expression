#!/usr/bin/env python
import sys

# usage 
if len(sys.argv) < 1:
        print "\n This program needs a file as an argument"

infle1 = sys.argv[1]
outfile = open("uniprotIDs_"+infle1,"w")

print infle1
print outfile

for line in open(infle1):
	if ('query' not in line):
		annot = (line.split()[1]).split("|")[1]
		#id = annot.spilt("|")[1]
		outfile.write(annot+"\n")

