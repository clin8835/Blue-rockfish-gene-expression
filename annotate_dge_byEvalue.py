#!/usr/bin/env python
import sys
import numpy as np


##################################################################################################
# usage 
#         This program needs a file as an argument
#     
#	Usage  -files <dge File first> <then annotation File>
#
#	example:  python annotate_dge_byEvalue.py -files gopher_ConDo_diffExpr.P0.01_C1.matrix.log2.centered.dat ../trinity_2sample_run_20180528.Trinity.swisprot.fasta

# Parse args.
for i in range(len(sys.argv)):
        if sys.argv[i]== "-files":
                infle1 = sys.argv[i+1]
                infle2 = sys.argv[i+2]


outfile = open(infle1+"_annot_byEval.txt","w")
outfile.write("query id\tsubject id\t%_identity\talignment_length\tmismatches\tgap_opens\tq.start\tq.end\ts.start\ts.end\tevalue\tbit_score\n")
print infle1
print infle2
#print outfile

# if the contig is the same as the one before it

prevContig=""
for line in open (infle1):
        #annot = (line.split()[1]).split("|")[1]
        #id = annot.spilt("|")[1]
	c = list()
	dgeContig = line.split()[0]
	minEval=[]
        for row in open(infle2):
		if "query" in row:
			continue
		isoform = (row.split()[0]).split("_")[4]
		anotContig = "_".join((row.split()[0]).split("_")[:4])
		if (dgeContig == anotContig):
			c.append(row.split())

########################################################################################################################################################
# query id        subject id      %_identity      alignment_length        mismatches      gap_opens    q.start q. end  s.start s.end   evalue  bit_score
#	0	1			2		3			4		5	6	7	8	9	10	11  
#	0			1			2	3	4	5	6	7	8	9	10
# TRINITY_DN81774_c0_g1_i1	sp|P02676|FIBB_BOVIN	51.3	76	35	1	232	5	49	122	1.0e-13	77.0

# if there is more than one isoform of the same DGE gene select the isoform with the lowest e-value,
#	longest contig length, and if there is still multiple isoforms check if they blast to different genes
#		1) if they are the same gene, choose lowest isoform number
#		2) if they are different genes, they are selected by longest alignment length
#		3) we should probably include %identity into this analysis, currently it is not included (7/30/16)
########################################################################################################################################################

	if (len(c)>1):
		print "hit list c:", "\n", c
		evalues = list()
		#pID = list()
		for i in range(len(c)):
			#print "%ID", c[i][2]		
			#print  "evalue", c[i][10]
			evalues.append(float(c[i][10]))
			#pID.append(c[i][2])
		minE = min(evalues)
		# if the minE value is _.0 it removes the trailing 0, this put its back to compare correctly to the c list
		fixE = str(minE).split('e')
		if ("." not in fixE[0]):
			num = str(fixE[0]+".0e"+fixE[1])
			minE=num
		#maxpID = max(pID)
		print "min evalue:", minE
		#print "max %ID:", max(pID)
		
		bestHit= list()

		# go through the list of contig isoforms and put the one(s) that equal the lowest e-value into another list
		for e in range(len(c)):
			if (str(minE) in c[e][10]):
				bestHit.append(c[e][:])
		# if there is only one match to the lowest e-val, write it to the outfile
		if (len(bestHit)==1):
			print "bestHit:", "length:", len(bestHit), "\n", "\t".join(bestHit[0][:])
			outfile.write("\t".join(bestHit[0][:])+"\n")
		# if there is more than one contig isoform with the lowest evalue, compare them by alignment length
		elif (len(bestHit)>1):
			print "best Hit Length:", len(bestHit), "\n", bestHit
			#align_len = list()
			pID = list()
			betterHit = list()
			# if the gene hits are exactly the same, write the one with the lowest isoform number
			for b in range(len(bestHit)):
				#if (bestHit[b][2] == bestHit[b+1][2]):	
				#	betterHit.append(bestHit[b][:])
				pID.append(bestHit[b][2])
			max_pID = max(pID)
			print "max_pID:", max_pID
			for p in range(len(bestHit)):
				if (bestHit[p][2] == str(max_pID)):
					betterHit.append(bestHit[p][:])
			# in the event they are not the same gene, 
			# compare by the length of the alignment
			#if (len(betterHit)==0):
				# make a list of the alignment length's
			#	for n in range(len(bestHit)):
			#		align_len.append(bestHit[n][3])
				# find the max length
			#	longAlign = max(align_len)
				# put the longest hit in the betterHit list
			#	for z in range(len(bestHit)):
			#		if (bestHit[z][3] == str(longAlign)):
			#			betterHit.append(bestHit[z][:])
			# write that better hit to the outfile
			print "betterHit:", betterHit
			outfile.write("\t".join(betterHit[0][:])+"\n")
	# if there is only one isoform in the annotation file, then just write that to the out file
	elif (len(c)==1):
		outfile.write("\t".join(c[0][:])+"\n")
		outstring =  c[0][0].split()[:3] + c[0][1:]
		print "One Hit c join:", "\t".join(c[0][:])




		#isoform = (row.split()[0]).split("_")[4]
		#if dgeContig :
		#	minEval.append(row)
		#	print "duplicate"	
		
		#if line.split()[0] in row:
                #if (isoform == "i1") and (line.split()[0] in row):
		#	prevContig = (row.split()[0]).split("_")[:4]
#			rowString = "\t".join(row.split()[1:])
#			list = row.split()
#                        outfile.write(row)
#			#outfile.write(line.split()[0]+"\t"+rowString+"\n")
#                        print line.split()[0], rowString
#                        #print line.split()[0], "\t".join(row.split()[1:])
#			print row
outfile.close()


