import sys
import os
import re
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62
from os import environ, listdir, rename
from os.path import join, exists, abspath, pardir, basename, isdir, dirname, isfile

def open_file(filename):
	lines = []
	try:
		f = open(filename, 'r')
	except:
		print >> sys.stderr, "Error: Can't open file " + filename
		sys.exit()
	else:
		for line in f:
			lines.append(line.rstrip('\t\n'))
		f.close()	
		return lines

def write_file(filename, lines):

	output = open(filename, "w")
	for line in lines:
		output.write(line + "\n")
	output.close()

def open_aln_file(filename):

	lines = open_file(filename)

	q = ""
	t = ""
	for line in lines:
		if re.search("_q", line):
			label, seq = line.split("\t")
			q += seq
		if re.search("_t", line):
			label, seq = line.split("\t")
			t += seq

	n = len(q)
	i = 0
	q_aln = ""
	t_aln = ""
	while i < n:
		if q[i] != "-" and t[i] != "-":
			q_aln += q[i]
			t_aln += t[i]
		i += 1

	return q_aln, t_aln

def get_new_seq(t_seq, aln_seq):

	new_seq = ""
	count = 0
	for i in aln_seq: 
		if i != "-":
			new_seq += t_seq[count] 
			count += 1
		else:
			new_seq += "-"

	return new_seq

def read_fasta(filename):

	lines = open_file(filename)

	labels = []
	seqs = []
	 
	for line in lines:
		if line[0] == ">":
			labels.append(line)
		else:
			seqs.append(line)

	return labels, seqs  

query = sys.argv[1]
alns = {}
for file in listdir("."):
	if file.endswith(".txt"):
		path = join("./", file)
		q_aln, t_aln = open_aln_file(path)
		alns[file] = (q_aln, t_aln)

qseq = "%s.seq" % (query)
labels, seqs = read_fasta(qseq)
qq_aln = "%s_%s.txt" % (query, query)
seq1 = seqs[0]
#seq1 = alns[qq_aln][0]

print ">%s" % (query)
print seq1
for key in alns.keys():
	if key != qq_aln:
		q_aln, t_aln = alns[key]
		seq2 = q_aln

		print ">%s" % (key.rstrip(".txt").split("_")[1])

		#alignments = pairwise2.align.globalds(seq1, seq2, blosum62, -10, -0.1)
		alignments = pairwise2.align.localds(seq1, seq2, blosum62, -5, -0.5)
		seq1_aln, seq2_aln, score, pos1, pos2 = alignments[0]

		new_t_aln = get_new_seq(t_aln, seq2_aln)

		print new_t_aln




