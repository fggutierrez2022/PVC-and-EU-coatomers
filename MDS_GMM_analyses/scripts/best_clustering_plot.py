"""
best_clustering_plot.py

This script plots the relative distance between a group of protein strutures
based on structural comparisons, estimating from the data the best clustering 
of proteins.

Autor: Fernando Gutiérrez
email:	fgutierrez1988@gmail.com
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from os.path import isfile, basename
import pandas as pd
import collections
from operator import itemgetter
from sklearn.manifold import MDS
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_samples, silhouette_score

# setting the seaborn style plot
sns.set(style="dark")

def open_file(filename):
	"""Opens a text file and extracts their lines"""
	
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
	"""Writes a new file"""

	output = open(filename, "w")
	for line in lines:
		output.write(line + "\n")
	output.close()

def extract_scores(results_path):
	"""Extract the results of the comparisons performed with MOMA2"""

	lines = open_file(results_path)
	
	scores = {}
	eq_qt = []
	labels = []
	for line in lines[1:]:
		tok = line.split("\t")
		# extract the alignment length, relative similarity, structural overlap, b score, the percentage of sequence identity and the ponderate RMSD 
		eq, sr, so, b, seqi, rmsdp = int(tok[2]), float(tok[3]), float(tok[4]), float(tok[5]), float(tok[6]), float(tok[7])
		qlabel, tlabel = tok[8], tok[9]
		
		if qlabel not in labels:
			labels.append(qlabel)
		if tlabel not in labels:
			labels.append(tlabel)
		
		scores[(qlabel, tlabel)] = (eq, sr, so, b, rmsdp)
		scores[(tlabel, qlabel)] = (eq, sr, so, b, rmsdp)

		if qlabel != tlabel:
			eq_qt.append(eq)
	
	max_eq = max(eq_qt)
	
	return labels, scores, max_eq

def dist_bscore(qname, tname, scores):
	"""Calculate a distance metric from the bscores values"""

	b_qt = scores[(qname, tname)][3]
	b_qq = scores[(qname, qname)][3]
	b_tt = scores[(tname, tname)][3]
	
	d = 1 - (2 * (b_qt)/(b_qq + b_tt))
	
	return float("%.3f" % (d))
	

def create_diff_tbl(labels, scores, max_eq, score_t="bscore"):
	"""Create a distance table based on the similarity score selected"""

	n = len(labels)
	
	dt = {}
	i = 0
	while i < n:
		j = 0
		name = labels[i]
		distances = [name]
		while j < n:
			qname = labels[i]
			tname = labels[j]
			if score_t == "bscore": # bscore
				d = dist_bscore(qname, tname, scores)
				distances.append(str(d))
			elif score_t == "sr": # relative similarity
				d = 1 - scores[(qname, tname)][1]
				distances.append(str(d))
			else:
				sys.exit()
			j += 1
		dt[name] = distances
		i += 1
	
	return dt

def create_csv(dt, labels, filename):
	"""Create a csv file from the distance table"""
	
	csv_lines = []
	
	line = "\t".join([''] + labels)
	csv_lines.append(line)
	for label in labels:
		line = "\t".join(dt[label])
		csv_lines.append(line)
	write_file(filename, csv_lines)

def plot_gmm(gmm, X, labels, title, svg_path):
	"""Plot the groups predicted by the GMM method"""
	
	plt.figure(figsize=(15,8))
	clusters = gmm.fit(X).predict(X)
	plt.scatter(X[:, 0], X[:, 1], c=clusters, s=40, cmap='viridis', zorder=2)

	plt.title(title)
	for label, x, y in zip(labels, X[:, 0], X[:, 1]):
		plt.annotate(label, xy = (x, y), xytext = (-15,8), textcoords = 'offset points')
	plt.savefig(svg_path, dpi=350, format="svg")
	plt.close()

def running_MDS(X, seed):
	"""Run the MDS method with 2 dimensions"""
	
	mds_sklearn = MDS(n_components=2, dissimilarity='precomputed', random_state=seed, eps=1e-5, max_iter=1000)
	x_sklearn = mds_sklearn.fit_transform(X)
	return mds_sklearn.stress_, x_sklearn

def getting_seeds(ini, ter):
	"""Generate a list of seeds"""

	seeds = []  
	for i in xrange(ini, ter, 1):
		seeds.append(i)

	return seeds

def check_clustering_with_outgroup(outgroup, labels, clf_groups):
	"""Check if the outgroup is well classified in an only group"""

	n = len(labels)
	i = 0
	sel_group = None
	while i < n:
		label = labels[i]
		group = clf_groups[i]
		
		# get the classification of the outgroup by clf
		if label == outgroup:
			sel_group = group
		i += 1
	
	# if the outgroup was classified in an only group, everything is ok
	counter = collections.Counter(clf_groups)
	if counter[sel_group] == 1:
		return True
	else:
		return False

def running_GMM(X, seed, labels, outgroup=None):
	"""Run the GMM method with 3 to 9 numbers of mixture components"""

	range_n_clusters = [3, 4, 5, 6, 7, 8, 9]
	
	n_clusters_sel = 0
	sil_avg_max = -1.0

	for n_clusters in range_n_clusters:
	 	# Fit a Gaussian mixture with the coordinates determined with the MDS method
    		clusterer = GaussianMixture(n_components=n_clusters, covariance_type='full', random_state=seed)
    		cluster_labels = clusterer.fit_predict(X)
    		silhouette_avg = silhouette_score(X, cluster_labels)
    		
    		if outgroup != None:
	    		if not check_clustering_with_outgroup(outgroup, labels, cluster_labels):
	    			continue
	    	
    		# Select the best clustering with the highest silhouette score
    		if silhouette_avg > sil_avg_max:
    			sil_avg_max = silhouette_avg
    			n_clusters_sel = n_clusters
    	
    	return n_clusters_sel, sil_avg_max

def main():

	# checking input
	if len(sys.argv) != 4 and len(sys.argv) != 5:
		print "Usage: python best_clustering_plot.py <file with results> <bscore|sr> <title> <outgroup>"
		sys.exit()

	results_path = sys.argv[1]
	score_option = sys.argv[2]
	title = sys.argv[3]
	
	outgroup = None
	if len(sys.argv) == 5:
		outgroup = sys.argv[4]
	
	if not isfile(results_path):
		print "%s does not exists!!!" % (results_path)
		sys.exit()
		
	if score_option != "bscore" and score_option != "sr":
		print "%s is not valid!!!" % score_option
		print "use 'bscore' or 'sr'..."
		sys.exit()

	# extract scores to create a distance table according to the type of score selected
	csv_filename = "diff_tbl_%s.csv" % (score_option)
	labels, scores, max_eq = extract_scores(results_path)
	dt = create_diff_tbl(labels, scores, max_eq, score_option)
	create_csv(dt, labels, csv_filename)
	
	# distance table to dataframe
	df = pd.read_csv(csv_filename, sep="\t", header = 0, index_col = 0)
	X = df.to_numpy()

	# determine the best MDS plot
	seeds = getting_seeds(3000, 4000)
	print "Searching the best MDS plot with a minor stress value:"
	results_MDS = []
	counter = 1
	for seed in seeds:
		stress, _ = running_MDS(X, seed)
		print counter, seed, stress
		results_MDS.append((seed, stress))
		counter += 1
	
	results_MDS_sorted = sorted(results_MDS, key=itemgetter(1), reverse=False)
	seed_MDS_sel, stress_sel = results_MDS_sorted[0]
	
	_ , x_sklearn = running_MDS(X, seed_MDS_sel)
	
	# determine the best clusters
	seeds = getting_seeds(1, 1000)
	results_GMM = []
	counter = 1
	for clf_seed in seeds:
		n_clusters, sil = running_GMM(x_sklearn, clf_seed, labels, outgroup)
		print counter, clf_seed, sil, n_clusters
		results_GMM.append((clf_seed, sil, n_clusters))
		counter += 1

	# Print the best five plots
	print "\nResults"
	print "MDS (seed, stress):", seed_MDS_sel, "%.5f" % stress_sel
	
	results_GMM_sorted = sorted(results_GMM, key=itemgetter(1,2), reverse=True)
	
	counter = 1
	bef_group = -1.0
	for res in results_GMM_sorted:
		seed_sel, sil_sel, n_sel = res
		if counter == 6:
			break
		if sil_sel != -1.0 and n_sel != bef_group:
			gmm = GaussianMixture(n_components=n_sel, covariance_type='full', random_state=seed_sel)
			plot_gmm(gmm, x_sklearn, labels, "%s (MDS and GMM, n: %d, sil: %.4f)" % (title, n_sel, sil_sel), "%s_%d.svg" % (basename(results_path).split(".")[0],counter))
			print "best clustering GMM (n groups, seed, avg_sil): %d" % n_sel, seed_sel, sil_sel
			bef_group = n_sel
			counter += 1
	
if __name__ == "__main__":
	main()
