import sys
import os
from shutil import copyfile
from os import environ, listdir, rename, mkdir
from os.path import join, exists, abspath, pardir, basename, isdir, dirname, isfile


if not isdir("./blocks_ali"):
	mkdir("./blocks_ali")
if not isdir("./best_superposition"):
	mkdir("./best_superposition")
if not isdir("./matrix_aln"):
	mkdir("./matrix_aln")

for folder in listdir("."):
	path = join(".", folder)
	
	if isdir(path) and path != "./blocks_ali" and path != "./best_superposition" and path != "./matrix_aln":
		tok = basename(path).split("-")
		q, t = tok[0], tok[1]
		
		src_file1 = join(path, "blocks_ali.txt")
		dst_file1 = join("./blocks_ali", "%s_%s.txt" % (q,t)) 

		src_file2 = join(path, "best_combination/best_alignment.p1m")
		dst_file2 = join("./best_superposition", "%s_%s.p1m" % (q,t))

		src_file3 = join(path, "blocks_selected.png")
		dst_file3 = join("./matrix_aln", "%s_%s.png" % (q,t))

		print src_file1, dst_file1
		print src_file2, dst_file2
		print src_file3, dst_file3
		copyfile(src_file1, dst_file1)
		copyfile(src_file2, dst_file2)
		copyfile(src_file3, dst_file3)


