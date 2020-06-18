#!/usr/bin/env python

# author: Xiao-Yang Zhi; email: xyzhi@ynu.edu.cn
# please perpare the genome file and blastn db (by makeblastdb), before running this script.

import os
import sys

def usage():
	print("\npython3 seqHFC.py genome.fna blastn_database\n")

def remove_return(content):
	content_str = ""
	for i in content:
		if ">" in i:
			content_str+= i.strip()+">"
		else:
			content_str += i.strip()
	content_new = content_str.split(">")
	content_new.pop(0)
	content_new2 = list()
	for i in range(0,int(len(content_new)/2)):
		content_new2.append(">"+content_new[2*i]+"\n")
		content_new2.append(content_new[2*i+1]+"\n")
	return content_new2

def cuter(seq,stap=100):
	seqname = seq[0].strip()
	seqseq = seq[1].strip()
	unit = int(len(seqseq) / stap)
	group = list()
	for i in range(0,unit):
		group.append(seqname+"_"+str(stap*i)+"-"+str(stap*(i+1))+"\n")
		group.append(seqseq[stap*i:stap*(i+1)]+"\n")
	return group

def genome_fragmentation(genome):
	segments = list()
	for i in range(0,int(len(genome)/2)):
		seqname = genome[2*i]
		sequence = genome[2*i+1]
		cuttedseq = cuter((seqname,sequence))
		for c in cuttedseq:
			segments.append(c)
	return segments

def blast(query, blastdb):
	tmp_file_name = "/tmp/g"+str(os.getpid())
	file = open(tmp_file_name, 'w')
	for i in query:
		file.write(i)
	file.close()
	## blastn search
	hists = os.popen("blastn -outfmt '6 qacc qstart qend sallacc sstart send length mismatch gapopen pident' -num_threads 60 -evalue 0.00001 -max_target_seqs 1 -perc_identity 50 -query %s -db %s" % (tmp_file_name, blastdb)).readlines()
	os.remove(tmp_file_name)
	return hists

def hfc(hists, scaf_lens):
	scafolds = list()
	for hit in hists:
		sca = hit.split("_")[0]
		if sca not in scafolds:
			scafolds.append(sca)
	for sca in scafolds:
		# avoid one query with multiple hits
		query_deduplicates = list() 
		for hit in hists:
			if hit.split("_")[0] == sca:
				if hit.split("\t")[0] not in query_deduplicates:
					query_deduplicates.append(hit.split("\t")[0])
		sca_hit_count = len(query_deduplicates)
		percentage = float(sca_hit_count * 100) / scaf_lens[sca]
		print("%s\t%d\t%.3f" % (sca, sca_hit_count*100, percentage))

if __name__ == '__main__':
	try:
		genomeFastaFile = sys.argv[1:][0]
		blastnDatabase = sys.argv[1:][1]
	except IndexError:
		usage()
		exit(0)
	
	try:
		file = open(genomeFastaFile)
		genomeSeq = remove_return(file.readlines())
		file.close()
	except FileNotFoundError:
		print("Genome file (fasta) cannot be open!\n")
		usage()
		exit(0)
	#get scafolds' length
	scaf_lens = {}
	for i in range(0,int(len(genomeSeq)/2)):
		scaf_lens[genomeSeq[2*i].strip().replace(">","")] = len(genomeSeq[2*i+1].strip())
	#genome fragmentation
	gfragments = genome_fragmentation(genomeSeq)
	#blast search
	blastn_hists = blast(gfragments, blastnDatabase)
	#calculate percentage of homologous fragments
	hfc(blastn_hists,scaf_lens)



