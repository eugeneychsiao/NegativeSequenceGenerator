import numpy as np
import pandas as pd
from pybedtools import BedTool
import os
import argparse
import random

#References
hg19_without_N = "reference/hg19_rm_Ns.bed"
hg38_without_N = "reference/hg38_rm_Ns.bed"
hg19_without_N_repeat = "reference/hg19_rm_repeat_Ns.bed"
hg38_without_N_repeat = "reference/hg38_rm_repeat_Ns.bed"
hg19_bed_file = "reference/hg19.fa.bed"
hg38_bed_file = "reference/hg38.fa.bed"
hg19_genome = "reference/human.hg19.genome"
hg38_genome = "reference/human.hg38.genome"
hg19_fasta = "../../../Downloads/hg19.fa"
hg38_fasta = "../../../Downloads/hg38.fa"




def rewrite_bed(filename, sequence_length):
	f = BedTool(filename)
	positive_peak_counter = 0
	kept_counter = 0
	with open('positive.bed', 'w') as g:
		for line in f:
			positive_peak_counter+=1
			mid_point = int(int((int(line[1]) + int(line[2]))) / 2)
			original_len = int(line[2]) - int(line[1])
			if mid_point - int(sequence_length / 2) < 0:
				pass
			else:
				if sequence_length % 2 == 0:
					line[1] = mid_point - int(sequence_length / 2)
					line [2] = mid_point + int(sequence_length / 2)
					g.write('\t'.join(line) + "\n") 
					kept_counter += 1
				else:
					line[1] = mid_point - int(sequence_length / 2)
					line [2] = mid_point + int(sequence_length / 2) + 1
					g.write('\t'.join(line) + "\n")
					kept_counter += 1

	return positive_peak_counter, kept_counter


# rewrite_bed("a.bed" , 1000)


def create_bed_and_fasta(hg =38):
	if os.path.exists('positive_fasta.fa') is False:
		if hg == 19:
			os.system('bedtools getfasta -fi %s -bed positive.bed -fo positive_fasta.fa' %(hg19_fasta))        #fasta file with positive peak sequences
			# print("19")
		else: 
			os.system('bedtools getfasta -fi %s -bed positive.bed -fo positive_fasta.fa' %(hg38_fasta)) 
			# print("38")   
	if os.path.exists('all_excluded.bed') is False:
		if hg ==19:
			os.system('bedtools subtract -a %s -b %s > N_regions.bed' %(hg19_bed_file, hg19_without_N))
			os.system('multiIntersectBed -i positive.bed N_regions.bed > all_excluded.bed') 
			os.system('bedtools sort -i all_excluded.bed > sorted_all_excluded.bed')
			os.system('bedtools merge -i sorted_all_excluded.bed > merged_all_excluded.bed')
			# print("19")
		else: 
			os.system('bedtools subtract -a %s -b %s > N_regions.bed' %(hg38_bed_file, hg38_without_N))
			os.system('bedtools subtract -a %s -b %s > N_regions.bed' %("N_regions.bed", "reference/RXRA.bed"))
			os.system('multiIntersectBed -i positive.bed N_regions.bed > all_excluded.bed') 
			os.system('bedtools sort -i all_excluded.bed > sorted_all_excluded.bed')          
			os.system('bedtools merge -i sorted_all_excluded.bed > merged_all_excluded.bed')
			# print("38")

# create_bed_and_fasta(19)

def create_negative_bed( pos_num, ratio , sequence_length,hg=38):
	total_negative_samples = round(pos_num * float(ratio))
	# print(total_negative_samples)
	if os.path.exists('random_negative.bed') is False:
		if hg == 19:
			os.system('bedtools random -g %s -n %s -l %s > random_negative.bed' %( hg19_genome,total_negative_samples, sequence_length))        
			# print("19")
		else: 
			os.system('bedtools random -g %s -n %s -l %s > random_negative.bed' %( hg38_genome,total_negative_samples, sequence_length)) 
			# print("38")   

	if os.path.exists('negative.bed') is False:
		if hg ==19:
			os.system('shuffleBed -i random_negative.bed -g %s -excl merged_all_excluded.bed > negative.bed' %(hg19_genome))
			# print("19")
		else: 
			os.system('shuffleBed -i random_negative.bed -g %s -excl merged_all_excluded.bed > negative.bed' %(hg38_genome))
			# print("38")
	if os.path.exists('negative_fasta.fa') is False:
		if hg == 19:
			os.system('bedtools getfasta -fi %s -bed negative.bed -fo negative_fasta.fa' %(hg19_fasta))        #fasta file with positive peak sequences
			# print("19")
		else: 
			os.system('bedtools getfasta -fi %s -bed negative.bed -fo negative_fasta.fa'%(hg38_fasta)) 
			# print("38") 


# create_negative_bed( 10 , 10, 100)


def install_genome_fa(hg=38):
	if hg ==19:
		if os.path.exists('reference/hg19.fa') is False:
			print("Downloading hg19 fasta file.")
			os.system('wget -P reference/ http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz')
			os.system('gzip -d reference/hg19.fa.gz')


	else: 
		if os.path.exists('reference/hg38.fa') is False:
			print("Downloading hg38 fasta file.")
			os.system('wget -P reference/ http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz')
			os.system('gzip -d reference/hg38.fa.gz')


def summary(filename, pos_count, kept_count, hg, sequencelen):
	negative_counter = 0
	f = BedTool(filename)
	with open('Summary.txt', 'w') as g:
		for line in f:
			negative_counter += 1
		g.write("Genome version: hg%s \n" %(hg))
		g.write("Sequence length: %s \n" %(sequencelen))
		g.write("%s (of %s) positive data-points kept \n" %(kept_count, pos_count))
		g.write("%s negative data-points generated" %(negative_counter))






if __name__ == '__main__':
	parser=argparse.ArgumentParser()
	parser.add_argument("filename",help="Enter .bed file")
	parser.add_argument("genome",help="Enter hg19 or hg38")
	parser.add_argument("seqlen",help="Enter sequence length")
	parser.add_argument("ratio",help="Enter positive/negative ratio")
	args=parser.parse_args()


	pos_count,kept_count = rewrite_bed(args.filename, int(args.seqlen))
	if args.genome=='hg19':
		# install_genome_fa(19)
		create_bed_and_fasta(hg = 19)
		create_negative_bed(int(pos_count), args.ratio, int(args.seqlen), hg = 19)
		summary('negative.bed', pos_count, kept_count, 19, args.seqlen)

	else: 
		# install_genome_fa(38)
		create_bed_and_fasta(hg = 38)
		create_negative_bed(int(pos_count), args.ratio, int(args.seqlen))
		summary('negative.bed', pos_count, kept_count, 38, args.seqlen)



	os.system('rm -r positive.bed')
	os.system('rm -r N_regions.bed')
	os.system('rm -r negative.bed')
	os.system('rm -r random_negative.bed')
	os.system('rm -r all_excluded.bed')
	os.system('rm -r merged_all_excluded.bed')
	os.system('rm -r sorted_all_excluded.bed')





