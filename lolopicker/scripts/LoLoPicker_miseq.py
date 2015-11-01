#!/usr/bin/python -u
from __future__ import division
import pysam
import pysamstats
import re
import numpy
from scipy import stats
from scipy.stats import binom
from scipy.cluster.vq import kmeans2, whiten

sample=[]
inputlist = open("/gs/project/vdu-032-aa/jzhang/DNET/Miseq_samplelist.txt")
for line in ( raw.strip().split() for raw in inputlist ):
	sample.append(line)
inputlist.close()

def process_stats(chr, pos, ref_seq, alt):
	print chr, pos, alt
	c_alf_array=[]; c_INFO_array=[]; judge_array=[]
	for bam in sample:
		c_ref_cov=0; c_variant_cov=0;
		samfile = pysam.AlignmentFile(bam[0])
		if "exon" in bam[0]:
                        chr = "chr8"
                else:
                        chr = "8"
		for pileupcolumn in samfile.pileup(chr, pos, pos+1, truncate=True):
			if pileupcolumn.pos ==  pos:
				for pileupread in pileupcolumn.pileups:
					if pileupread.alignment.mapping_quality >= 30 and pileupread.alignment.query_qualities[pileupread.query_position] >= 30:
						if pileupread.alignment.query_sequence[pileupread.query_position] == ref_seq:
							 c_ref_cov += 1
						elif pileupread.alignment.query_sequence[pileupread.query_position] == alt:
							c_variant_cov += 1
				c_total = c_ref_cov+c_variant_cov
				if c_total >= 10:
					c_alf = c_variant_cov/c_total
					allINFO = c_ref_cov, c_variant_cov, c_alf, bam[1]
					c_INFO_array.append(allINFO)
					c_alf_array.append(c_alf)
	
	if all(i < 0.1 for i in c_alf_array):
		return "NA"
	elif all(i >= 0.1 for i in c_alf_array):
		return c_INFO_array	

	else:
		output=[]
		while True:
			km_results = kmeans2(c_alf_array, 2, iter=10, thresh=0)
			km = km_results[0]
			cluster1 = km[0]; cluster2 = km[1]
			if cluster1 > cluster2:
				if cluster1 < 1 and cluster2 > 0:
					break
			if cluster1 < cluster2:
				if cluster2 < 1 and cluster1 > 0:
					break
		if cluster1 > cluster2:
			j=0
			for f in c_alf_array:
				if km_results[1][j] == 0:
					output.append(c_INFO_array[j])
				j += 1
		if cluster1 < cluster2:
			j=0
			for f in c_alf_array:
				if km_results[1][j] == 1:
                                	output.append(c_INFO_array[j])
				j += 1
		return output


bed = "/gs/project/vdu-032-aa/jzhang/DNET/FGFR1.bed"

ref = pysam.FastaFile("/gs/project/vdu-032-aa/ngs/repository/references/lib/hg19_broad_reference/Homo_sapiens_assembly19.fasta")

for bedline in ( raw.strip().split() for raw in open(bed)):
	foundReg=[]
	for bam in sample:
		samfile = pysam.AlignmentFile(bam[0])
		if "exon" in bam[0]:
			ref = pysam.FastaFile("/gs/project/vdu-032-aa/ngs/repository/references/lib/hg19_chr_files/hg19_wRandomsNew.fa")
			chr = "chr8"
		else: 
			ref = pysam.FastaFile("/gs/project/vdu-032-aa/ngs/repository/references/lib/hg19_broad_reference/Homo_sapiens_assembly19.fasta")
			chr = "8"
		for rec in pysamstats.stat_variation(samfile, ref, chrom=chr, start=int(bedline[1]), end=int(bedline[2])):
			if rec['reads_pp'] >= 100:
				if rec['mismatches_pp']/rec['reads_pp'] > 0.01:
					chr = rec['chrom']; pos = rec['pos']; ref_seq = rec['ref']
					variant_A_cov=0; variant_G_cov=0; variant_C_cov=0; variant_T_cov=0
					ref_cov = 0
					for pileupcolumn in samfile.pileup(chr, pos, pos+1, truncate=True):		
						if pileupcolumn.pos ==  pos:
							print rec['mismatches_pp'], rec['reads_pp'], bam[0]
							for pileupread in pileupcolumn.pileups:
								if pileupread.alignment.query_qualities[pileupread.query_position] >= 30:
									if pileupread.alignment.query_sequence[pileupread.query_position] == ref_seq:
										ref_cov += 1
									elif pileupread.alignment.query_sequence[pileupread.query_position] == "A":
										variant_A_cov += 1
									elif pileupread.alignment.query_sequence[pileupread.query_position] == "T":
										variant_T_cov += 1
									elif pileupread.alignment.query_sequence[pileupread.query_position] == "C":
										variant_C_cov += 1
									elif pileupread.alignment.query_sequence[pileupread.query_position] == "G":
										variant_G_cov += 1	
							total = variant_A_cov +  variant_T_cov + variant_C_cov + variant_G_cov + ref_cov
							print total, variant_C_cov, bam[0]
							if total >= 10:
								t_alf_A = variant_A_cov/total
								t_alf_T = variant_T_cov/total
								t_alf_G = variant_G_cov/total
								t_alf_C = variant_C_cov/total
						
								if t_alf_A >= 0.01:
									try:
										foundReg.index(str(chr)+"\t"+str(pileupcolumn.pos)+"\t"+ "A")
                                                			except ValueError:
                                                        			foundReg.append(str(chr)+"\t"+str(pileupcolumn.pos)+"\t"+ "A")
										results = process_stats(chr, pos, ref_seq, "A")
										if results != "NA":
											print chr, pos, "A", results
								elif t_alf_G >= 0.01:
									try:
										foundReg.index(str(chr)+"\t"+str(pileupcolumn.pos)+"\t"+ "G")
									except ValueError:
										foundReg.append(str(chr)+"\t"+str(pileupcolumn.pos)+"\t"+ "G")
										results = process_stats(chr, pos, ref_seq, "G") 
										if results != "NA":
											print chr, pos, "G", results
								elif t_alf_T >= 0.01:
									try:
										foundReg.index(str(chr)+"\t"+str(pileupcolumn.pos)+"\t"+ "T")
									except ValueError:
										foundReg.append(str(chr)+"\t"+str(pileupcolumn.pos)+"\t"+ "T")
										results = process_stats(chr, pos, ref_seq, "T")
										if results != "NA":
											print chr, pos, "T", results
								elif t_alf_C >= 0.01:
									try:	
										foundReg.index(str(chr)+"\t"+str(pileupcolumn.pos)+"\t"+ "C")
									except ValueError:
										foundReg.append(str(chr)+"\t"+str(pileupcolumn.pos)+"\t"+ "C")
										results = process_stats(chr, pos, ref_seq, "C")
										if results != "NA":
											print chr, pos, "C", results
										

									
									

