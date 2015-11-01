#!/usr/bin/python -u
from __future__ import division
import sys, getopt
import re
import pysam
import pysamstats

def main(argv):
	if len(sys.argv) < 5:
		print 'usage: LVpicker.py -t <tumorfile> -n <normalfile> -r <reference> -b <bedfile> -o <outputpath>'
		sys.exit(1)
	try:
		opts, args = getopt.getopt(argv,"hr:b:t:n:o:", ["help","tumorfile=", "normalfile=", "reference=", "bedfile=", "outputpath="])
	except getopt.GetoptError:
		print 'usage: LVpicker.py -t <tumorfile> -n <normalfile> -r <reference> -b <bedfile> -o <outputpath>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
        	 	print 'usage: LVpicker.py -t <tumorfile> -n <normalfile> -r <reference> -b <bedfile> -o <outputpath>'
         		sys.exit()
		elif opt in ("-t", "--tumorfile"):
			tumorfile = arg
		elif opt in ("-n", "--normalfile"):
			normalfile = arg
		elif opt in ("-r", "--reference"):
			reference = arg
		elif opt in ("-b", "--bedfile"):
			bed = arg
      		elif opt in ("-o", "--outputpath"):
			outputpath = arg
			tempfile = outputpath+"/raw_somatic_varants.txt"
	return tumorfile, normalfile, reference, bed, outputpath

if __name__ == '__main__':
        (tumorfile, normalfile, reference, bed, outputpath) = main(sys.argv[1:])
        ref = pysam.FastaFile(reference)
	t_samfile = pysam.AlignmentFile(tumorfile)
	n_samfile = pysam.AlignmentFile(normalfile)
	tempfile = outputpath+"/raw_somatic_varants.txt"

	def process_reads(columns, columns_pos, refseq, filter_pos) :
		alt_A=0; alt_C=0; alt_G=0; alt_T=0; refcount=0
		foundReadName=[]; altReadPosE_f_A=[]; altReadPosE_r_A=[]; altReadPosE_f_G=[]; altReadPosE_r_G=[]; altReadPosE_f_C=[]; altReadPosE_r_C=[]; altReadPosE_f_T=[]; altReadPosE_r_T=[]; altReadPosS_f_A=[]; altReadPosS_r_A=[]; altReadPosS_f_G=[]; altReadPosS_r_G=[]; altReadPosS_f_C=[]; altReadPosS_r_C=[]; altReadPosS_f_T=[]; altReadPosS_r_T=[];
		pileupcolumn = columns; refseq = ref_seq
		if pileupcolumn.pos == columns_pos: 
			for pileupread in pileupcolumn.pileups:
				if pileupread.alignment.is_proper_pair and not pileupread.alignment.is_duplicate:
					if not pileupread.is_del and not pileupread.is_refskip:
						if pileupread.alignment.mapping_quality >= 30 and pileupread.alignment.query_qualities[pileupread.query_position] >= 30:
							altseq = pileupread.alignment.query_sequence[pileupread.query_position]
							try:
								foundReadName.index(pileupread.alignment.query_name)
							except ValueError:
								foundReadName.append(pileupread.alignment.query_name)
								if pileupread.alignment.query_sequence[pileupread.query_position] == ref_seq.upper():
									refcount += 1
								else:
									if altseq == "A":
										alt_A += 1
									elif altseq == "G":
										alt_G += 1
									elif altseq == "C":
										alt_C += 1
									elif altseq == "T":
										alt_T += 1
							pos = pileupread.alignment.qend - pileupread.query_position
							if filter_pos == "Y":
								if altseq == "A":
									if not pileupread.alignment.is_reverse:
										altReadPosE_f_A.append(pos)
										altReadPosS_f_A.append(pileupread.query_position)
									else:
										altReadPosE_r_A.append(pos)
										altReadPosS_r_A.append(pileupread.query_position)
								elif altseq == "G":
									if not pileupread.alignment.is_reverse:
										altReadPosE_f_G.append(pos)
										altReadPosS_f_G.append(pileupread.query_position)
									else:
										altReadPosE_r_G.append(pos)
										altReadPosS_r_G.append(pileupread.query_position)
								elif altseq == "C":
									if not pileupread.alignment.is_reverse:
										altReadPosE_f_C.append(pos)
										altReadPosS_f_C.append(pileupread.query_position)
									else:
										altReadPosE_r_C.append(pos)
										altReadPosS_r_C.append(pileupread.query_position)
								elif altseq == "T":
									if not pileupread.alignment.is_reverse:
										altReadPosE_f_T.append(pos)
										altReadPosS_f_T.append(pileupread.query_position)
									else:
										altReadPosE_r_T.append(pos)
										altReadPosS_r_T.append(pileupread.query_position)

			if filter_pos == "Y":
				return alt_A, alt_C, alt_G, alt_T, refcount, altReadPosE_f_A, altReadPosE_r_A, altReadPosE_f_G, altReadPosE_r_G, altReadPosE_f_C, altReadPosE_r_C, altReadPosE_f_T, altReadPosE_r_T, altReadPosS_f_A, altReadPosS_r_A, altReadPosS_f_G, altReadPosS_r_G, altReadPosS_f_C, altReadPosS_r_C, altReadPosS_f_T, altReadPosS_r_T
			if filter_pos == "N":
				return alt_A, alt_C, alt_G, alt_T, refcount

	def filter_germline(chr, t_columns_pos, ref_seq, alt_base):
		n_alt_A=0; n_alt_C=0; n_alt_G=0; n_alt_T=0; n_refcount=0
		for rec in pysamstats.stat_variation(n_samfile, ref, chrom=str(chr), start=int(t_columns_pos), end=int(t_columns_pos)+1):
			chr = rec['chrom']; pos = rec['pos']; ref_seq = rec['ref']
			if pos == t_columns_pos:
				if rec['reads_pp'] > 4:
					if float(rec['mismatches_pp']/rec['reads_pp']) <= 0.10:
						for n_columns in n_samfile.pileup(chr, int(t_columns_pos), int(t_columns_pos)+1, truncate=True):
							if n_columns.pos == t_columns_pos:
								(n_alt_A, n_alt_C, n_alt_G, n_alt_T, n_refcount) = process_reads(n_columns, n_columns.pos, ref_seq, "N");
					else:
						n_alt_A = rec['mismatches_pp']; n_alt_C = rec['mismatches_pp']; n_alt_G = rec['mismatches_pp']; n_alt_T = rec['mismatches_pp']; n_refcount = rec['reads_pp']-rec['mismatches_pp']
				else:
					n_alt_A=0; n_alt_G=0; n_alt_C=0; n_alt_T=0; n_refcount=rec['reads_pp']
 
		if alt_base == "A":
			return n_alt_A, n_refcount
              	if alt_base == "G":
                      	return n_alt_G, n_refcount
                if alt_base == "T":
                       	return n_alt_T, n_refcount
                if alt_base == "C":
                       	return n_alt_C, n_refcount

	foundReg=[];	raws=[]
	cov = 0
	ftemp = open(tempfile, 'w')	
	print >>ftemp, "chr\tpos\tref\talt\tt_ref\tt_alt\tn_ref\tn_alt\tjudge"

	for bedline in ( raw.strip().split() for raw in open(bed)):
 	     	for rec in pysamstats.stat_variation(t_samfile, ref, chrom=str(bedline[0]), start=int(bedline[1]), end=int(bedline[2])):
			if rec['reads_pp'] > 4:
				if rec['mismatches_pp'] > 2 and rec['insertions'] <= 3 and rec['deletions'] <= 3:
					chr = rec['chrom']; pos = rec['pos'];
					ref_seq = rec['ref']
					for t_columns in t_samfile.pileup(chr, int(pos), int(pos)+1, truncate=True):
						try:
							foundReg.index(str(chr)+"\t"+str(t_columns.pos))
                    				except ValueError:
                                			foundReg.append(str(chr)+"\t"+str(t_columns.pos))
							(t_alt_A, t_alt_C, t_alt_G, t_alt_T, t_refcount, altReadPosE_f_A, altReadPosE_r_A, altReadPosE_f_G, altReadPosE_r_G, altReadPosE_f_C, altReadPosE_r_C, altReadPosE_f_T, altReadPosE_r_T, altReadPosS_f_A, altReadPosS_r_A, altReadPosS_f_G, altReadPosS_r_G, altReadPosS_f_C, altReadPosS_r_C, altReadPosS_f_T, altReadPosS_r_T) = process_reads(t_columns, t_columns.pos, ref_seq, "Y")
							totalcount = t_refcount + t_alt_A + t_alt_G + t_alt_C + t_alt_T
							if totalcount > 4:
								if t_alt_A > 2 and t_alt_A/totalcount > 0.02:
									if (all(5 >= i for i in altReadPosS_f_A) and all(5 >= i for i in altReadPosS_r_A)) or (all(5 >= i for i in altReadPosE_f_A) and all(5 >= i for i in altReadPosE_r_A)):
										raw_calls = (chr+"\t"+str(pos)+"\t"+ref_seq.upper()+"\t"+"A"+"\t"+str(t_refcount)+"\t"+str(t_alt_A)+"\t"+"0"+"\t"+"0"+"\t"+"clustered_pos")
										print chr, pos, altReadPosS_r_A, altReadPosS_f_A, altReadPosE_f_A, altReadPosE_r_A, "A", "clustered_pos"
									else:
										(n_alt_A, n_refcount) = filter_germline(chr, t_columns.pos, ref_seq, "A")
										if n_alt_A < 2 and n_refcount >= 5:
											raw_calls = (chr+"\t"+str(pos)+"\t"+ref_seq.upper()+"\t"+"A"+"\t"+str(t_refcount)+"\t"+str(t_alt_A)+"\t"+str(n_refcount)+"\t"+str(n_alt_A)+"\t"+"pass_to_test")
										else:
											raw_calls = (chr+"\t"+str(pos)+"\t"+ref_seq.upper()+"\t"+"A"+"\t"+str(t_refcount)+"\t"+str(t_alt_A)+"\t"+str(n_refcount)+"\t"+str(n_alt_A)+"\t"+"possible_germline")
								elif t_alt_T > 2 and t_alt_T/totalcount > 0.02: 
									if (all(5 >= i for i in altReadPosS_f_T) and all(5 >= i for i in altReadPosS_r_T)) or (all(5 >= i for i in altReadPosE_f_T) and all(5 >= i for i in altReadPosE_r_T)):
										raw_calls = (chr+"\t"+str(pos)+"\t"+ref_seq.upper()+"\t"+"T"+"\t"+str(t_refcount)+"\t"+str(t_alt_T)+"\t"+"0"+"\t"+"0"+"\t"+"clustered_pos")
										print chr, pos, altReadPosS_r_T, altReadPosS_f_T, altReadPosE_f_T, altReadPosE_r_T, "T", "clustered_pos"
									else:
										(n_alt_T, n_refcount) = filter_germline(chr, t_columns.pos, ref_seq, "T")	
										if n_alt_T < 2 and n_refcount >= 5:
											raw_calls = (chr+"\t"+str(pos)+"\t"+ref_seq.upper()+"\t"+"T"+"\t"+str(t_refcount)+"\t"+str(t_alt_T)+"\t"+str(n_refcount)+"\t"+str(n_alt_T)+"\t"+"pass_to_test")
										else: 
											raw_calls = (chr+"\t"+str(pos)+"\t"+ref_seq.upper()+"\t"+"T"+"\t"+str(t_refcount)+"\t"+str(t_alt_T)+"\t"+str(n_refcount)+"\t"+str(n_alt_T)+"\t"+"possible_germline")
								elif t_alt_G > 2 and t_alt_G/totalcount > 0.02: 
									if (all(5 >= i for i in altReadPosS_f_G) and all(5 >= i for i in altReadPosS_r_G)) or (all(5 >= i for i in altReadPosE_f_G) and all(5 >= i for i in altReadPosE_r_G)):
										raw_calls = (chr+"\t"+str(pos)+"\t"+ref_seq.upper()+"\t"+"G"+"\t"+str(t_refcount)+"\t"+str(t_alt_G)+"\t"+"0"+"\t"+"0"+"\t"+"clustered_pos")
										print chr, pos, altReadPosS_r_G, altReadPosS_f_G, altReadPosE_f_G, altReadPosE_r_G, "G", "clustered_pos"
									else:
										(n_alt_G, n_refcount) = filter_germline(chr, t_columns.pos, ref_seq, "G")
										if n_alt_G < 2 and n_refcount >=5:
											raw_calls = (chr+"\t"+str(pos)+"\t"+ref_seq.upper()+"\t"+"G"+"\t"+str(t_refcount)+"\t"+str(t_alt_G)+"\t"+str(n_refcount)+"\t"+str(n_alt_G)+"\t"+"pass_to_test")
										else:
											raw_calls = (chr+"\t"+str(pos)+"\t"+ref_seq.upper()+"\t"+"G"+"\t"+str(t_refcount)+"\t"+str(t_alt_G)+"\t"+str(n_refcount)+"\t"+str(n_alt_G)+"\t"+"possible_germline")
								elif t_alt_C > 2 and t_alt_C/totalcount > 0.02:
									if (all(5 >= i for i in altReadPosS_f_C) and all(5 >= i for i in altReadPosS_r_C)) or (all(5 >= i for i in altReadPosE_f_C) and all(5 >= i for i in altReadPosE_r_C)):
							 			raw_calls = (chr+"\t"+str(pos)+"\t"+ref_seq.upper()+"\t"+"C"+"\t"+str(t_refcount)+"\t"+str(t_alt_C)+"\t"+"0"+"\t"+"0"+"\t"+"clustered_pos")
                                                                                print chr, pos, altReadPosS_r_C, altReadPosS_f_C, altReadPosE_f_C, altReadPosE_r_C, "G", "clustered_pos"
									else:
										(n_alt_C, n_refcount) = filter_germline(chr, t_columns.pos, ref_seq, "C")
										if n_alt_C < 2 and n_refcount >=5:
											raw_calls = (chr+"\t"+str(pos)+"\t"+ref_seq.upper()+"\t"+"C"+"\t"+str(t_refcount)+"\t"+str(t_alt_C)+"\t"+str(n_refcount)+"\t"+str(n_alt_C)+"\t"+"pass_to_test")
										else: 
											raw_calls = (chr+"\t"+str(pos)+"\t"+ref_seq.upper()+"\t"+"C"+"\t"+str(t_refcount)+"\t"+str(t_alt_C)+"\t"+str(n_refcount)+"\t"+str(n_alt_C)+"\t"+"possible_germline")
								else:
									raw_calls = (chr+"\t"+str(pos)+"\t"+ref_seq.upper()+"\t"+"N"+"\t"+str(t_refcount)+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t"+"low_mutant_reads")
							else:
								raw_calls = (chr+"\t"+str(pos)+"\t"+ref_seq.upper()+"\t"+"N"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t"+"low_total_reads")
							raws.append(raw_calls)
							print >>ftemp, raw_calls

	ftemp.close()

if __name__ == '__main__':
        main(sys.argv[1:])
