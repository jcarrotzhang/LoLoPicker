#!/usr/bin/python -u
from __future__ import division
import sys, getopt
import re
import pysam
import pysamstats

def main(argv):
	base_quality = 20
	mapping_quality = 20
	tumor_alt_cutoff = 1
	tumor_alf_cutoff = 0
	normal_alt_cutoff = 3
	normal_alf_cutoff = 0
	phasing_mode = "Y"

	if len(sys.argv) < 5:
		print 'usage: python LoLoPicker_somatic.py -t <tumor.bam> -n <normal.bam> -r <reference.fa> -b <intervals.bed> -o <outputpath>'
		sys.exit(1)
	try:
		opts, args = getopt.getopt(argv,"ht:n:r:b:o:B:m:f:c:F:C:P", ["help","tumorfile=", "normalfile=", "reference=", "bedfile=", "outputpath=", "basequality=", "mappingquality=", "normalalteredratio=", "normalalteredreads=", "tumoralteredratio=", "tumoralteredreads=", "phasing_mode="])
	except getopt.error:
		print 'usage: python LoLoPicker_somatic.py -t <tumor.bam> -n <normal.bam> -r <reference.fa> -b <intervals.bed> -o <outputpath>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
        	 	print 'usage: LoLoPicker_somatic.py -t <tumor.bam> -n <normal.bam> -r <reference.fa> -b <intervals.bed> -o <outputpath>'
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
		elif opt in ("-B", "--basequality"):
                        base_quality = arg
		elif opt in ("-m", "--mappingquality"):
                        mapping_quality = arg
		elif opt in ("-f", "--normalalteredratio"):
			normal_alf_cutoff = arg
		elif opt in ("-c", "--normalalteredreads"):
                        normal_alt_cutoff = arg
		elif opt in ("-F", "--tumoralteredratio"):
			tumor_alf_cutoff = arg
		elif opt in ("-c", "--tumoralteredreads"):
			tumor_alt_cutoff = arg
		elif opt in ("-P", "--phasing_mode"):
			phasing_mode = arg

	return tumorfile, normalfile, reference, bed, outputpath, int(base_quality), int(mapping_quality), float(normal_alf_cutoff), int(normal_alt_cutoff), int(tumor_alt_cutoff), float(tumor_alf_cutoff), str(phasing_mode) 

if __name__ == '__main__':
        (tumorfile, normalfile, reference, bed, outputpath, base_quality, mapping_quality, normal_alf_cutoff, normal_alt_cutoff, tumor_alt_cutoff, tumor_alf_cutoff, phasing_mode) = main(sys.argv[1:])
    	#print tumor_alt_cutoff, tumor_alf_cutoff
	ref = pysam.FastaFile(reference)
	t_samfile = pysam.AlignmentFile(tumorfile)
	n_samfile = pysam.AlignmentFile(normalfile)
	tempfile = outputpath+"/raw_somatic_varants.txt"
	def process_reads(columns, columns_pos, refseq, filter_pos, phasing_info):
		alt_A=0; alt_C=0; alt_G=0; alt_T=0; refcount=0
		foundReadName=[]; foundReadName2=[]; altReadPosE_f_A=[]; altReadPosE_r_A=[]; altReadPosE_f_G=[]; altReadPosE_r_G=[]; altReadPosE_f_C=[]; altReadPosE_r_C=[]; altReadPosE_f_T=[]; altReadPosE_r_T=[]; altReadPosS_f_A=[]; altReadPosS_r_A=[]; altReadPosS_f_G=[]; altReadPosS_r_G=[]; altReadPosS_f_C=[]; altReadPosS_r_C=[]; altReadPosS_f_T=[]; altReadPosS_r_T=[]; hap_A=[]; hap_G=[]; hap_C=[]; hap_T=[]; hap_ref=[] 
		pileupcolumn = columns; refseq = ref_seq
		if pileupcolumn.pos == columns_pos: 
			for pileupread in pileupcolumn.pileups:
				if pileupread.alignment.is_proper_pair and not pileupread.alignment.is_duplicate:
					if not pileupread.is_del and not pileupread.is_refskip:
						if pileupread.alignment.mapping_quality >= mapping_quality and pileupread.alignment.query_qualities[pileupread.query_position] >= base_quality:
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

							if phasing_info == "Y":
								if pileupread.alignment.has_tag('HP'):
									if pileupread.alignment.get_tag('PC') >= 30:
										try:
											foundReadName2.index(pileupread.alignment.query_name)
										except ValueError:
											foundReadName2.append(pileupread.alignment.query_name)
											hap = pileupread.alignment.get_tag('HP')
											if pileupread.alignment.query_sequence[pileupread.query_position] == ref_seq.upper():
												hap_ref.append(hap)
											else:
												if altseq == "A":
													hap_A.append(hap)
												elif altseq == "G":
													hap_G.append(hap)
												elif altseq == "C":
													hap_C.append(hap)
												elif altseq == "T":
													hap_T.append(hap)
												

			if filter_pos == "Y":
				return alt_A, alt_C, alt_G, alt_T, refcount, altReadPosE_f_A, altReadPosE_r_A, altReadPosE_f_G, altReadPosE_r_G, altReadPosE_f_C, altReadPosE_r_C, altReadPosE_f_T, altReadPosE_r_T, altReadPosS_f_A, altReadPosS_r_A, altReadPosS_f_G, altReadPosS_r_G, altReadPosS_f_C, altReadPosS_r_C, altReadPosS_f_T, altReadPosS_r_T, hap_A, hap_G, hap_C, hap_T, hap_ref
			if filter_pos == "N":
				return alt_A, alt_C, alt_G, alt_T, refcount, hap_A, hap_G, hap_C, hap_T



	def filter_germline(chr, t_columns_pos, ref_seq, alt_base, phasing_info):
		n_alt_A=0; n_alt_C=0; n_alt_G=0; n_alt_T=0; n_refcount=0; n_hap_A=[]; n_hap_G=[]; n_hap_C=[]; n_hap_T=[];
		for rec in pysamstats.stat_variation(n_samfile, ref, chrom=str(chr), start=int(t_columns_pos), end=int(t_columns_pos)+1):
			chr = rec['chrom']; pos = rec['pos']; ref_seq = rec['ref']
			if pos == t_columns_pos:
				if rec['reads_pp'] > 4:
					if float(rec['mismatches_pp']/rec['reads_pp']) <= 0.10:
						for n_columns in n_samfile.pileup(chr, int(t_columns_pos), int(t_columns_pos)+1, truncate=True):
							if n_columns.pos == t_columns_pos:
								(n_alt_A, n_alt_C, n_alt_G, n_alt_T, n_refcount, n_hap_A, n_hap_G, n_hap_C, n_hap_T) = process_reads(n_columns, n_columns.pos, ref_seq, "N", phasing_mode);
					else:
						n_alt_A = rec['mismatches_pp']; n_alt_C = rec['mismatches_pp']; n_alt_G = rec['mismatches_pp']; n_alt_T = rec['mismatches_pp']; n_refcount = rec['reads_pp']-rec['mismatches_pp']
				else:
					n_alt_A=0; n_alt_G=0; n_alt_C=0; n_alt_T=0; n_refcount=rec['reads_pp']
 


		if alt_base == "A":
			return n_alt_A, n_refcount, n_hap_A
		if alt_base == "G":
			return n_alt_G, n_refcount, n_hap_G
		if alt_base == "T":
			return n_alt_T, n_refcount, n_hap_T
		if alt_base == "C":
			return n_alt_C, n_refcount, n_hap_C




	foundReg=[];	raws=[]
	cov = 0
	ftemp = open(tempfile, 'w')	
	print >>ftemp, "#chr\tpos\tref\talt\tt_ref\tt_alt\tn_ref\tn_alt\tjudge"
	print "calling raw variants..."
	for bedline in ( raw.strip().split() for raw in open(bed)):
 	     	for rec in pysamstats.stat_variation(t_samfile, ref, chrom=str(bedline[0]), start=int(bedline[1]), end=int(bedline[2])):
			if rec['reads_pp'] > 4:
				if (rec['mismatches_pp'] > 2) and (rec['insertions'] <= 3) and (rec['deletions'] <= 3):
					chr = rec['chrom']; pos = rec['pos'];
					ref_seq = rec['ref']
					for t_columns in t_samfile.pileup(chr, int(pos), int(pos)+1, truncate=True):
						try:
							foundReg.index(str(chr)+"\t"+str(t_columns.pos))
                    				except ValueError:
                                			foundReg.append(str(chr)+"\t"+str(t_columns.pos))
							(t_alt_A, t_alt_C, t_alt_G, t_alt_T, t_refcount, altReadPosE_f_A, altReadPosE_r_A, altReadPosE_f_G, altReadPosE_r_G, altReadPosE_f_C, altReadPosE_r_C, altReadPosE_f_T, altReadPosE_r_T, altReadPosS_f_A, altReadPosS_r_A, altReadPosS_f_G, altReadPosS_r_G, altReadPosS_f_C, altReadPosS_r_C, altReadPosS_f_T, altReadPosS_r_T, t_hap_A, t_hap_G, t_hap_C, t_hap_T, hap_ref) = process_reads(t_columns, t_columns.pos, ref_seq, "Y", phasing_mode)
							totalcount = t_refcount + t_alt_A + t_alt_G + t_alt_C + t_alt_T
							if totalcount > 4:
								if (t_alt_A > tumor_alt_cutoff) and ((t_alt_A/totalcount) > tumor_alf_cutoff):
									if (all(5 >= i for i in altReadPosS_f_A) and all(5 >= i for i in altReadPosS_r_A)) or (all(5 >= i for i in altReadPosE_f_A) and all(5 >= i for i in altReadPosE_r_A)):
										raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"A"+"\t"+str(t_refcount)+"\t"+str(t_alt_A)+"\t"+"0"+"\t"+"0"+"\t"+"clustered_pos")
									else:
										(n_alt_A, n_refcount, n_hap_A) = filter_germline(chr, t_columns.pos, ref_seq, "A", phasing_mode)
										if (n_refcount + n_alt_A) < 3:
											raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"A"+"\t"+str(t_refcount)+"\t"+str(t_alt_A)+"\t"+str(n_refcount)+"\t"+str(n_alt_A)+"\t"+"germline_uncovered")
										else: 
											if (n_alt_A < normal_alt_cutoff) or ((n_alt_A/(n_refcount + n_alt_A)) < normal_alf_cutoff):
											#	if not all(x==n_hap_A[0] for x in n_hap_A):
													if all(x==1 for x in t_hap_A):
														if (t_alt_A < 5) and (len(t_hap_A) > 1):
															t_hap_ref = len(hap_ref)
															t_hap_alt_A = t_hap_A.count(1)
															print pos, t_hap_ref, t_hap_alt_A, t_refcount, t_alt_A 
															raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"A"+"\t"+str(t_hap_ref)+"\t"+str(t_hap_alt_A)+"\t"+str(n_refcount)+"\t"+str(n_alt_A)+"\t"+"pass_to_test")
														else:
															raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"A"+"\t"+str(t_refcount)+"\t"+str(t_alt_A)+"\t"+str(n_refcount)+"\t"+str(n_alt_A)+"\t"+"pass_to_test")
													elif all(x==2 for x in t_hap_A):
														if (t_alt_A < 5) and (len(t_hap_A) > 1):
															t_hap_ref = len(hap_ref)
															t_hap_alt_A = t_hap_A.count(2)
															raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"A"+"\t"+str(t_hap_ref)+"\t"+str(t_hap_alt_A)+"\t"+str(n_refcount)+"\t"+str(n_alt_A)+"\t"+"pass_to_test")
														else:
															raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"A"+"\t"+str(t_refcount)+"\t"+str(t_alt_A)+"\t"+str(n_refcount)+"\t"+str(n_alt_A)+"\t"+"pass_to_test")
													else:
														if len(t_hap_A) < 2:
															raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"A"+"\t"+str(t_refcount)+"\t"+str(t_alt_A)+"\t"+str(n_refcount)+"\t"+str(n_alt_A)+"\t"+"pass_to_test")
														else:
															#print t_hap_A, pos , len(t_hap_A)
															raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"A"+"\t"+str(t_refcount)+"\t"+str(t_alt_A)+"\t"+str(n_refcount)+"\t"+str(n_alt_A)+"\t"+"phasing_failed")
											#	else:
											#		if len(n_hap_A) < 2:
											#			raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"A"+"\t"+str(t_refcount)+"\t"+str(t_alt_A)+"\t"+str(n_refcount)+"\t"+str(n_alt_A)+"\t"+"pass_to_test")
											#		else:
											#			print "normal", n_hap_A
											#			raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"A"+"\t"+str(t_refcount)+"\t"+str(t_alt_A)+"\t"+str(n_refcount)+"\t"+str(n_alt_A)+"\t"+"possible_germline")
											else:
												raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"A"+"\t"+str(t_refcount)+"\t"+str(t_alt_A)+"\t"+str(n_refcount)+"\t"+str(n_alt_A)+"\t"+"possible_germline")

								elif (t_alt_T > tumor_alt_cutoff) and ((t_alt_T/totalcount) > tumor_alf_cutoff): 
									if (all(5 >= i for i in altReadPosS_f_T) and all(5 >= i for i in altReadPosS_r_T)) or (all(5 >= i for i in altReadPosE_f_T) and all(5 >= i for i in altReadPosE_r_T)):
										raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"T"+"\t"+str(t_refcount)+"\t"+str(t_alt_T)+"\t"+"0"+"\t"+"0"+"\t"+"clustered_pos")
									else:
										(n_alt_T, n_refcount, n_hap_T) = filter_germline(chr, t_columns.pos, ref_seq, "T", phasing_mode)	
										if (n_refcount + n_alt_T) < 3:
											raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"T"+"\t"+str(t_refcount)+"\t"+str(t_alt_T)+"\t"+str(n_refcount)+"\t"+str(n_alt_T)+"\t"+"germline_uncovered")
										else:
											if (n_alt_T < normal_alt_cutoff) or ((n_alt_T/(n_refcount + n_alt_T)) < normal_alf_cutoff):
											#	if not all(x==n_hap_T[0] for x in n_hap_T):
													if  all(x==1 for x in t_hap_T):
														if (t_alt_T < 5) and (len(t_hap_T) > 1):
															t_hap_ref = len(hap_ref)		
															t_hap_alt_T = t_hap_T.count(1)
															#if pos == 19012132:
															print pos, t_hap_T, t_hap_alt_T, t_alt_T, hap_ref, t_hap_ref, t_refcount, n_refcount
															raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"T"+"\t"+str(t_hap_ref)+"\t"+str(t_hap_alt_T)+"\t"+str(n_refcount)+"\t"+str(n_alt_T)+"\t"+"pass_to_test")	
														else:			
															raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"T"+"\t"+str(t_refcount)+"\t"+str(t_alt_T)+"\t"+str(n_refcount)+"\t"+str(n_alt_T)+"\t"+"pass_to_test")

													elif all(x==2 for x in t_hap_T):
                                                                                                               	if (t_alt_T < 5) and (len(t_hap_T) > 1):
                                                                                                                	t_hap_ref = len(hap_ref)
                                                                                                                        t_hap_alt_T = t_hap_T.count(2)
                                                                                                                        raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"T"+"\t"+str(t_hap_ref)+"\t"+str(t_hap_alt_T)+"\t"+str(n_refcount)+"\t"+str(n_alt_T)+"\t"+"pass_to_test")
														else:
															raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"T"+"\t"+str(t_refcount)+"\t"+str(t_alt_T)+"\t"+str(n_refcount)+"\t"+str(n_alt_T)+"\t"+"pass_to_test")

													else:
														if len(t_hap_T) < 2:
															raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"T"+"\t"+str(t_refcount)+"\t"+str(t_alt_T)+"\t"+str(n_refcount)+"\t"+str(n_alt_T)+"\t"+"pass_to_test")
														else:
															#print t_hap_T, pos , len(t_hap_T)
															raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"T"+"\t"+str(t_refcount)+"\t"+str(t_alt_T)+"\t"+str(n_refcount)+"\t"+str(n_alt_T)+"\t"+"phasing_failed")
											#	else:
											#		if len(n_hap_T) < 2:
											#			raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"T"+"\t"+str(t_refcount)+"\t"+str(t_alt_T)+"\t"+str(n_refcount)+"\t"+str(n_alt_T)+"\t"+"pass_to_test")
											#		else:
											#		raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"T"+"\t"+str(t_refcount)+"\t"+str(t_alt_T)+"\t"+str(n_refcount)+"\t"+str(n_alt_T)+"\t"+"possible_germline")
											else:
												raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"T"+"\t"+str(t_refcount)+"\t"+str(t_alt_T)+"\t"+str(n_refcount)+"\t"+str(n_alt_T)+"\t"+"possible_germline")


								elif (t_alt_G > tumor_alt_cutoff) and ((t_alt_G/totalcount) > tumor_alf_cutoff): 
									if (all(5 >= i for i in altReadPosS_f_G) and all(5 >= i for i in altReadPosS_r_G)) or (all(5 >= i for i in altReadPosE_f_G) and all(5 >= i for i in altReadPosE_r_G)):
										raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"G"+"\t"+str(t_refcount)+"\t"+str(t_alt_G)+"\t"+"0"+"\t"+"0"+"\t"+"clustered_pos")
									else:
										(n_alt_G, n_refcount, n_hap_G) = filter_germline(chr, t_columns.pos, ref_seq, "G", phasing_mode)
										if (n_refcount + n_alt_G) < 3:
											raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"G"+"\t"+str(t_refcount)+"\t"+str(t_alt_G)+"\t"+str(n_refcount)+"\t"+str(n_alt_G)+"\t"+"germline_uncovered")
										else:
											if (n_alt_G < normal_alt_cutoff) or ((n_alt_G/(n_refcount + n_alt_G)) < normal_alf_cutoff):
										#		if not all(x==n_hap_G[0] for x in n_hap_G): 
													if all(x==1 for x in t_hap_G):
														if (t_alt_G < 5) and (len(t_hap_G) > 1):
															t_hap_ref = len(hap_ref)
															t_hap_alt_G = t_hap_G.count(1)
															raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"G"+"\t"+str(t_hap_ref)+"\t"+str(t_hap_alt_G)+"\t"+str(n_refcount)+"\t"+str(n_alt_G)+"\t"+"pass_to_test")
														else:
															raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"G"+"\t"+str(t_refcount)+"\t"+str(t_alt_G)+"\t"+str(n_refcount)+"\t"+str(n_alt_G)+"\t"+"pass_to_test")

													elif all(x==2 for x in t_hap_G):
														if (t_alt_G < 5) and (len(t_hap_G) > 1):
															t_hap_ref = len(hap_ref)
															t_hap_alt_G = t_hap_G.count(2)
															raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"G"+"\t"+str(t_hap_ref)+"\t"+str(t_hap_alt_G)+"\t"+str(n_refcount)+"\t"+str(n_alt_G)+"\t"+"pass_to_test")
														else:
															raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"G"+"\t"+str(t_refcount)+"\t"+str(t_alt_G)+"\t"+str(n_refcount)+"\t"+str(n_alt_G)+"\t"+"pass_to_test")

													else:   
														if len(t_hap_G) < 2:
															raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"G"+"\t"+str(t_refcount)+"\t"+str(t_alt_G)+"\t"+str(n_refcount)+"\t"+str(n_alt_G)+"\t"+"pass_to_test")
														else:
															#print t_hap_G, pos , len(t_hap_G)
															raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"G"+"\t"+str(t_refcount)+"\t"+str(t_alt_G)+"\t"+str(n_refcount)+"\t"+str(n_alt_G)+"\t"+"phasing_failed")
										#		else:
										#			if len(n_hap_G) < 2:
										#				raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"G"+"\t"+str(t_refcount)+"\t"+str(t_alt_G)+"\t"+str(n_refcount)+"\t"+str(n_alt_G)+"\t"+"pass_to_test")
										#			else:
										#				raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"G"+"\t"+str(t_refcount)+"\t"+str(t_alt_G)+"\t"+str(n_refcount)+"\t"+str(n_alt_G)+"\t"+"possible_germline")
											else:
												raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"G"+"\t"+str(t_refcount)+"\t"+str(t_alt_G)+"\t"+str(n_refcount)+"\t"+str(n_alt_G)+"\t"+"possible_germline")

								elif (t_alt_C > tumor_alt_cutoff) and ((t_alt_C/totalcount) > tumor_alf_cutoff):
									if (all(5 >= i for i in altReadPosS_f_C) and all(5 >= i for i in altReadPosS_r_C)) or (all(5 >= i for i in altReadPosE_f_C) and all(5 >= i for i in altReadPosE_r_C)):
							 			raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"C"+"\t"+str(t_refcount)+"\t"+str(t_alt_C)+"\t"+"0"+"\t"+"0"+"\t"+"clustered_pos")
									else:
										(n_alt_C, n_refcount, n_hap_C) = filter_germline(chr, t_columns.pos, ref_seq, "C", phasing_mode)
										if (n_refcount + n_alt_C) < 3:
											raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"C"+"\t"+str(t_refcount)+"\t"+str(t_alt_C)+"\t"+str(n_refcount)+"\t"+str(n_alt_C)+"\t"+"germline_uncovered")
										else:
											if (n_alt_C < normal_alt_cutoff) or ((n_alt_C/(n_refcount + n_alt_C)) < normal_alf_cutoff):
										#		if not all(x==n_hap_C[0] for x in n_hap_C):
													if  all(x==1 for x in t_hap_C):
														if (t_alt_C < 5) and (len(t_hap_C) > 1):
															t_hap_ref = len(hap_ref)
															t_hap_alt_C = t_hap_C.count(1)
															raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"C"+"\t"+str(t_hap_ref)+"\t"+str(t_hap_alt_C)+"\t"+str(n_refcount)+"\t"+str(n_alt_C)+"\t"+"pass_to_test")
														else:
															raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"C"+"\t"+str(t_refcount)+"\t"+str(t_alt_C)+"\t"+str(n_refcount)+"\t"+str(n_alt_C)+"\t"+"pass_to_test")

													elif all(x==2 for x in t_hap_C):
														if (t_alt_C < 5) and (len(t_hap_C) > 1):
															 t_hap_ref = len(hap_ref)
															 t_hap_alt_C = t_hap_C.count(2)
															 raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"C"+"\t"+str(t_hap_ref)+"\t"+str(t_hap_alt_C)+"\t"+str(n_refcount)+"\t"+str(n_alt_C)+"\t"+"pass_to_test")

														else:
															raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"C"+"\t"+str(t_refcount)+"\t"+str(t_alt_C)+"\t"+str(n_refcount)+"\t"+str(n_alt_C)+"\t"+"pass_to_test")
													else:
														if len(t_hap_C) < 2:
															raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"C"+"\t"+str(t_refcount)+"\t"+str(t_alt_C)+"\t"+str(n_refcount)+"\t"+str(n_alt_C)+"\t"+"pass_to_test")
														else:
															#print t_hap_C, pos , len(t_hap_C)
															raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"C"+"\t"+str(t_refcount)+"\t"+str(t_alt_C)+"\t"+str(n_refcount)+"\t"+str(n_alt_C)+"\t"+"phasing_failed")

										#		else: 
										#			if len(n_hap_C) < 2:
										#				raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"C"+"\t"+str(t_refcount)+"\t"+str(t_alt_C)+"\t"+str(n_refcount)+"\t"+str(n_alt_C)+"\t"+"pass_to_test")	
												#	else:
												#		raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"G"+"\t"+str(t_refcount)+"\t"+str(t_alt_G)+"\t"+str(n_refcount)+"\t"+str(n_alt_G)+"\t"+"possible_germline")
											else:
												raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"C"+"\t"+str(t_refcount)+"\t"+str(t_alt_C)+"\t"+str(n_refcount)+"\t"+str(n_alt_C)+"\t"+"possible_germline")
												
								else:
									raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"N"+"\t"+str(t_refcount)+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t"+"low_mutant_reads")
							else:
								raw_calls = (chr+"\t"+str(pos+1)+"\t"+ref_seq.upper()+"\t"+"N"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t"+"low_total_reads")
							raws.append(raw_calls)
							print >>ftemp, raw_calls

	ftemp.close()
	print "Done calling raw variants!"

if __name__ == '__main__':
        main(sys.argv[1:])

