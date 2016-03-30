# LoLoPicker

#Installation
#####git clone https://github.com/jcarrotzhang/LoLoPicker
#####cd LoLoPicker
#####python setup.py install
######dependencies pysam (>=0.8.4) pysamstats
######Please note that LoLoPicker is no longer available on PyPI

#Step one: calling raw, somatic variants using matched tumor/normal

python LoLoPicker_somatic.py -t tumor.bam -n normal.bam -r reference.fa -b interval.bed (e.g. CCDS_in_bed_format) -o outputpath 
######(options: --basequality only_count_reads_with_base_quality_above_cutoff --mappingquality only_count_reads_with_mapping_quality_above_cutoff --normalalteredreads keep_variants_where_number_of_altered_reads_in_normal_less_than_cutoff)

#Step two: inspecting your control cohort

python LoLoPicker_control.py -l samplelist.txt -r reference.fa -o outputpath
######(options: --basequality --mappingquality -n thread)

#####please provide your control panel in samplelist.txt using the following tab-delimited format:
######Bam_file_of_each_control      \t      control_sampleID

#Step three: performing core stats

python LoLoPicker_stats.py -o outputpath 

######(options: --basecov --genome)

#####For analyzing whole-genome sequencing data, please split your job by genomic intervals (e.g. chromosomes) and merge all your control_stats.txt files before going to step three.

