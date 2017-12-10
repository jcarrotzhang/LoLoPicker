# LoLoPicker


### Haplotype filter adapting to 10X genomics (HP) is added in the recent update.

* use --phasing_mode Y to activate HP filter. 
* detailed description will appear in a manuscript of 10X genomic sequencing.


## Installation
```
git clone https://github.com/jcarrotzhang/LoLoPicker
cd LoLoPicker
python setup.py install
```
#### dependencies pysam 0.8.4; pysamstats
```
pip install pysam==0.8.4
pip install pysamstats
```
#### Please note that LoLoPicker is no longer available on PyPI

## Step one: calling raw, somatic variants using matched tumor/normal
```
python LoLoPicker_somatic.py -t tumor.bam -n normal.bam -r reference.fa -b interval.bed -o outputpath 
```
##### options: #####
* --basequality: only_count_reads_with_base_quality_above_cutoff 
* --mappingquality: only_count_reads_with_mapping_quality_above_cutoff 
* --tumoralteredreads: keep_variants_where_number_of_altered_reads_in_tumor_more_than_cutoff
* --normalalteredreads: keep_variants_where_number_of_altered_reads_in_normal_less_than_cutoff

## Step two: inspecting your control cohort
```
python LoLoPicker_control.py -l samplelist.txt -r reference.fa -o outputpath
```
##### options: ##### 
* --basequality 
* --mappingquality
* -n: number of threads

###### please provide your control panel in samplelist.txt using the following tab-delimited format:
```
Bam_file_of_each_control      control_sampleID 
```

# Step three: performing core stats
```
python LoLoPicker_stats.py -o outputpath 
```
##### options: #####
* --genome: for_analyzing_WGS_data 
* --SNPcutoff: keep_variants_present_in_number_of_normal_samples_less_than_cutoff 
* --intervalsize: size_of_the_targeted_region_of_your_experiment

##### Note: 
* For analyzing whole-genome sequencing data, please split your job by genomic intervals (e.g. chromosomes) and merge all your control_stats.txt files before going to step three. 
* For analyzing data from targeted re-sequencing, --intervalsize option is required. If the size of the targeted region of your experiment is 1000 base-pair, you should use interval size as 1000.

