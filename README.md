# LoLoPicker

#Installation
######git clone https://github.com/jcarrotzhang/LoLoPicker 
######cd LoLoPicker 
######python setup.py install 
######dependencies pysam (>=0.8.4) pysamstats



#Step one: calling raw, somatic variants using matched tumor/normal 

python LoLoPicker_somatic.py -t tumorfile -n normalfile -r reference -b bedfile -o outputpath

#Step two: inspecting your control cohort

python LoLoPicker_control.py -l samplelist -r reference -n thread -o outputpath

#Step three: performing core stats

python LoLoPicker_stats.py -o outputpath

