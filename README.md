# LoLoPicker

#Installation
#######git clone https://github.com/jcarrotzhang/LoLoPicker 
#######cd LoLoPicker 
#######python setup.py install 


#Step one: calling raw, somatic variants using matched tumor/normal 

LoLoPicker_somatic.py -t tumorfile -n normalfile -r reference -b bedfile -o outputpath

#Step two: inspecting your control cohort

LoLoPicker_control.py -l samplelist -r reference -n thread -o outputpath

#Step three: performing core stats

LoLoPicker_stats.py -o outputpath

