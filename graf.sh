
# generate distance from plink files
./graf-master/bin/graf -plink /data/tangw3/NCI_MD_Prostate_Cancer_Study/delivery/subjects -pop test.txt

# generate ancestries from distance files with self-identified
./graf-master/bin/PlotPopulations.pl test.txt test.results.txt -spf spf.txt
./graf-master/bin/PlotPopulations.pl test.txt test.spf.pops.png -spf spf.txt -pops 1,3,2,4

#generate relateness files
./graf-master/bin/graf -plink /data/tangw3/NCI_MD_Prostate_Cancer_Study/delivery/subjects
#generate relateness plots
./graf-master/bin/PlotGraf.pl graf_rel_20201103_1657.txt hgmr.60.20.png 1 -xmax 60 -ymax 20
