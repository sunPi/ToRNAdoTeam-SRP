import sys
import os

arguments = []

for s in range(len(sys.argv)):
    arguments.append(sys.argv[s])
    #print(arguments[s])

def process(string):
    p = os.popen(string)
    out = p.read()
    print(out)

process('module load sratoolkit')
process('module load fastqc')
process('module load R/4.0.0')


process('prefetch -v SRR12621861 SRR12621862 SRR12621863 SRR12621864 SRR12621865 SRR12621866 SRR12621867 SRR12621868 SRR12621869 SRR12621870 SRR12621871 SRR12621872')
process('for i in *.sra; do; fasterq-dump -vpx $i; done')
process('fasterq-dump -vpx SRR12621861 SRR12621862 SRR12621863 SRR12621864 SRR12621865 SRR12621866 SRR12621867 SRR12621868 SRR12621869 SRR12621870 SRR12621871 SRR12621872')
process('for i in .fastqc; do; fastqc $i; done')
process('fastqc SRR12621861 SRR12621862 SRR12621863 SRR12621864 SRR12621865 SRR12621866 SRR12621867 SRR12621868 SRR12621869 SRR12621870 SRR12621871 SRR12621872')


#Use Linux shell to install neccessary LINUX packages for R packages
process('sudo apt install libgfortran4')
process('sudo apt install libudunits2-dev')
process('sudo apt install libgdal-dev')

#Install BioConductor
#https://cole-trapnell-lab.github.io/monocle3/docs/installation/

#Install Seurat
#In R: > install.packages('Seurat')

#Download UMI cell count matrix
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157783/suppl/GSE157783_IPDCO_hg_midbrain_UMI.tar.gz

