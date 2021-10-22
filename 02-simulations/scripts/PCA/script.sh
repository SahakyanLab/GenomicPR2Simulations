#!/bins/

# obtain working directory
pwd=$(pwd)
ncpu=2

# make directories for saving files
mkdir ../../figures/PCA

for species in "eukaryotes" "prokaryotes" "viruses"
do
    Rscript PCA.R $pwd $ncpu $species
done