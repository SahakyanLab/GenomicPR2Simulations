#!/bins/

# obtain working directory
pwd="$(pwd)/"
runs=25000000
ncpu=4
save_as="pdf"

# run simulations and generate plots
Rscript Process.R $pwd $runs $ncpu $save_as