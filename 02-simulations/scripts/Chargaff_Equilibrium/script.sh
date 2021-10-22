#!/bins/

# obtain working directory
pwd=$(pwd)
ncpu=2

# make directories for saving files
mkdir ../../data/Chargaff_Equilibrium
mkdir ../../data/Raw/Michael_Lynch
mkdir ../../figures/Chargaff_Equilibrium
mkdir ../../figures/Chargaff_Equilibrium/Michael_Lynch

# run chargaff equilibrium scripts
echo "Running Chargaff equilibrium simulation..."
Rscript Simulation_Run.R $pwd $ncpu

# obtain chargaff equilibrium plots
for save_file in "png" "pdf"
do
    for minimal in "TRUE" "FALSE"
    do
        echo "Obtainig Chargaff equilibration plots in" $save_file "format..."
        Rscript ChargaffEquilibrationPlots.R $pwd $save_file $minimal "FALSE"
    done
done

# # obtain mutation rates from michael lynch paper
# echo "Obtaining mutation rates from michael lynch paper..."
# Rscript Obtain_Lynch_Rates.R $pwd
# 
# run equilibrium simulation with lynch mutation rates
echo "Running Chargaff equilibrium simulation with lynch mutation rates..."
Rscript Lynch_Simulation_Run.R $pwd $ncpu
echo "Done!"

# obtain chargaff equilibrium plots
for save_file in "png" "pdf"
do
    for minimal in "TRUE" "FALSE"
    do
        echo "Obtainig chargaff equilibration plots for Lynch rates in" $save_file "format..."
        Rscript ChargaffEquilibrationPlots.R $pwd $save_file $minimal "TRUE"
    done
done