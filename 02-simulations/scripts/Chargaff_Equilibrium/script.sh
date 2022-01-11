#!/bins/
pwd=$(pwd)
ncpu=2
save_file="pdf"
minimal="FALSE"

mkdir ../../data/Chargaff_Equilibrium
mkdir ../../data/Raw/Michael_Lynch
mkdir ../../figures/Chargaff_Equilibrium
mkdir ../../figures/Chargaff_Equilibrium/Michael_Lynch

# chargaff equilibrium scripts
echo "Running Chargaff equilibrium simulation..."
Rscript Simulation_Run.R $pwd $ncpu

# chargaff equilibrium plots
echo "Obtainig Chargaff equilibration plots in" $save_file "format..."
Rscript ChargaffEquilibrationPlots.R $pwd $save_file $minimal "FALSE"

# mutation rates from michael lynch paper
echo "Obtaining mutation rates from michael lynch paper..."
Rscript Obtain_Lynch_Rates.R $pwd

# equilibrium simulation with lynch mutation rates
echo "Running Chargaff equilibrium simulation with lynch mutation rates..."
Rscript Lynch_Simulation_Run.R $pwd $ncpu
echo "Done!"

# chargaff equilibrium plots
echo "Obtainig chargaff equilibration plots for Lynch rates in" $save_file "format..."
Rscript ChargaffEquilibrationPlots.R $pwd $save_file $minimal "TRUE"