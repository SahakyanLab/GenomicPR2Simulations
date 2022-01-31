#!/bins/
pwd=$(pwd)
ncpu=$1
save_file="pdf"
trek_scale="FALSE"

# create folders for storage
mkdir -p ../../{data,figures}/Chargaff_Equilibrium
mkdir -p ../../figures/Chargaff_Equilibrium/Theoretical_GC
mkdir -p ../../data/Raw/Michael_Lynch

# chargaff equilibrium scripts
echo "Running Chargaff equilibrium simulation..."
Rscript Simulation_Run.R $pwd $ncpu

# chargaff equilibrium plots
echo "Obtainig Chargaff equilibration plots in" $save_file "format..."
Rscript ChargaffEquilibrationPlots.R $pwd $save_file "FALSE"

# mutation rates from michael lynch paper
echo "Obtaining mutation rates from michael lynch paper..."
Rscript Obtain_Lynch_Rates.R $pwd $trek_scale
Rscript GCcontentcurve.R $pwd $save_file

# equilibrium simulation with lynch mutation rates
echo "Running Chargaff equilibrium simulation with lynch mutation rates..."
Rscript Lynch_Simulation_Run.R $pwd $ncpu

# chargaff equilibrium plots
echo "Obtainig chargaff equilibration plots for Lynch rates in" $save_file "format..."
Rscript ChargaffEquilibrationPlots.R $pwd $save_file "TRUE"