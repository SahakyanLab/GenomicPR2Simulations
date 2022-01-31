#!/bins/

pwd=$(pwd)
ncpu=$1
runs=10000000
trek_scale="TRUE"
save_file="pdf"

# folders for storage
mkdir -p ../data/Simulation/Test
mkdir -p ../figures/{Test,Train}

# simulation to obtain test data for Eureqa evaluation
Rscript Simulation_Run.R $pwd $ncpu $runs

# filter results of test runs
Rscript FilterSimulationBatches.R $pwd

# test data through symbolic expressions from Eureqa
for test in "TRUE" "FALSE"
do
    Rscript EureqaEquations.R $pwd $test
done

# plots comparing test vs. training results
for test in "TRUE" "FALSE"
do
    Rscript SymbolicRegressionPlots.R $pwd $save_file $test
done

# rate constant pred vs. true of the 17 species
Rscript ExperimentalRateConst.R $pwd $trek_scale