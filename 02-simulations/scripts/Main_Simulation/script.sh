#!/bins/

# obtain working directory
pwd=$(pwd)
ncpu=$1
parallel=$2
save_file="png" # png or pdf
runs=25000000

# make directories for saving files
mkdir -p ../../data/Main_Simulation/{Non_Symmetric-uniform,Non_Symmetric-normal,Strand_Symmetric-normal}
mkdir -p ../../figures/Main_Simulation/{Non_Symmetric-uniform,Non_Symmetric-normal,Strand_Symmetric-normal}/Scaling_{1,2,5,10}

if [ "$parallel" == "TRUE" ]
    then
        echo "Running simulation in parallel execution..."
        # run non-symmetric simulations
        for dist in "uniform" "normal"
        do
            echo "Running all Non_Symmetric" "$dist" "simulations..."
            Rscript Simulation_Run.R $pwd $ncpu $runs "Non_Symmetric" $dist 1 "FALSE" &
            Rscript Simulation_Run.R $pwd $ncpu $runs "Non_Symmetric" $dist 2 "FALSE" &
            Rscript Simulation_Run.R $pwd $ncpu $runs "Non_Symmetric" $dist 5 "FALSE" &
            Rscript Simulation_Run.R $pwd $ncpu $runs "Non_Symmetric" $dist 10 "FALSE" &
        done

        # run strand-symmetric simulations
        echo "Running all Strand_Symmetric" "$dist" "simulations..."
        Rscript Simulation_Run.R $pwd $ncpu $runs "Strand_Symmetric" normal 1 "TRUE" &
        Rscript Simulation_Run.R $pwd $ncpu $runs "Strand_Symmetric" normal 2 "TRUE" &
        Rscript Simulation_Run.R $pwd $ncpu $runs "Strand_Symmetric" normal 5 "TRUE" &
        Rscript Simulation_Run.R $pwd $ncpu $runs "Strand_Symmetric" normal 10 "TRUE" &
    else
        echo "Running simulation in sequential execution..."
        # run non-symmetric simulations
        for dist in "uniform" "normal"
        do
            for scaling in 1 2 5 10
            do
                echo "Running Non_Symmetric" "$dist" "simulation with scaling factor" "$scaling" "..."
                Rscript Simulation_Run.R $pwd $ncpu $runs "Non_Symmetric" $dist $scaling "FALSE"
            done
        done

        # run strand-symmetric simulations
        for scaling in 1 2 5 10
        do
            echo "Running Strand_Symmetric normal simulation with scaling factor" "$scaling" "..."
            Rscript Simulation_Run.R $pwd $ncpu $runs "Strand_Symmetric" normal $scaling "TRUE"
        done
fi

# plots for the main simulations
for scaling in 1 2 5 10
do
   for Tolerance in "TRUE" "FALSE"
   do
      for dist in "uniform" "normal"
         echo "Obtaining plots for Non-symmetric $dist results with Tolerance = $Tolerance..."
         Rscript Non_Symmetric.R $pwd $save_file $scaling $dist $Tolerance
   done

   echo "Obtaining plots for Strand-symmetric normal results..."
   Rscript Strand_Symmetric-normal.R $pwd $save_file $scaling $Tolerance

   echo "Creating animation of strand symmetry skews evolution..."
   Rscript Strand_Symmetric-normal-evolution.R $pwd $scaling
done