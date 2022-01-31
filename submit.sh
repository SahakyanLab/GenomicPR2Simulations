#!/bins/
parallel="FALSE"
pca="FALSE"
ncpu=4

# run all files in the genome composition folder
cd ./01-genome_composition/scripts/
bash script.sh
cd ../../

# run simulation scripts
cd ./02-simulations/scripts/Main_Simulation/
bash script.sh $ncpu $parallel
cd ../../

# run chargaff equilibrium scripts
cd ./02-simulations/scripts/Chargaff_Equilibrium/
bash script.sh $ncpu
cd ../../

# run xgboost machine learning model
cd ./02-simulations/scripts/XGBoost/
bash script.sh $ncpu
cd ../../

# run principal component analysis (optional)
if [ "$pca" == "TRUE" ]
    then 
        cd ./02-simulations/scripts/PCA/
        bash script.sh $ncpu
        cd ../../
fi

# run analysis of symbolic regression results
cd ./03-symbolic_regression/scripts/
bash script.sh $ncpu
cd ../../