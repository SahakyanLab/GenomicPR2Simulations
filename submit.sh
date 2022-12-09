#!/bins/

# run all files in the genome composition folder
cd ./01-genome_composition/scripts/
bash script.sh
cd ../../

# run simulation scripts
cd ./02-simulations/scripts/
bash script.sh
cd ../../

# run machine learning model
cd ./03-machine_learning/scripts/
bash script.sh
cd ../../

# run analysis of symbolic regression results
cd ./04-symbolic_regression/scripts/
bash script.sh
cd ../../