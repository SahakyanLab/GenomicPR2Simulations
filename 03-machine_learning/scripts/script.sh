#!/bins/

# obtain working directory
pwd="$(pwd)/"
seed=2022
cv=6
cvrep=1
ml_seed=1234
avg_ml_seed=123

# PCA if you want to...
for species in "prokaryotes" "eukaryotes" "viruses"
do
    Rscript PCA.R $pwd $species
done

# run xgboost model
species="eukaryotes"
Rscript Process.R "$pwd/" $species $seed $ml_seed $avg_ml_seed $cv $cvrep