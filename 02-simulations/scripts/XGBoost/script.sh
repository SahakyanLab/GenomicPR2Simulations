#!/bins/

# obtain working directory
pwd="$(pwd)/"
ncpu=2

# make directories for saving files
mkdir ../../data/XGBoost
mkdir ../../figures/XGBoost

# hyper-parameter optimisation of xgboost ML model
echo "Optimising xgboost hyper-parameters..."
Rscript xgboost_hyper-opt.R $pwd $ncpu

# hyper-parameter optimisation of xgboost ML model
echo "Obtaining average xgboost model feature importance..."
Rscript xgboost_average.R $pwd $ncpu