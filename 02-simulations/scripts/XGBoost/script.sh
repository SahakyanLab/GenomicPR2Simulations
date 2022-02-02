#!/bins/

# obtain working directory
pwd="$(pwd)/"
ncpu=$1

# make directories for saving files
mkdir -p ../../{data,figures}/XGBoost
mkdir -p ../../../03-symbolic_regression/data/Train

# hyper-parameter optimisation of xgboost ML model
echo "Optimising xgboost hyper-parameters..."
Rscript xgboost_hyper-opt.R $pwd $ncpu

# hyper-parameter optimisation of xgboost ML model
echo "Obtaining average xgboost model feature importance..."
Rscript xgboost_average.R $pwd $ncpu