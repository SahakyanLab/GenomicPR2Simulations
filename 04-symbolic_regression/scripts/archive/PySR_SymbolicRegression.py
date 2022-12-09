import numpy as np
from pysr import pysr, best
import pandas as pd
from sympy import *
from sklearn import preprocessing
from sklearn.model_selection import train_test_split

# import data set
df = pd.read_csv("../../../data/simulation/PR2_compliance/Non_Symmetric-uniform_species-eukaryotes.csv")
df = df.iloc[:,1:]

# split data
for ind in range(df.shape[1]):

    var = df.columns[ind]
    X = df.loc[:, df.columns != var]
    y = df.loc[:, df.columns == var].to_numpy()
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.3, random_state = 42)

    # normalise columns
    min_max_scaler = preprocessing.MinMaxScaler()

    # train set
    np_scaled = min_max_scaler.fit_transform(X_train)
    df_normalised_X_train = pd.DataFrame(np_scaled, columns = X_train.columns)

    # test set
    np_scaled = min_max_scaler.fit_transform(X_test)
    df_normalised_X_test = pd.DataFrame(np_scaled, columns = X_test.columns)

    # Learn equations
    equations = pysr(
        df_normalised_X_train, y_train, 
        niterations = 4, 
        procs = 8,
        denoise = True,
        fast_cycle = True,
        ncyclesperiteration = 1000,
        optimizer_algorithm = 'NelderMead',
        maxsize = 40, 
        equation_file = f'../../data/symbolic_regression/{var}.csv', 
        binary_operators = [
        '+', '-', '*'
    ])
    equations.to_csv(f'../../data/symbolic_regression/{var}.csv',index = False)

    best_expression = equations.iloc[equations.MSE.argmin()]
    ex1 = best_expression.sympy_format
    ex2 = ex1
    for a in preorder_traversal(ex1):
        if isinstance(a, Float):
            ex2 = ex2.subs(a, round(a, 1))

    results = pd.DataFrame({
        'y_true': y_true.reshape(-1), 
        'y_pred': y_pred
    })

    # save results
    results.to_csv(f'../../data/symbolic_regression/{var}_prediction.csv', index = False)

    # save hyperparameters 
    with open("SR_hyperparams.txt", "a") as f:
        if var == 'kag':
            f.write(
                f'{var}: \n'
                f'Full =  {ex1}\n'
                f'Rounded = {ex2}\n'
                'niterations         = 3\n'
                'procs               = 8\n'
                'denoise             = True\n'
                'fast_cycle          = True\n'
                'ncyclesperiteration = 1000\n'
                'optimizer_algorithm = NelderMead\n'
                'maxsize             = 40\n'
                'binary_operators    = [+, -, *]\n'
                'unary_operators     = None'
            )
        else:
            f.write(
                f'\n\n{var}:\n'
                f'Full =  {ex1}\n'
                f'Rounded = {ex2}\n'
                'niterations         = 3\n'
                'procs               = 8\n'
                'denoise             = True\n'
                'fast_cycle          = True\n'
                'ncyclesperiteration = 1000\n'
                'optimizer_algorithm = NelderMead\n'
                'maxsize             = 40\n'
                'binary_operators    = [+, -, *]\n'
                'unary_operators     = None'
            )
