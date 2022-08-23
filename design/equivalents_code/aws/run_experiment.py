import boto3
import collections
from datetime import datetime
import itertools
import numpy as np
import os
import pandas as pd
import pickle
import scipy.linalg
import scipy.optimize
import scipy.stats
import sys

BUCKET_NAME = 'equivalents'

def _solve_lambda_n(Sigma, n, tol=1e-10):
    if np.isclose(n, np.linalg.matrix_rank(Sigma)):
        return 0
    else:
        _, S, _ = scipy.linalg.svd(Sigma, full_matrices=False)
        exp_subset_size_func = lambda lam: np.sum(S / (S + lam))

        opt_result = scipy.optimize.minimize_scalar(
            lambda lam: (exp_subset_size_func(lam) - n)**2,
            method='bounded',
            bounds=(0, d * max(S) / n),
            tol=tol,
        )
        return opt_result.x

def run_experiment(param):
    return {
        'param': param,
    }

def get_current_results():
    s3 = boto3.resource('s3')
    data = []
    for obj in s3.Bucket(BUCKET_NAME).objects.all():
        result = pickle.loads(obj.get()['Body'].read())
        data.append(result)
    return pd.DataFrame(data)

if __name__ == '__main__':
    start = datetime.now()
    rng = np.random.RandomState(42)

    param_grid = [0, 1, 2, 3]

    if os.environ.get('AWS_BATCH_JOB_ARRAY_INDEX') is None:
        # local, exit early printing size, see `build_and_push.sh`
        print(len(param_grid))
        sys.exit(0)
    print(os.environ.get('AWS_BATCH_JOB_ARRAY_INDEX'))

    param = param_grid[int(os.environ.get('AWS_BATCH_JOB_ARRAY_INDEX'))]

    current_results = get_current_results()
    try:
        previous_matching_runs = (
            (current_results['param'] == param)
        ).sum()
        if previous_matching_runs != 0:
            print("Previous run found, exiting")
    except:
        pass

    print("No previous runs, executing")

    result = run_experiment(param)
    print(result)

    result_filename = f"result-param={param}.pkl"
    with open(result_filename, "wb") as f:
        pickle.dump(result, f)

    s3 = boto3.client('s3')
    s3.upload_file(result_filename, BUCKET_NAME, result_filename)
