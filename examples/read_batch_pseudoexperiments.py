"""
Compute final pval from t statistics of pseudoexperiments run in batch.
"""

import sys
import os
import time
import math
import argparse


def read_batch_pseudoexperiments(dirname, t_obs):
    num_greater, num_tot = 0, 0
    for fn in os.listdir(dirname):
        fn = os.path.join(dirname, fn)
        if fn.split('.')[-2].startswith('o') or fn.split('.')[-1].startswith('o'):
            with open(fn) as f:
                all_t = [float(t) for t in f.readlines()]
            num_greater += len([t for t in all_t if t > t_obs])
            num_tot += len(all_t)

    pval = float(num_greater)/float(num_tot)
    err = math.sqrt(pval * (1.0 - pval)/num_tot)

    return pval, err, num_tot


def main(args):
    timestamp = int(time.time())
    pval, err, n = read_batch_pseudoexperiments(args.directory, args.t_obs)
    readable = "with t_obs = %s, p-val = %s \pm %s from %s pseudoexperiments" % (args.t_obs, pval, err, n)
    print readable

    if not args.nowrite:
        with open(os.path.join(args.directory, "final_pval_%s.txt" % timestamp), "w") as f:
            print >>f, readable


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Read pseudoexperiments t statistics and "
                                     "compute final p-value")
    parser.add_argument('--dir', type=os.path.abspath, required=True, dest='directory',
                        help='directory containing batch output')
    parser.add_argument('--t-obs', type=float, required=True,
                        help='observed t statistic')
    parser.add_argument('--nowrite', action="store_true", required=False,
                       help='Do not write an output file with the p-val')
    args = parser.parse_args()
    main(args)
