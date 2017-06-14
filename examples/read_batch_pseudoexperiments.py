"""
Compute final pval from t statistics of pseudoexperiments run in batch.
Run like:
$ python read_batch_pseudoexperiments.py <batch_output_dir> <t_obs>
"""

import sys
import os
import math

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
    return "%s \pm %s" % (pval, err)


if __name__ == '__main__':
    print read_batch_pseudoexperiments(sys.argv[1], float(sys.argv[2]))
