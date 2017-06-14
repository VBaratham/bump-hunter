import os
import sys
import re
import math
import time
import json
import argparse
import subprocess
import tempfile

from ROOT import TFile

from read_batch_pseudoexperiments import read_batch_pseudoexperiments

"""
Run many pseudoexperiments, compile the results
"""

def create_script(scriptname, args):
    with open(scriptname, "w") as f:
        print >>f, "#!/bin/bash"
        print >>f, "#$ -V"
        print >>f, ""
        print >>f, '\n'.join('# %s' % line for line in json.dumps(args.__dict__, sort_keys=True,
                                                                  indent=4).split('\n'))
        print >>f, '# %s' % os.path.split(os.getcwd())[-1]
        print >>f, ""
        arg_str = '--rootfile %s --name %s --num %s --formula "%s"' % (args.rootfile, args.name,
                                                                       args.num_per_host, args.formula)
        print >>f, "python ${BUMPHUNTER_LIB}/examples/pseudoexperiments.py %s" % arg_str


def submit_jobs(scriptname, args):
    """
    Submit batch jobs described by args in the current directory
    """
    create_script(scriptname, args)
    qsub_out = subprocess.check_output("qsub -t 1-%s:1 %s" % (args.num_hosts, scriptname), shell=True)
    try:
        return re.search("Your job-array (\d+)", qsub_out).group(1)
    except AttributeError:
        return None


def wait_jobs(jobnum):
    while True:
        num_left = subprocess.check_output('qstat -u vbaratha | grep "%s" | wc -l' % jobnum, shell=True)
        if num_left == 0:
            break
        time.sleep(60)

    
def main(args):
    timestamp = int(time.time())
    
    dirname = "pseudoexperiments_%s" % timestamp
    os.mkdir(dirname)
    os.chdir(dirname)

    scriptname = "pseudoexperiments.sh"
    jobnum = submit_jobs(scriptname, args)
    if not jobnum:
        print "Could not read job number from qsub command"
        sys.exit(1)
    if args.t_obs: # User wants to calculate and display final p-val
        wait_jobs(jobnum)
        print read_batch_pseudoexperiments(os.getcwd(), args.t_obs)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run pseudoexperiments for a given histogram in batch")
    parser.add_argument('--rootfile', type=os.path.abspath, required=True,
                        help='name of root file containing data histogram')
    parser.add_argument('--t-obs', type=float,
                        help='observed t statistic. Only include if you want this process to '
                        'calculate final pval (can always calculate it later using '
                        'read_batch_pseudoexperiments.py')
    parser.add_argument('--name', type=str, default='signal',
                        help='name of histogram in root file')
    parser.add_argument('--num-hosts', type=int, default=100,
                        help='number of hosts to run on')
    parser.add_argument('--num-per-host', type=int, default=100,
                        help='number of pseudoexperiments to run on each host')
    parser.add_argument('--formula', type=str, default="[0]*exp(-[1] - [2]*x - [3]*y)",
                        help='formula to fit')

    args = parser.parse_args()
    main(args)
