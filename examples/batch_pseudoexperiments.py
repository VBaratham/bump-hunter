import os
import sys
import math
import time
import argparse
import subprocess
import tempfile

from ROOT import TFile

"""
Run many pseudoexperiments, compile the results
"""

def create_script(scriptname, args):
    with open(scriptname, "w") as f:
        print >>f, "#!/bin/bash"
        print >>f, "#$ -V"
        arg_str = '%s --name %s --num %s --formula "%s"' % (args.rootfile, args.name,
                                                            args.num_per_host, args.formula)
        print >>f, "python ${BUMPHUNTER_LIB}/examples/pseudoexperiments.py %s" % arg_str


def submit_jobs(scriptname, args):
    """
    Submit batch jobs described by args in the current directory
    """
    create_script(scriptname, args)

    for i in range(args.num_hosts):
        os.system("qsub %s" % scriptname)

            
def read_jobs(scriptname, args):
    while True:
        # num_left = int(subprocess.check_output('qstat -u vbaratha | grep "%s" | wc -l' % scriptname))
        # if num_left != 0:
        #     break
        if len(os.listdir(os.getcwd())) == 2 * args.num_hosts + 1: # 2 output files/host + 1 bash script
            break
        time.sleep(60)

    num_greater, num_tot = 0, 0
    for fn in os.listdir(os.getcwd()):
        if fn.split('.')[-1].startswith('o'):
            with open(fn) as f:
                all_t = [float(t) for t in f.readlines()]
            num_greater += len([t for t in all_t if t > args.t_obs])
            num_tot += len(all_t)

    pval = float(num_greater)/float(num_tot),
    err = math.sqrt(pval * (1.0 - pval)/num_tot)
    return pval, err


def main(args):
    timestamp = int(time.time())
    cwd = os.getcwd()
    
    # dirname = tempfile.mkdtemp()
    dirname = "pseudoexperiments_%s" % timestamp
    os.mkdir(dirname)
    os.chdir(dirname)

    scriptname = "pseudoexperiments.sh"
    submit_jobs(scriptname, args)
    print read_jobs(scriptname, args)
    
    os.chdir(cwd)
    # os.removedirs(dirname)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run pseudoexperiments for a given histogram")
    parser.add_argument('rootfile', type=os.path.abspath, 
                        help='name of root file containing data histogram')
    parser.add_argument('--t-obs', type=float, required=True,
                        help='observed t statistic')
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
