import os
import sys
import argparse
from datetime import datetime

from ROOT import TFile, TF2

from bumphunter.bumphunter import BumpHunter

"""
Run one batch of pseudoexperiments, spitting t statistics to stdout
"""

def get_fit_fcn(histo, args):
    fit_fcn = TF2(
        "expo2", args.formula,
        histo.GetXaxis().GetXmin(),
        histo.GetXaxis().GetXmax(),
        histo.GetYaxis().GetXmin(),
        histo.GetYaxis().GetXmax(),
    )
    fit_fcn.SetNpx(histo.GetNbinsX())
    fit_fcn.SetNpy(histo.GetNbinsY())

    fit_fcn.SetParameter(0, 1000)
    fit_fcn.SetParameter(1, 1)
    fit_fcn.SetParameter(2, 0.2)
    fit_fcn.SetParameter(3, 0.2)

    return fit_fcn


def main(args):
    f = TFile.Open(args.rootfile)
    bh = BumpHunter.from_rootfile(f)
    fit_fcn = get_fit_fcn(bh.histo, args)
    bh.pseudoexperiments(args.num, fit_fcn)
    f.Close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run pseudoexperiments for a given histogram")
    parser.add_argument('--rootfile', type=os.path.abspath, required=True,
                        help='name of root file containing data histogram')
    parser.add_argument('--num', type=int, default=100,
                        help='number of pseudoexperiments to run')
    parser.add_argument('--formula', type=str, default="[0]*exp(-[1] - [2]*x - [3]*y)",
                        help='formula to fit')
    # TODO: allow passing in initial parameter values

    args = parser.parse_args()
    main(args)
