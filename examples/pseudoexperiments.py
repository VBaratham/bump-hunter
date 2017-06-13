import sys
import argparse
from datetime import datetime

from ROOT import TFile, TF2

from bumphunter.bumphunter import BumpHunter2D

"""
Run one batch of pseudoexperiments, spitting t statistics to stdout
"""

def get_fit_fcn(histo):
    fit_fcn = TF2(
        "expo2", "[0]*exp(-[1] - [2]*x - [3]*y)",
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
    hist = f.Get(args.name)
    fit_fcn = get_fit_fcn(hist)
    bh = BumpHunter2D(hist, fit_fcn) # TODO: need same config that was used to calculate t_obs
    bh.pseudoexperiments_simple(args.num, fit_fcn)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run pseudoexperiments for a given histogram")
    parser.add_argument('rootfile', type=str,
                        help='name of root file containing data histogram')
    parser.add_argument('--name', type=str, default='signal',
                        help='name of histogram in root file')
    parser.add_argument('--num', type=int, default=100,
                        help='number of pseudoexperiments to run')
    parser.add_argument('--formula', type=str, default="[0]*exp(-[1] - [2]*x - [3]*y)",
                        help='formula to fit')

    args = parser.parse_args()
    main(args)
