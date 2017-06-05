"""
Creates a 2d histogram consisting of a Gaussian signal and exponentially decaying background,
runs 2d bumphunter, and shows results in ROOT windows
"""

import sys
import argparse
import math
from datetime import datetime

from ROOT import TH2F, TF2, TRandom3, TCanvas, TFile

from bumphunter.bumphunter import BumpHunter2D
from bumphunter.utils import MultiOutstream


def make_histo(title, args, include_signal=True):
    """
    Make and return a histogram describe by the above parameters
    """
    histo = TH2F(title, "2D BumpHunter Test %s" % title, args.nbins, 0,
                 args.nbins*args.binsize, args.nbins, 0, args.nbins*args.binsize)
    histo.GetXaxis().SetTitle("X")
    histo.GetXaxis().SetTitleOffset(2)
    histo.GetYaxis().SetTitle("Y")
    histo.GetYaxis().SetTitleOffset(2)
    
    rnd = TRandom3()
    rnd.SetSeed(0)
    for i in range(args.nbkg):
        histo.Fill(rnd.Exp(args.bkg_mean), rnd.Exp(args.bkg_mean))
    if include_signal:
        for i in range(args.nsig):
            x, y = rnd.Gaus(args.sig_x, args.sig_spread_x), rnd.Gaus(args.sig_y, args.sig_spread_y)
            # print x, y
            histo.Fill(x, y)

    return histo


def make_histo_1d(title, args, include_signal=True):
    """
    Make and return a 1D TH2F. All the data will have a y coordinate of args.binsize/2.0
    (to keep the data in the first bin)
    but the fit needs at least 4 points in y, so we make it 4 bins wide along y
    """
    histo = TH2F(title, "1D BumpHunter2D Test %s" % title, args.nbins, 0,
                 args.nbins*args.binsize, 4, 0, 4*args.binsize)
    histo.GetXaxis().SetTitle("X")
    histo.GetXaxis().SetTitleOffset(2)
    histo.GetYaxis().SetTitle("Y")
    histo.GetYaxis().SetTitleOffset(2)
    
    rnd = TRandom3()
    rnd.SetSeed(0)
    for i in range(args.nbkg):
        histo.Fill(rnd.Exp(args.bkg_mean), args.binsize/2.0)
    if include_signal:
        for i in range(args.nsig):
            x, y = rnd.Gaus(args.sig_x, args.sig_spread_x), args.binsize/2.0
            # print x, y
            histo.Fill(x, y)

    return histo


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
    timestamp = datetime.now().isoformat()

    histo = make_histo('signal', args, include_signal=True)
    bh = BumpHunter2D(histo, fit_fcn=get_fit_fcn(histo))

    f = TFile("runs/%s.root" % timestamp, "NEW")

    c2 = TCanvas("c2")
    bh.histo.Draw("LEGO2 HIST")
    bh.histo.Write()
    # c3 = TCanvas("c3")
    # bh.bkg_histo.Draw("LEGO2 HIST")
    bh.bkg_histo.Write()

    f.Close()

    with open("runs/%s.stdout" % timestamp, "w") as outfile:
        out = MultiOutstream(sys.stdout, outfile)

        import json
        print >>out, json.dumps(args.__dict__, indent=4, sort_keys=True)

        print >>out, "\nRunning BumpHunter on test data..."

        t, best_p, best_leftedge, best_width = bh.get_best_bump()

        print >>out, "Done"
        print >>out, "----------------------------------"
        # print "P-value: p = %s" % best_p # this is not very meaningful
        print >>out, "test statistic: t = %s" % t
        print >>out, "Left edge: (%s, %s):" % ( (best_leftedge[0] - 1) * args.binsize,
                                                (best_leftedge[1] - 1) * args.binsize )
        print >>out, "Width: (%s, %s):" % (best_width[0]*args.binsize, best_width[1]*args.binsize)
        print >>out, ""
        print >>out, "Running %s pseudoexperiments..." % args.num_pseudo
        print >>out, ""

        final_pval, err = bh.pseudoexperiments(
            args.num_pseudo, target_pval=args.pval, progress_out=out)

        print >>out, ""
        print >>out, "Final p-value: p = %s \pm %s" % (final_pval, err)

    # Hang the program so the user can look at the output
    raw_input("Press enter to quit")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run BumpHunter on randomly generated 2D data")
    parser.add_argument('--nbkg', type=int, default=100000,
                        help='Number of background events to generate')
    parser.add_argument('--nsig', type=int, default=200,
                        help='Number of signal events to generate')
    parser.add_argument('--nbins', type=int, default=40,
                        help='Number of bins along each dimension')
    parser.add_argument('--binsize', type=float, default=0.5,
                        help='Size of each bin along each dimension')
    parser.add_argument('--bkg-mean', type=float, default=8,
                        help='Mean of background exponential distribution')
    parser.add_argument('--sig-x', type=float, default=5,
                        help='x position of signal peak')
    parser.add_argument('--sig-y', type=float, default=12,
                        help='y position of signal peak')
    parser.add_argument('--sig-spread-x', type=float, default=.6,
                        help='stdev along x axis of signal peak')
    parser.add_argument('--sig-spread-y', type=float, default=.2,
                        help='stdev along y axis of signal peak')
    parser.add_argument('--num-pseudo', type=int, default=10,
                        help='max number of pseudoexperiments to run to estimate final p-val')
    parser.add_argument('--pval', type=float, required=False,
                        help='run pseudoexperiments until this pval is reached, not exceeding '
                        '--num-pseudo')

    args = parser.parse_args()
    main(args)
