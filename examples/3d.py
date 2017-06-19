"""
Creates a 3d histogram consisting of a Gaussian signal and exponentially decaying background,
runs 3d bumphunter, and shows results in ROOT windows
"""

import sys
import argparse
import math
from datetime import datetime

from ROOT import TH3F, TF3, TRandom3, TCanvas, TFile, gRandom

from bumphunter.bumphunter import BumpHunter3D
from bumphunter.utils import MultiOutstream


def make_histo(title, args, include_signal=True):
    """
    Make and return a histogram describe by the above parameters
    """
    histo = TH3F(
        title, "3D BumpHunter Test %s" % title,
        args.nbins, 0, args.nbins*args.binsize,
        args.nbins, 0, args.nbins*args.binsize,
        args.nbins, 0, args.nbins*args.binsize,
    )
    histo.GetXaxis().SetTitle("X")
    histo.GetXaxis().SetTitleOffset(2)
    histo.GetYaxis().SetTitle("Y")
    histo.GetYaxis().SetTitleOffset(2)
    histo.GetZaxis().SetTitle("Z")
    histo.GetZaxis().SetTitleOffset(2)
    
    for i in range(args.nbkg):
        histo.Fill(gRandom.Exp(args.bkg_mean),
                   gRandom.Exp(args.bkg_mean),
                   gRandom.Exp(args.bkg_mean))
    if include_signal:
        for i in range(args.nsig):
            x = gRandom.Gaus(args.sig_x, args.sig_spread_x)
            y = gRandom.Gaus(args.sig_y, args.sig_spread_y)
            z = gRandom.Gaus(args.sig_z, args.sig_spread_z)
            histo.Fill(x, y, z)
    return histo


def get_fit_fcn(histo):
    fit_fcn = TF3(
        "expo3", "[0]*exp(-[1] - [2]*x - [3]*y - [4]*z)",
        histo.GetXaxis().GetXmin(),
        histo.GetXaxis().GetXmax(),
        histo.GetYaxis().GetXmin(),
        histo.GetYaxis().GetXmax(),
        histo.GetZaxis().GetXmin(),
        histo.GetZaxis().GetXmax(),
    )
    fit_fcn.SetNpx(histo.GetNbinsX())
    fit_fcn.SetNpy(histo.GetNbinsY())
    fit_fcn.SetNpz(histo.GetNbinsZ())

    fit_fcn.SetParameter(0, 1000)
    fit_fcn.SetParameter(1, 1)
    fit_fcn.SetParameter(2, 0.2)
    fit_fcn.SetParameter(3, 0.2)
    fit_fcn.SetParameter(4, 0.2)

    return fit_fcn


def main(args):
    timestamp = datetime.now().isoformat()


    histo = make_histo('signal', args, include_signal=True)
    bh = BumpHunter3D(histo, fit_fcn=get_fit_fcn(histo))


    f = TFile("runs/%s.root" % timestamp, "NEW")

    c2 = TCanvas("c2")
    bh.histo.Draw("BOX")
    c3 = TCanvas("c3")
    bh.bkg_histo.Draw("BOX")
    bh.write_rootfile()

    f.Close()

    with open("runs/%s.stdout" % timestamp, "w") as outfile:
        out = MultiOutstream(sys.stdout, outfile)

        print >>out, timestamp

        import json
        print >>out, json.dumps(args.__dict__, indent=4, sort_keys=True)

        print >>out, "\nRunning BumpHunter on test data..."

        t, best_p, best_leftedge, best_width = bh.get_best_bump()

        print >>out, "Done"
        print >>out, "----------------------------------"
        # print "P-value: p = %s" % best_p # this is not very meaningful
        print >>out, "test statistic: t = %s" % t
        print >>out, "Left edge: (%s, %s, %s):" % ( (best_leftedge[0] - 1) * args.binsize,
                                                    (best_leftedge[1] - 1) * args.binsize,
                                                    (best_leftedge[2] - 1) * args.binsize )
        print >>out, "Width: (%s, %s, %s):" % (best_width[0]*args.binsize, best_width[1]*args.binsize, best_width[2]*args.binsize)
        print >>out, ""

    # Hang the program so the user can look at the output
    raw_input("Press enter to quit")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run BumpHunter on randomly generated 2D data")
    parser.add_argument('--nbkg', type=int, default=1000000,
                        help='Number of background events to generate')
    parser.add_argument('--nsig', type=int, default=200,
                        help='Number of signal events to generate')
    parser.add_argument('--nbins', type=int, default=20,
                        help='Number of bins along each dimension')
    parser.add_argument('--binsize', type=float, default=0.5,
                        help='Size of each bin along each dimension')
    parser.add_argument('--bkg-mean', type=float, default=4,
                        help='Mean of background exponential distribution')
    parser.add_argument('--sig-x', type=float, default=7,
                        help='x position of signal peak')
    parser.add_argument('--sig-y', type=float, default=7,
                        help='y position of signal peak')
    parser.add_argument('--sig-z', type=float, default=7,
                        help='z position of signal peak')
    parser.add_argument('--sig-spread-x', type=float, default=.6,
                        help='stdev along x axis of signal peak')
    parser.add_argument('--sig-spread-y', type=float, default=.2,
                        help='stdev along y axis of signal peak')
    parser.add_argument('--sig-spread-z', type=float, default=.2,
                        help='stdev along z axis of signal peak')
    parser.add_argument('--num-pseudo', type=int, default=10,
                        help='max number of pseudoexperiments to run to estimate final p-val')
    parser.add_argument('--pval', type=float, required=False,
                        help='run pseudoexperiments until this pval is reached, not exceeding '
                        '--num-pseudo')

    args = parser.parse_args()
    main(args)
