import sys
import argparse
import math
from datetime import datetime

from ROOT import TH1D, TH2F, TF2, TRandom3, TCanvas, TFile

from bumphunter.bumphunter import BumpHunterConfig, BumpHunter, BumpHunter2D
from bumphunter.utils import MultiOutstream

class BumpHunter1D(BumpHunter2D):
    """
    Special case of the 2D BumpHunter for 1D data. The 2D fit requires at least 4 points
    along each dimensions (and accordingly produces a bkg histo with 5 y bins),
    so this class hacks everything needed to deal with that. Expects a histo with all data
    in the first bin
    """

    def window_widths(self):
        max_x = int(math.floor(self.histo.GetNbinsX()/2))
        return ((x, 2) for x in range(max_x))


    def window(self, leftedge, width):
        ## The histo has data in all 5 y bins, but take only y=1 to avoid quintuple-counting
        return tuple((leftedge[0] + x, 1) for x in range(width[0]))



# End class BumpHunter1D    


def make_histo_1d(title, args, include_signal=True):
    """
    Make and return a 1D TH2F. All the data will have a y coordinate of args.binsize/2.0
    (to keep the data in the first bin)
    but the fit needs at least 4 points in y, so we make it 5 bins wide along y
    """
    histo = TH2F(title, "1D BumpHunter2D Test %s" % title, args.nbins, 0,
                 args.nbins*args.binsize, 5, 0, 5*args.binsize)
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


def get_fit_fcn_1d(histo):
    fit_fcn = TF2(
        "expo1", "[0]*exp(-[1] - [2]*x + 0*y)",
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

    return fit_fcn


def main(args):
    timestamp = datetime.now().isoformat()

    histo = make_histo_1d('signal', args, include_signal=True)
    bh = BumpHunter1D(
        histo,
        fit_fcn=get_fit_fcn_1d(histo),
        config=BumpHunterConfig(
            allow_oob_y=True,
            sideband_def=BumpHunterConfig.SIDEBAND_DEF_NONE
        )
    )

    pvals = TH1D("pvals", "P-val for each window", 40, 0, 1)

    f = TFile("runs/%s.root" % timestamp, "NEW")

    c2 = TCanvas("c2")
    bh.histo.Draw("LEGO2 HIST")
    bh.histo.Write()
    c3 = TCanvas("c3")
    bh.bkg_histo.Draw("LEGO2 HIST")
    bh.bkg_histo.Write()

    with open("runs/%s.stdout" % timestamp, "w") as outfile:
        out = MultiOutstream(sys.stdout, outfile)

        import json
        print >>out, json.dumps(args.__dict__, indent=4, sort_keys=True)

        print >>out, "\nRunning BumpHunter on test data..."

        t, best_p, best_leftedge, best_width = bh.get_best_bump(pvals)

        print >>out, "Done"
        print >>out, "----------------------------------"
        # print "P-value: p = %s" % best_p # this is not very meaningful
        print >>out, "test statistic: t = %s" % t
        print >>out, "Left edge: (%s, %s):" % ( (best_leftedge[0] - 1) * args.binsize,
                                                (best_leftedge[1] - 1) * args.binsize )
        print >>out, "Width: (%s, %s):" % (best_width[0]*args.binsize, best_width[1]*args.binsize)
        print >>out, ""
        print >>out, "All pvals:\n",
        print >>out, '\n'.join(str(x) for x in sorted(bh.pvals()))

        if args.num_pseudo > 0:
            print >>out, "Running %s pseudoexperiments..." % args.num_pseudo
            print >>out, ""

            final_pval, err = bh.pseudoexperiments(args.num_pseudo, progress_out=out)

            print >>out, ""
            print >>out, "Final p-value: p = %s \pm %s" % (final_pval, err)


    c4 = TCanvas("c4")
    pvals.Draw()
    pvals.Write()
    f.Close()

    # Hang the program so the user can look at the output
    raw_input("Press enter to quit")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run BumpHunter on randomly generated 2D data")
    parser.add_argument('--nbkg', type=int, default=20000,
                        help='Number of background events to generate')
    parser.add_argument('--nsig', type=int, default=600,
                        help='Number of signal events to generate')
    parser.add_argument('--nbins', type=int, default=50,
                        help='Number of bins along each dimension')
    parser.add_argument('--binsize', type=float, default=0.5,
                        help='Size of each bin along each dimension')
    parser.add_argument('--bkg-mean', type=float, default=8,
                        help='Mean of background exponential distribution')
    parser.add_argument('--sig-x', type=float, default=12,
                        help='x position of signal peak')
    parser.add_argument('--sig-spread-x', type=float, default=.7,
                        help='stdev along x axis of signal peak')
    parser.add_argument('--num-pseudo', type=int, default=100,
                        help='max number of pseudoexperiments to run to estimate final p-val')
    parser.add_argument('--pval', type=float, default=0.1,
                        help='run pseudoexperiments until this pval is reached, not exceeding '
                        '--num-pseudo')
    # parser.add_argument('--rootfile', type=str,
    #                     help='take the "signal" histo from this root file instead of generating fake data')

    args = parser.parse_args()
    main(args)


