from datetime import datetime
import sys
import math
import argparse

from ROOT import TH2F, TH1F, TF2, TRandom3, TFile

from bumphunter.bumphunter import BumpHunter2D, BumpHunterConfig
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

def make_histo_1d(histo2d, args):
    "Collapse a 2d histo into the first bin"
    histo = TH2F("signal1d", "1D projection of signal", args.nbins, 0,
                 args.nbins*args.binsize, 5, 0, 5*args.binsize)
    histo.GetXaxis().SetTitle("X")
    histo.GetXaxis().SetTitleOffset(2)
    histo.GetYaxis().SetTitle("Y")
    histo.GetYaxis().SetTitleOffset(2)

    histo_1d = histo2d.ProjectionX("signalprojection")
    
    for i in range(histo_1d.GetNbinsX()):
        histo.SetBinContent(histo.GetBin(i+1, 1),
                            histo_1d.GetBinContent(histo_1d.GetBin(i+1)))

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

    histo = make_histo('signal', args, include_signal=True)
    histo1d = make_histo_1d(histo, args)
    bh = BumpHunter2D(histo, fit_fcn=get_fit_fcn(histo))
    bh1d = BumpHunter1D(histo1d, fit_fcn=get_fit_fcn_1d(histo1d),config=BumpHunterConfig(
        allow_oob_y=True, sideband_def=BumpHunterConfig.SIDEBAND_DEF_NONE
    ))

    f = TFile("runs/%s.root" % timestamp, "NEW")
    bh.histo.Write()
    bh.bkg_histo.Write()
    bh1d.histo.Write()
    bh1d.bkg_histo.Write()
    f.Close()

    with open("runs/%s.stdout" % timestamp, "w") as outfile:
        out = MultiOutstream(sys.stdout, outfile)

        import json
        print >>out, json.dumps(args.__dict__, indent=4, sort_keys=True)

        print >>out, "Running BumpHunter1D on test data..."

        t, best_p, best_leftedge, best_width = bh1d.get_best_bump()

        print >>out, "Done"
        print >>out, "----------------------------------"
        print >>out, "test statistic: t = %s" % t
        print >>out, "Left edge: (%s, %s):" % ( (best_leftedge[0] - 1) * args.binsize,
                                                (best_leftedge[1] - 1) * args.binsize )
        print >>out, "Width: (%s, %s):" % (best_width[0]*args.binsize, best_width[1]*args.binsize)
        print >>out, ""

        print >>out, "Running %s pseudoexperiments..." % args.num_pseudo
        print >>out, ""

        final_pval, err = bh1d.pseudoexperiments(args.num_pseudo, progress_out=out)

        print >>out, ""
        print >>out, "Final p-value: p = %s \pm %s" % (final_pval, err)


        print >>out, ""
        print >>out, "========================================"
        print >>out, ""


        print >>out, "\nRunning BumpHunter2D on test data..."

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


        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run BumpHunter on randomly generated 2D data")
    parser.add_argument('--nbkg', type=int, default=20000,
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
    parser.add_argument('--num-pseudo', type=int, default=50,
                        help='max number of pseudoexperiments to run to estimate final p-val')
    parser.add_argument('--pval', type=float, required=False,
                        help='run pseudoexperiments until this pval is reached, not exceeding '
                        '--num-pseudo')

    args = parser.parse_args()
    main(args)
