"""
Creates a 2d histogram consisting of a Gaussian signal and exponentially decaying background,
runs 2d bumphunter, and shows results in ROOT windows
"""

import sys
import math
from ROOT import TH2F, TF2, TRandom3, TCanvas

NUM_BINS = 40  # per dimension
BIN_SIZE = 0.5 # square bins for simplicity
NUM_BKG_EVENTS = 100000
NUM_SIG_EVENTS = 100
BKG_MEAN = 8
SIG_CENTER_X = 5
SIG_CENTER_Y = 12
SIG_SPREAD_X = .5
SIG_SPREAD_Y = .2
NUM_PSEUDOEXPERIMENTS = 10


def make_histo(title, include_signal=True):
    """
    Make and return a histogram describe by the above parameters
    """
    histo = TH2F(title, "2D BumpHunter Test %s" % title, NUM_BINS, 0, NUM_BINS*BIN_SIZE, NUM_BINS, 0, NUM_BINS*BIN_SIZE)
    histo.GetXaxis().SetTitle("X")
    histo.GetXaxis().SetTitleOffset(2)
    histo.GetYaxis().SetTitle("Y")
    histo.GetYaxis().SetTitleOffset(2)
    
    rnd = TRandom3()
    rnd.SetSeed(0)
    for i in range(NUM_BKG_EVENTS):
        histo.Fill(rnd.Exp(BKG_MEAN), rnd.Exp(BKG_MEAN))
    if include_signal:
        for i in range(NUM_SIG_EVENTS):
            x, y = rnd.Gaus(SIG_CENTER_X, SIG_SPREAD_X), rnd.Gaus(SIG_CENTER_Y, SIG_SPREAD_Y)
            # print x, y
            histo.Fill(x, y)

    return histo


if __name__ == '__main__':
    histo = make_histo('signal')
    # histo.Draw("LEGO2")

    # bkg_histo = make_histo('bkg', include_signal=False)
    # c2 = TCanvas('c2')
    # bkg_histo.Draw("LEGO2")

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

    from bumphunter import BumpHunter2D
    bh = BumpHunter2D(histo, fit_fcn=fit_fcn)

    c2 = TCanvas("c2")
    bh.histo.Draw("LEGO2 HIST")
    c3 = TCanvas("c3")
    bh.bkg_histo.Draw("LEGO2 HIST")

    t, best_p, best_center, best_width = bh.get_best_bump()

    print "Ran 1 BumpHunter on the test data:"
    print "----------------------------------"
    print "P-value: p = %s" % best_p
    print "test statistic: t = %s" % t
    print "Center: (%s, %s):" % (best_center[0]*BIN_SIZE, best_center[1]*BIN_SIZE)
    print "Width: (%s, %s):" % (best_width[0]*BIN_SIZE, best_width[1]*BIN_SIZE)

    print "Running %s pseudoexperiments..." % NUM_PSEUDOEXPERIMENTS
    bh.pseudoexperiments(NUM_PSEUDOEXPERIMENTS, progress_out=sys.stdout)
    with open('pseudoexperiments_t.txt', 'w') as out:
        print >>out, '\n'.join(str(t) for t in bh.pseudoexperiments_t)
    final_pval, err = bh.final_pval()
    print "Final p-value: p = %s \pm %s" % (final_pval, err)

    # Hang the program so the user can look at the output
    raw_input("Press enter to quit")
