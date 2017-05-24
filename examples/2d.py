"""
Creates a 2d histogram consisting of a Gaussian signal and exponentially decaying background,
runs 2d bumphunter, and shows results in ROOT windows
"""

from ROOT import TH2F, TRandom3, TCanvas

NUM_BINS = 40  # per dimension
BIN_SIZE = 0.5
NUM_BKG_EVENTS = 100000
NUM_SIG_EVENTS = 300
BKG_MEAN = 8
SIG_CENTER_X = 5
SIG_CENTER_Y = 12
SIG_SPREAD_X = .5
SIG_SPREAD_Y = .2

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
    histo.Draw("LEGO2")

    # The current behavior of BumpHunter2D.bkg_histo() is stupid,
    # so create a background histo to pass in
    bkg_histo = make_histo('bkg', include_signal=False)
    c2 = TCanvas('c2')
    bkg_histo.Draw("LEGO2")

    from bumphunter import BumpHunter2D
    b = BumpHunter2D(histo, bkg_histo=bkg_histo)
    print b.get_statistic()

    # Hang the program so the user can look at the output
    raw_input("Press enter to quit")
