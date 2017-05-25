"""
Multivariate implementation of BumpHunter hypothesis test

Author: Vyassa Baratham <vbaratham@berkeley.edu>
"""

import math
from itertools import product
from operator import itemgetter

from ROOT import Math, TMath, TF2

class BumpHunter(object):
    """
    Superclass for BumpHunters. Currently only contains static constants (enums that I'm
    too lazy to make into real enums) and helper functions. In the future, may contain an
    implementation of the 1D BumpHunter (or maybe that should go in a new class
    BumpHunter1D - I'm not going to think too hard about it right now)
    """
    WINDOW_DEF_RECTANGLE = "WINDOW_DEF_RECTANGLE"
    SIDEBAND_DEF_RECTANGLE = "SIDEBAND_DEF_RECTANGLE"
    SIDEBAND_DEF_NONE = "SIDEBAND_DEF_NONE"

    @classmethod
    def integrate_histo_with_error(cls, histo, bins):
        """
        Add the given bins from self.histo and their errors, return both sums.

        bins - collection of (binx, biny) tuples of bins to integrate
        """
        return (
            sum(histo.GetBinContent(*b) for b in bins),
            math.sqrt(sum(histo.GetBinError(*b)*histo.GetBinError(*b) for b in bins))
        )

    @classmethod
    def mean_poisson_pval(cls, d, b, b_error, conv_width=3, step_size=1):
        """
        Mostly transcribed from
        https://svnweb.cern.ch/trac/atlasoff/browser/Trigger/TrigFTK/SuPlot/trunk/src/bumphunter/StatisticsAnalysis.C
        (including the docstring)

        Convolve a Gaussian (non-negative part) with a Poisson. Background is b+-b_error, and
        we need the mean Poisson probability to observe at least d (if d >=b, else at most d), given this PDF for b.
        The way is to cut the PDF or b into segments, whose area is exactly calculable, take the
        Poisson probability at the center of each gaussian slice, and average the probabilities
        using the area of each slice as weight.

        But here, currently, we are guaranteed to have d > b so some of this simplifies.

        d - data counts
        b - background counts
        b_error - error in bkg counts
        conv_width - range of l in convolution loop (I think this is how many sigmas of to convolve)
        step_size - step size for convolution

        TODO: I think there might be a better way to do this using ROOT.TMath
        """
        if b_error == 0:
            return TMath.Gamma(d, b) # Guaranteed to have d > b, so this is equivalent to commonFunctions.h:503

        # TODO: Pythonify the following
        mean, total_weight = 0.0, 0.0
        l = -conv_width
        while l <= conv_width:
            bcenter = max(0, b + l*b_error)
            this_slice_weight = Math.normal_cdf(l + 0.5*step_size) - Math.normal_cdf(l - 0.5*step_size)
            this_pval = TMath.Gamma(d, bcenter) # Guaranteed to have d > b, so this is equivalent to StatisticsAnalysis.C:1710
            mean += this_pval*this_slice_weight
            total_weight += this_slice_weight
            l += step_size
        return mean / total_weight


class BumpHunter2D(BumpHunter):
    """Implementation of 2D BumpHunter"""
    
    def __init__(
            self,
            histo,
            bkg_histo=None,
            window_def=BumpHunter.WINDOW_DEF_RECTANGLE,
            sideband_def=BumpHunter.SIDEBAND_DEF_RECTANGLE,
            sideband_req=1e-3,
            deviation_from_sq=4,
    ):
        """
        histo - TH2 containing the data to analyze
        bkg_histo - background histo (null hypothesis). If not specified, uses the result of
                    BumpHunter2D.bkg_histo() as bkg_histo
        window_def - how to define the window shape (one of BumpHunter.WINDOW_DEF_*)
        sideband_def - how to define the sideband shape (one of BumpHunter.SIDEBAND_DEF_*)
        sideband_req - the maximum pval for the sideband that we accept (higher: more
                       strict prohibition of excesses in the sideband relative to bkg)
        deviation_from_sq - we use windows of width that are almost square (ie, xwidth = ywidth).
                            This param indicates how far off we allow it to be
                            (ie, the max abs(xwidth-ywidth))
        """
        self.histo = histo
        self.bkg_histo = bkg_histo or BumpHunter2D.bkg_histo(histo)
        self.window_def = window_def
        self.sideband_def = sideband_def
        self.sideband_req = sideband_req
        self.deviation_from_sq = deviation_from_sq


    @classmethod
    def bkg_histo(cls, data_histo):
        """
        Return a histogram to be used as background (null hypothesis). Fits
        a decaying exponential to the data and returns a histogram representing
        the fit.
        Eventually it should ideally use either the process described
        in sec 2.1.1 of arxiv.org/pdf/1101.0390.pdf, or the one in sec 8.1 of
        https://cds.cern.ch/record/2151829/files/ATL-COM-PHYS-2016-471.pdf,
        or at least just fit an exponential to the data histo and fill a new histo
        with that distribution.
        """
        fcn = TF2(
            "expo2", "[0]*exp(-[1] - [2]*x - [3]*y)",
            data_histo.GetXaxis().GetXmin(),
            data_histo.GetXaxis().GetXmax(),
            data_histo.GetYaxis().GetXmin(),
            data_histo.GetYaxis().GetXmax(),
        )
        fcn.SetNpx(data_histo.GetNbinsX())
        fcn.SetNpy(data_histo.GetNbinsY())
        fcn.SetParameter(0, 1000)
        fcn.SetParameter(1, 1)
        fcn.SetParameter(2, 0.2)
        fcn.SetParameter(3, 0.2)
        
        fit = data_histo.Fit(fcn)
        bkg_histo = fcn.CreateHistogram()
        # TODO: set title, axis labels, etc.

        return bkg_histo


    def window_widths(self):
        """
        Return an iterable containing (x, y) tuples of window widths to scan.
        Default behavior is to start with square windows, then
        skew every window by adding and subtracting each integer in
        [1, self.deviation_from_sq) from the xwidth of each window, leaving ywidth
        untouched, returning only those windows with widths in [1, floor(N/2)].
        To see why we do this, look at the following, which shows the main diagonal
        (square windows) as X's, the "nearby" (slightly skewed windows) as 0's, and
        far away (heavily skewed windows) as dots:

          1 2 3 4 5 6 7 8 9
          -----------------
        1| X 0 0 0 . . . . .
        2| 0 X 0 0 0 . . . .
        3| 0 0 X 0 0 0 . . .
        4| 0 0 0 X 0 0 0 . .
        5| . 0 0 0 X 0 0 0 .
        6| . . 0 0 0 X 0 0 0 
        7| . . . 0 0 0 X 0 0
        8| . . . . 0 0 0 X 0
        9| . . . . . 0 0 0 X

        # TODO: allow skewing further in one direction than the other
        """
        # TODO: this could probably be prettier
        max_x = int(math.floor(self.histo.GetNbinsX()/2))
        max_y = int(math.floor(self.histo.GetNbinsY()/2))
        for xwidth, ywidth in zip(
            xrange(1, max_x),
            xrange(1, max_y),
        ):
            yield (xwidth, ywidth)
            for i in range(1, self.deviation_from_sq):
                if xwidth - i >= 1:
                    yield (xwidth - i, ywidth)
                if xwidth + i <= max_x:
                    yield (xwidth + i, ywidth)


    def window(self, center, width):
        """
        Return a collection of (binx, biny) tuples corresponding to the window
        with the given center location and width. Exactly how the window is defined
        depends on self.window_def.

        center - (x, y) tuple of window coordinate
        width - (x, y) tuple of window width
        """
        if self.window_def == BumpHunter.WINDOW_DEF_RECTANGLE:
            return tuple(
                (center[0] + x, center[1] + y)
                for x in range(-width[0], width[0])
                for y in range(-width[1], width[1])
            )        
        else:
            raise ValueError("Unrecognized window definition (self.window_def). "
                             "Use one of BumpHunter.WINDOW_DEF_*")


    def sideband_width(self, window_width):
        """
        Return (widthx, widthy) of the sideband to be used for a particular window size.
        This is a single-valued function of central window width. TODO: is this ok?
        Do we need to use different sideband widths for different center locations
        with the same window_width?

        window_width - width of central window
        """
        return tuple(max(1, int(math.floor(x/2))) for x in window_width)


    def sideband(self, center, window_width, sideband_width, window):
        """
        Return a collection of (binx, biny) tuples corresponding
        to the sideband of a particular central window. Exactly how the sideband is
        defined depends on self.sideband_def

        center - coordinates of the center of the window
        window_width - x/y widths of the central window
        sideband_width - x/y widths of the sideband
        """
        if self.sideband_def != BumpHunter.SIDEBAND_DEF_NONE and self.window_def != BumpHunter.WINDOW_DEF_RECTANGLE:
            raise NotImplementedError("Sidebands not yet implemented for "
                                      "non-rectangular central window")
        if self.sideband_def == BumpHunter.SIDEBAND_DEF_RECTANGLE:
            tot_width = window_width + sideband_width
            return tuple(
                cell for cell in set(
                    (center[0] + x, center[1] + y)
                    for x in range(-tot_width[0], tot_width[0])
                    for y in range(-tot_width[1], tot_width[1])
                ) if cell not in window
            )  ## TODO: kinda inefficient^
        if self.sideband_def == BumpHunter.SIDEBAND_DEF_NONE:
            return []
        else:
            raise ValueError("Unrecognized sideband definition (self.sideband_def). "
                             "Use one of BumpHunter.SIDEBAND_DEF_*")


    def step_size(self, window_width):
        """
        Return the step size along each dimension to use for the given window width
        """
        return tuple(max(1, int(math.floor(x/2))) for x in window_width)


    def windows_and_sidebands(self, window_width, sideband_width):
        """
        Yield tuples of (center, window, sideband), where window
        and sideband are collections of (binx, biny) tuples. Each window
        must be defined so that the window+sideband is in the histogram.
        """
        minctr = window_width + sideband_width
        maxctr = (self.histo.GetNbinsX() - minctr[0], self.histo.GetNbinsY() - minctr[1])
        step_size = self.step_size(window_width)
        for center in product(range(minctr[0], maxctr[0], step_size[0]),
                              range(minctr[1], maxctr[1], step_size[1])):
            window = self.window(center, window_width)
            sideband = self.sideband(center, window_width, sideband_width, window)
            yield (center, window, sideband)


    def get_statistic(self):
        """
        Compute the BumpHunter test statistic
        """
        best_p, best_center, best_width = self.get_best_bump()
        return -math.log(best_p)


    def get_best_bump(self):
        best_p, best_center, best_width = min(self.pvals(), key=itemgetter(0))
        return best_p, best_center, best_width


    def pvals(self):
        """
        Iterate over widths and center locations, compute
        the background (null hypothesis),and calculate p-values for each. Yield the
        ones that are interesting (a p-val is interesting if it is less than 1),
        along with their center and width
        """
        for window_width in self.window_widths():
            sideband_width = self.sideband_width(window_width)
            for center, window, sideband in self.windows_and_sidebands(window_width, sideband_width):
                # Integrate histograms over the central window and sidebands
                # Note we will not be able to use TH2.Integrate() for non rectangular windows,
                # so we don't bother using it here even though it would presumably be faster

                # counts/error of data in the central window
                dc, dc_error = BumpHunter.integrate_histo_with_error(self.histo, window)

                # counts/error of bkg in the central window
                bc, bc_error = BumpHunter.integrate_histo_with_error(self.bkg_histo, window)

                # counts/error of data in the sideband
                ds, ds_error = BumpHunter.integrate_histo_with_error(self.histo, sideband)

                # counts/error of bkg in the sideband
                bs, bs_error = BumpHunter.integrate_histo_with_error(self.bkg_histo, sideband)

                if bc == 0:
                    continue  # TODO: Is this necessary?
                if dc <= bc:
                    continue  # Not an excess

                # TODO: implement deficit detection

                window_p = BumpHunter.mean_poisson_pval(dc, bc, bc_error)
                sideband_p = BumpHunter.mean_poisson_pval(ds, bs, bs_error) # TODO: check if we are not using sidebands

                if sideband_p < self.sideband_req:
                    continue

                # We return window_p, rather than eq (17) in
                # arxiv.org/pdf/1101.0390.pdf, in light of the paragraph
                # after that equation
                yield window_p, center, window_width

