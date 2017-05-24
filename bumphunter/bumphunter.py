"""
Multivariate implementation of BumpHunter hypothesis test

Author: Vyassa Baratham <vbaratham@berkeley.edu>
"""

import math

from ROOT import TMath

class BumpHunter(object):
    """
    Superclass for BumpHunters. Currently only contains static constants (enums that I'm
    too lazy to make into real enums) and helper functions. In the future, may contain an
    implementation of the 1D BumpHunter (or maybe that should go in a new class
    BumpHunter1D, I'm not going to think too hard about it right now)
    """
    WINDOW_DEF_RECTANGLE = "WINDOW_DEF_RECTANGLE"
    SIDEBAND_DEF_RECTANGLE = "SIDEBAND_DEF_RECTANGLE"
    SIEBAND_DEF_NONE = "SIDEBAND_DEF_NONE"

    
    def integrate_histo_with_error(bins):
        """
        Add the given bins from self.histo and their errors, return both sums.

        bins - collection of Global bin numbers of bins to integrate
        """
        return (
            sum(self.histo.GetBinContent(b) for b in bins),
            math.sqrt(sum(self.histo.GetBinError(b)*self.histo.GetBinError(b) for b in bins))
        )


    def mean_poisson_pval(d, b, b_error, conv_width=3, step_size=1):
        """
        Mostly transcribed from
        https://svnweb.cern.ch/trac/atlasoff/browser/Trigger/TrigFTK/SuPlot/trunk/src/bumphunter/StatisticsAnalysis.C
        (including the docstring)

        Convolve a Gaussian (non-negative part) with a Poisson. Background is b+-deltaB, and
        we need the mean Poisson probability to observe at least d (if d >=b, else at most d), given this PDF for b.
        The way is to cut the PDF or b into segments, whose area is exactly calculable, take the
        Poisson probability at the center of each gaussian slice, and average the probabilities
        using the area of each slice as weight.

        d - data counts
        b - background counts
        b_error - error in bkg counts
        conv_width - range of l in convolution loop (I think this is how many sigmas of to convolve)
        step_size - step size for convolution

        TODO: I think there might be a better way to do this using ROOT.TMath
        """
        if b_error == 0:
            return TMath.Gamma(d, b) # Guaranteed to have d > b, so this is equivalent to commonFunctions.h:503

        # TODO: Pythonify
        mean, total_weight = 0.0, 0.0
        l = -conv_width
        while l <= conv_width:
            bcenter = max(0, b + l*b_error)
            this_slice_weight = TMath.normal_cdf(l + 0.5*step_size) - TMath.normal_cdf(l - 0.5*step_size)
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
    ):
        """
        histo - TH2 containing the data to analyze
        bkg_histo - background histo (null hypothesis) to be used for all "psuedoexperiments"
        window_def - how to define the window shape (one of BumpHunter.WINDOW_DEF_*)
        sideband_def - how to define the sideband shape (one of BumpHunter.SIDEBAND_DEF_*)
        sideband_req - the maximum pval for the sideband that we accept (higher: more
                       strict prohibition of excesses in the sideband relative to bkg)
        """
        self.histo = histo
        self.bkg_histo = bkg_histo or BumpHunter2D.bkg_histo(histo)
        self.window_def = window_def
        self.sideband_def = sideband_def
    
    def bkg_histo(data_histo):
        """
        Return a histogram to be used as background (null hypothesis). When
        Implemented, it should ideally use the process described
        in sec 2.1.1 of arxiv.org/pdf/1101.0390.pdf, or in sec 8.1 of
        https://cds.cern.ch/record/2151829/files/ATL-COM-PHYS-2016-471.pdf,
        or at least just fit an exponential to the data histo and fill a new histo
        with that distribution.
        """
        raise NotImplementedError("deriving background from data not yet implemented")

    def central_window_widths(self):
        """
        Return an iterable containing (x, y) tuples of central window widths to use.
        Implements default behavior of [1, floor(N)] for each dimension
        """
        return zip(
            xrange(1, math.floor(self.histo.GetNbinsX())),
            xrange(1, math.floor(self.histo.GetNbinsY()))
        )


    def central_window(self, center, width):
        """
        Return a collection of Global bin numbers from self.histo within the window at position @center,
        with width @width. Exactly how the window is defined depends on self.window_def

        center - (x, y) tuple of central window coordinate
        width - (x, y) tuple of central window width
        """
        if self.window_def == BumpHunter.WINDOW_DEF_RECTANGLE:
            return (
                self.histo.GetBin(center[0] + x, center[1] + y)
                for x in range(-width[0], width[0])
                for y in range(-width[1], width[1])
            )
        else:
            raise ValueError("Unrecognized window definition (self.window_def). "
                             "Use one of BumpHunter.WINDOW_DEF_*")


    def sideband_width(self, central_width):
        """
        Return the width along each dimension of the sideband to be used for a
        particular central window.
        This is a single-valued function of central window width. TODO: is this ok?
        Do we need to use different sideband widths for different center locations
        with the same central_width?

        central_width - width of central window
        """
        return tuple(max(1, math.floor(x/2)) for x in central_width)


    def sideband(self, center, central_width, sideband_width, central_window):
        """
        Return a collection of the Global bin numbers from self.histo corresponding
        to the sideband of a particular central window. Exactly how the sideband is
        defined depends on self.sideband_def

        center -
        central_width -
        sideband_width - 
        """
        if self.sideband_def != BumpHunter.SIDEBAND_DEF_NONE and self.window_def != BumpHunter.WINDOW_DEF_RECTANGLE:
            raise NotImplementedError("Sidebands not yet implemented for "
                                      "non-rectangular central window")
        if self.sideband_def == BumpHunter.SIDEBAND_DEF_RECTANGLE:
            tot_width = central_width + sideband_width
            return (
                cell for cell in set(
                    self.histo.GetBin(center[0] + x, center[1] + y)
                    for x in range(-tot_width[0], tot_width[0])
                    for y in range(-tot_width[1], tot_width[1])
                ) if cell not in central_window
            )  ## TODO: kinda inefficient^
        if self.sideband_def == BumpHunter.SIDEBAND_DEF_NONE:
            return []
        else:
            raise ValueError("Unrecognized sideband definition (self.sideband_def). "
                             "Use one of BumpHunter.SIDEBAND_DEF_*")


    def step_size(self, central_width):
        """
        Return the step size along each dimension to use for the given central width
        """
        return tuple(max(1, math.floor(x/2)) for x in central_width)


    def central_windows_and_sidebands(self, central_width, sideband_width):
        """
        Yield tuples of (center, central_window, sideband), where central_window
        and sideband are collections of Global bin numbers. Each window
        must be defined so that the window+sieband is in the histogram.
        """
        minctr = central_width + sideband_width
        maxctr = (self.histo.GetNbinsX(), self.histo.GetNbinsY()) - minctr
        step_size = self.step_size(central_width)
        for center in zip(range(minctr[0], maxctr[0], step_size[0]),
                          range(minctr[1], maxctr[1], step_size[1])):
            central_window = self.central_window(center, central_width)
            sideband = self.sideband(center, central_width, sideband_width, central_window)
            yield (center, central_window, sideband)


    def get_statistic(self):
        """
        Compute the BumpHunter test statistic
        """
        return -math.log(min(self.pvals))


    def pvals(self):
        """
        Iterate over widths and center locations, compute
        the background (null hypothesis),and calculate p-values for each. Yield the
        ones that are interesting (a p-val is interesting if it is less than 1), along with their center and width
        """
        for central_width in self.central_window_widths():
            sideband_width = self.sideband_width(central_width)
            for center, central_window, sideband in self.central_windows_and_sidebands(central_width, sideband_width):
                # At this point, central_window and sideband are collections of
                # Global bin numbers
                bkg_histo = self.bkg_histo #or self.bkg_histo(central_window, sideband)

                # Integrate histograms over the central window and sidebands
                # Note we will not be able to use TH2.Integrate() for non rectangular windows,
                # so we don't bother using it here even though it would presumably be faster
                dc, dc_error = BumpHunter.integrate_histo_with_error(self.histo, central_window)  # counts/error of data in the central window
                bc, bc_error = BumpHunter.integrate_histo_with_error(bkg_histo, central_window)  # counts/error of bkg in the central window
                ds, ds_error = BumpHunter.integrate_histo_with_error(self.histo, sideband)  # counts/error of data in the sideband
                bs, bs_error = BumpHunter.integrate_histo_with_error(bkg_histo, sideband)  # counts/error of bkg in the sideband

                if bkg_central_window_sum == 0:
                    continue  # TODO: Is this necessary?
                if data_central_window_sum <= bkg_central_window_sum:
                    continue  # Not an excess

                # TODO: implement deficit detection

                central_window_p = BumpHunter.mean_poisson_pval(dc, bc, bc_error)
                sideband_p = BumpHunter.mean_poisson_pval(ds, bs, bs_error) # TODO: check if we are not using sidebands

                if sideband_p < self.sideband_req:
                    continue

                # We return central_window_p, rather than eq (17) in
                # arxiv.org/pdf/1101.0390.pdf, in light of the paragraph
                # after that equation
                yield central_window_p

