"""
Multivariate implementation of BumpHunter hypothesis test

Author: Vyassa Baratham <vbaratham@berkeley.edu>
"""

import math

from ROOT import Math

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
            sum(self.histo.GetBinError(b) for b in bins)
        )
        

class BumpHunter2D(BumpHunter):
    

    def __init__(
            self, histo,
            bkg_histo=None,
            window_def=BumpHunter.WINDOW_DEF_RECTANGLE,
            sideband_def=BumpHunter.SIDEBAND_DEF_RECTANGLE
    ):
        """
        histo - TH2 containing the data to analyze
        bkg_histo - background histo (null hypothesis) to be used for all "psuedoexperiments"
        bkg_fcn - callable that takes arguments (data histo, central_window, sideband)
                  for each "pseudoexperiment" and 
        """
        self.histo = histo
        self.bkg_histo = bkg_histo
        self.bkg_fcn = bkg_fcn
        self.window_def = window_def
        self.sideband_def = sideband_def
    

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


    def bkg_histo(self, central_window, sideband):
        """
        Return a histogram to be used as background (null hypothesis)
        for examining the given central_window and sideband. If self.bkg_histo
        is defined, just uses that.

        central_window - list of Global bin numbers of the central window
        sideband - list of Global bin numbers of the sideband
        """
        if self.bkg_histo:
            return self.bkg_histo
        else:
            # TODO: Implement
            raise NotImplementedError("deriving background from data not yet implemented")


    def get_statistic(self):
        """
        Compute the BumpHunter test statistic
        """
        pass

    def get_p_vals(self):
        """
        Iterate over widths and center locations ("pseudoexperiments"), compute
        the background (null hypothesis),and calculate p-values for each. Yield the
        ones that are interesting (a p-val is interesting if it is less than 1)
        """
        for central_width in self.central_window_widths():
            sideband_width = self.sideband_width(central_width)
            for center, central_window, sideband in self.central_windows_and_sidebands(central_width, sideband_width):
                # At this point, central_window and sideband are collections of
                # Global bin numbers
                bkg_histo = self.bkg_histo or self.bkg_fcn(self.histo, central_window, sideband)
                data_central_window_sum, data_central_window_error = BumpHunter.integrate_histo_with_error(self.histo, central_window)
                bkg_central_window_sum, bkg_central_window_error = BumpHunter.integrate_histo_with_error(bkg_histo, central_window)
                data_sideband_sum, data_sideband_error = BumpHunter.integrate_histo_with_error(self.histo, sideband)
                bkg_sideband_sum, bkg_sideband_error = BumpHunter.integrate_histo_with_error(bkg_histo, sideband)

                if bkg_central_window_sum == 0:
                    continue  # TODO: Is this necessary?
                if data_central_window_sum <= bkg_central_window_sum:
                    continue  # Not an excess

                # TODO: implement deficit detection

                central_window_p = None
                sideband_p = None
