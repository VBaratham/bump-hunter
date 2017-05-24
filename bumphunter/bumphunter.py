"""
Multivariate implementation of BumpHunter hypothesis test

Author: Vyassa Baratham <vbaratham@berkeley.edu>
"""

import math
from itertools import izip, chain

class BumpHunter(object):
    """
    Superclass for BumpHunters. Currently only contains static constants (enums that I'm
    too lazy to make into real enums). In the future, may contain an implementatoin
    of the 1D BumpHunter (or maybe that should go in a new class BumpHunter1D, I'm not
    going to think too hard about it right now)
    """
    WINDOW_DEF_RECTANGLE = "WINDOW_DEF_RECTANGLE"
    SIDEBAND_DEF_RECTANGLE = "SIDEBAND_DEF_RECTANGLE"
    SIEBAND_DEF_NONE = "SIDEBAND_DEF_NONE"

class BumpHunter2D(BumpHunter):
    

    def __init__(
            self, histo,
            window_def=BumpHunter.WINDOW_DEF_RECTANGLE,
            sideband_def=BumpHunter.SIDEBAND_DEF_RECTANGLE
    ):
        """
        histo - TH2 containing the data to analyze
        """
        self.histo = histo
        self.window_def = window_def
        self.sideband_def = sideband_def
    

    def central_window_widths(self):
        """
        Return an iterable containing (x, y) tuples of central window widths to use.
        Implements default behavior of [1, floor(N)] for each dimension
        """
        return izip(
            xrange(1, math.floor(self.histo.GetNbinsX())),
            xrange(1, math.floor(self.histo.GetNbinsY()))
        )


    def central_window(self, center, width):
        """
        Return a collection of cells from self.histo within the window at position @center,
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
            raise ValueError("Unrecognized window definition (self.window_def). Use one of BumpHunter.WINDOW_DEF_*")


    def central_windows(self, width):
        """
        Return a collection of central windows with the given width (a central window
        is the output of self.central_window, a collection of cells from self.histo)
        """
        

    def sideband_width(self, central_width):
        """
        Return the width of the sideband to be used for a particular central window.
        In the BumpHunter algorithm, this should be a single-valued function of
        central window width.
        Implements default behavior of max(1, floor(central_width/2)) for each dimension

        central_width - width of central window
        """
        return tuple(max(1, math.floor(x/2)) for x in central_width)


    def sideband(self, center, central_width, sideband_width):
        """
        Return a collection of the cells from self.histo corresponding to the sideband of a particular
        central window. Exactly how the sideband is defined depends on self.sideband_def

        center -
        central_width -
        sideband_width - 
        """
        if self.sideband_def != BumpHunter.SIDEBAND_DEF_NONE and self.window_def != BumpHunter.WINDOW_DEF_RECTANGLE:
            raise NotImplementedError("Sidebands not yet implemented for non-rectangular central window")
        if self.sideband_def == BumpHunter.SIDEBAND_DEF_RECTANGLE:
            central_window = set(self.central_window(center, central_width))  ## TODO: don't recompute central_window here
            tot_width = central_width + sideband_width
            return (
                cell for cell in set(
                    self.histo.GetBin(center[0] + x, center[1] + y)
                    for x in range(-tot_width[0], tot_width[0])
                    for y in range(-tot_width[1], tot_width[1])
                ) if cell not in central_window
            )  ## TODO: inefficient^
        if self.sideband_def == BumpHunter.SIDEBAND_DEF_NONE:
            return []
        else:
            raise ValueError("Unrecognized sideband definition (self.sideband_def). Use one of BumpHunter.SIDEBAND_DEF_*")


    def GetStatistic(self):
        """
        Compute the BumpHunter test statistic
        """
        for width in self.central_window_widths():
            pass
            
