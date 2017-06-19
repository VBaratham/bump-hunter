"""
Multivariate implementation of BumpHunter hypothesis test

Author: Vyassa Baratham <vbaratham@berkeley.edu>
"""

import abc
import array
import sys
import math
import json
import itertools

from operator import add, itemgetter

from ROOT import TF1, TF2, TF3, TH1F, TH2F, TH3F, TTree

from .utils import mean_poisson_pval, poisson_pval

class BumpHunterConfig(object):
    """
    Configuration options for BumpHunters. This whole object is passed to the
    BumpHunter instances used to calculate test statistics for pseudoexperiments
    """

    WINDOW_DEF_RECTANGLE = "WINDOW_DEF_RECTANGLE"
    SIDEBAND_DEF_RECTANGLE = "SIDEBAND_DEF_RECTANGLE"
    SIDEBAND_DEF_X_ONLY = "SIDEBAND_DEF_X_ONLY"
    SIDEBAND_DEF_Y_ONLY = "SIDEBAND_DEF_Y_ONLY"
    SIDEBAND_DEF_NONE = "SIDEBAND_DEF_NONE"

    def __init__(
            self,
            window_def=WINDOW_DEF_RECTANGLE,
            sideband_def=SIDEBAND_DEF_RECTANGLE,
            sideband_req=1e-3,
            deviation_from_sq=2,
            allow_oob_x=False,
            allow_oob_y=False
    ):
        """
        window_def - how to define the window shape (one of BumpHunterConfig.WINDOW_DEF_*)
        sideband_def - how to define the sideband shape (one of BumpHunterConfig.SIDEBAND_DEF_*)
        sideband_req - the maximum pval for the sideband that we accept (higher: more
                       strict prohibition of excesses in the sideband relative to bkg)
        deviation_from_sq - we use windows that are almost square (ie, xwidth = ywidth).
                            This param indicates how far off we allow it to be
                            (ie, the max abs(xwidth-ywidth))
        allow_oob_{x,y} - if True, the sideband is allowed to be partially/completely
                          out of the histogram bounds in the {x,y} dimension
        """
        self.window_def = window_def
        self.sideband_def = sideband_def
        self.sideband_req = sideband_req
        self.deviation_from_sq = deviation_from_sq
        self.allow_oob_x = allow_oob_x
        self.allow_oob_y = allow_oob_y

    def to_dict(self):
        return {name: getattr(self, name) for name in [
            'window_def', 'sideband_def', 'deviation_from_sq', 'allow_oob_x', 'allow_oob_y',
        ]}


class BumpHunter(object):
    """
    Superclass for BumpHunters. Currently only contains helper functions.
    In the future, may contain an
    implementation of the 1D BumpHunter (or maybe that should go in a new class
    BumpHunter1D - I'm not going to think too hard about it right now)
    (see examples/1d.py for a hacky implementation of 1D BumpHunter)
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def _dim(self):
        raise NotImplementedError()

    @abc.abstractproperty
    def _histo_cls(self):
        raise NotImplementedError()

    @abc.abstractproperty
    def _fcn_cls(self):
        raise NotImplementedError()

    def __init__(self, histo, fit_fcn=None, bkg_histo=None, config=None):
        """
        histo - histogram to analyze
        bkg_histo - background histo (null hypothesis)
        fit_fcn - function to use for fitting to estimate background
        config - BumpHunterConfig instance containing config options
        """
        assert (bkg_histo or fit_fcn), "Need to supply bkg_histo or fit_fcn"

        if not isinstance(histo, self._histo_cls):
            raise TypeError(
                "For %s, histo must be a %s" % (self.__class__.__name__,
                                                self._histo_cls.__name__)
            )

        self.histo = histo
        self.bkg_histo = bkg_histo or self.make_bkg_histo(histo, fit_fcn)
        self.fit_fcn = fit_fcn

        self.config = config or BumpHunterConfig()

    @abc.abstractmethod
    def copy_histo_dimensions(self, h, name, title):
        """
        Create and return a new (empty) histo with the same dimensions/binning
        as self.histo, but with the given name and title
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def windows_and_sidebands(self):
        """
        Return tuples of (leftedge, width, window, sideband) for all windows to evaluate.
        window and sideband are d-tuples of (x, y, ...)
        """
        raise NotImplementedError()


    @staticmethod
    def integrate_histo_with_error(histo, bins):
        """
        Add the given bins from self.histo and their errors, return both sums.

        bins - collection of (binx, biny) tuples of bins to integrate
        """
        return (
            sum(histo.GetBinContent(*b) for b in bins),
            math.sqrt(sum(histo.GetBinError(*b)*histo.GetBinError(*b) for b in bins))
        )

    def prettify_axes(self, histo, x="X", y="Y", z="Z", offset=2):
        """
        Put labels on axes and offset them.

        x, y, z - axis labels to use (only pass the ones relevant to your dimensions)
        offset - offset to use
        """
        if self._dim >= 1:
            histo.GetXaxis().SetTitle(x)
            histo.GetXaxis().SetTitleOffset(offset)
        if self._dim >= 2:
            histo.GetYaxis().SetTitle(y)
            histo.GetYaxis().SetTitleOffset(offset)
        if self._dim >= 3:
            histo.GetZaxis().SetTitle(z)
            histo.GetZaxis().SetTitleOffset(offset)

    def make_bkg_histo(self, data_histo, fit_fcn):
        """
        Return a histogram to be used as background (null hypothesis). Fits
        a decaying exponential to the data and returns a histogram representing
        the fit.
        Eventually it should ideally use either the process described
        in sec 2.1.1 of arxiv.org/pdf/1101.0390.pdf, or the one in sec 8.1 of
        https://cds.cern.ch/record/2151829/files/ATL-COM-PHYS-2016-471.pdf

        data_histo - histogram of data to fit
        """
        fit = data_histo.Fit(fit_fcn, 'q') # 'q' = quiet (no output)
        bkg_histo = fit_fcn.CreateHistogram()
        bkg_histo.SetNameTitle("bkg", "Background Estimate")
        self.prettify_axes(bkg_histo)

        return bkg_histo

    def get_best_bump(self, histo=None):
        """
        Compute and store the BumpHunter test statistic

        histo - optional TH1 to store pvals for each window. Mainly used for debugging.
        """
        best_p, best_leftedge, best_width = min(self.pvals(histo=histo), key=itemgetter(0))
        t = -math.log(best_p)
        return t, best_p, best_leftedge, best_width

    def pvals(self, histo=None):
        """
        Iterate over widths and left edge locations, compute
        the background (null hypothesis),and calculate p-values for each. Yield the
        ones that are interesting (a p-val is interesting if it is less than 1),
        along with their left edge and width

        histo - optional TH1 to store pvals for each window. Mainly used for debugging.
        """
        for leftedge, width, window, sideband in self.windows_and_sidebands():
            # Integrate histograms over the central window and sidebands
            # Note we will not be able to use TH2.Integrate() for non rectangular windows,
            # so we don't bother using it here even though it would presumably be faster

            # counts/error of data in the central window
            dc, dc_error = self.integrate_histo_with_error(self.histo, window)

            # counts/error of bkg in the central window
            bc, bc_error = self.integrate_histo_with_error(self.bkg_histo, window)

            # counts/error of data in the sideband
            ds, ds_error = self.integrate_histo_with_error(self.histo, sideband)

            # counts/error of bkg in the sideband
            bs, bs_error = self.integrate_histo_with_error(self.bkg_histo, sideband)

            if bc == 0:
                continue  # TODO: Is this necessary?
            if dc <= bc:
                continue  # Not an excess

            # TODO: implement deficit detection

            window_p = mean_poisson_pval(dc, bc, bc_error)

            if self.config.sideband_def != BumpHunterConfig.SIDEBAND_DEF_NONE:
                sideband_p = mean_poisson_pval(ds, bs, bs_error)

                if sideband_p < self.config.sideband_req:
                    continue

            # if window_p == 0:
            #     import pdb; pdb.set_trace()

            if histo:
                histo.Fill(window_p)

            # We return window_p, rather than the true p-val eq (17) in
            # arxiv.org/pdf/1101.0390.pdf, in light of the paragraph
            # after that equation
            yield window_p, leftedge, width

    def pseudoexperiments(self, n, fit_fcn=None):
        """
        Run pseudoexperiments and yield t statistics
        """
        fit_fcn = fit_fcn or self.fit_fcn
        assert fit_fcn

        for i in range(n):
            yield self.one_pseudoexperiment(fit_fcn)

    def one_pseudoexperiment(self, fit_fcn):
        """
        Run one pseudoexperiment, return test statistic
        """
        pseudo_histo = self.copy_histo_dimensions(
            self.histo, "pseudodata", "2D BumpHunter pseudoexperiment")
        pseudo_histo.FillRandom(self.bkg_histo, self.histo.Integral())

        # Create a bumphunter instance with this BumpHunter's config attributes
        # Don't use copy.deepcopy() because we don't want to copy self.best_p,
        # self.pseudoexperiments_t, etc.
        # We want to re-run the constructor so it re-fits the background
        bh = self.__class__(pseudo_histo, fit_fcn=fit_fcn, config=self.config)

        t = bh.get_best_bump()[0]

        # TODO: Are these lines necessary?
        bh.histo.Delete()
        bh.bkg_histo.Delete()

        return t

    def write_rootfile(self):
        """
        Put the signal and bkg histos and the config as a json string
        into the currently open rootfile.

        rootfile - writable TFile
        """
        self.histo.Write()
        self.bkg_histo.Write()

        config_str = json.dumps(self.config.to_dict())
        config_array = array.array('B', config_str)
        
        t = TTree("configtree", "BumpHunter Config")
        t.Branch("config", config_array, "config[%s]/C" % len(config_array))
        t.Fill()
        t.Write()

    @classmethod
    def from_rootfile(cls, rootfile, signal=None, bkg=None, config=None):
        """
        Create a Bumphunter from an open root file.

        rootfile - TFile to read. The signal histo should be called 'signal', the bkg histo
                   should be called 'bkg', and the config should be stored as a json string
                   in a branch called "config" in a tree called "configtree"
        signal, bkg, config - override whats in the root file
        """
        signal = signal or rootfile.Get("signal")
        bkg = bkg or rootfile.Get("bkg")
        
        config_tree = rootfile.Get("configtree")
        if config_tree:
            config_str = iter(config_tree).next().config
            config = BumpHunterConfig(**json.loads(config_str[:-1]))
            # The tree seems to read out the null terminator, hence the [:-1] above^
        else:
            # Use default config for backwards compatibility with runs from before the
            # config was written into the rootfile (prior to 6/15/2017)
            config = BumpHunterConfig()

        bh_cls = BumpHunter1D if isinstance(signal, TH1F) else BumpHunter2D if isinstance(signal, TH2F) else BumpHunter3D # TODO: can just use `cls`?

        return bh_cls(signal, bkg_histo=bkg, config=config)


class BumpHunter1D(BumpHunter):
    @property
    def _dim(self):
        return 1
    
    @property
    def _histo_cls(self):
        return TH1F

    @property
    def _fcn_cls(self):
        return TF1

    def copy_histo_dimensions(self, h, name, title):
        return TH1F(
            name, title,
            h.GetNbinsX(), h.GetXaxis().GetXMin(), h.GetXaxis().GetXmax()
        )

    def windows_and_sidebands(self):
        pass # TODO

    def window_widths(self):
        max_x = int(math.floow(self.histo.GetNbinsX()/2))
        return xrange(1, max_x)

    def window(self, leftedge, width):
        pass

    
class BumpHunter2D(BumpHunter):
    """Implementation of 2D BumpHunter"""

    @property
    def _dim(self):
        return 2

    @property
    def _histo_cls(self):
        return TH2F

    @property
    def _fcn_cls(self):
        return TF2

    def copy_histo_dimensions(self, h, name, title):
        return TH2F(
            name, title,
            h.GetNbinsX(), h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax(),
            h.GetNbinsY(), h.GetYaxis().GetXmin(), h.GetYaxis().GetXmax(),
        )

    def windows_and_sidebands(self):
        """
        Yield tuples of (leftedge, window, sideband), where window
        and sideband are collections of (binx, biny) tuples. Each window
        must be defined so that the window+sideband is in the histogram, unless
        self.config.allow_oob_x/y is True, in which case the sideband may be partially or
        completely out of bounds.
        """
        return itertools.chain.from_iterable(
            self.windows_and_sidebands_for_width(window_width, self.sideband_width(window_width))
            for window_width in self.window_widths()
        )

    def window_widths(self):
        """
        Return an iterable containing (x, y) tuples of window widths to scan.
        Default behavior is to start with square windows, then
        skew every window by adding and subtracting each integer in
        [1, self.config.deviation_from_sq) from the xwidth of each window, leaving ywidth
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
        """
        max_x = int(math.floor(self.histo.GetNbinsX()/2))
        max_y = int(math.floor(self.histo.GetNbinsY()/2))
        for xwidth, ywidth in zip(
            xrange(1, max_x),
            xrange(1, max_y),
        ):
            yield (xwidth, ywidth)
            for i in range(1, self.config.deviation_from_sq):
                if xwidth - i >= 1:
                    yield (xwidth - i, ywidth)
                if xwidth + i <= max_x:
                    yield (xwidth + i, ywidth)

    def window(self, leftedge, width):
        """
        Return a collection of (binx, biny) tuples corresponding to the window
        with the given left edge location and width. Exactly how the window is defined
        depends on self.config.window_def.

        leftedge - (x, y) tuple of window's left edge coordinate
        width - (x, y) tuple of window width
        """
        if self.config.window_def == BumpHunterConfig.WINDOW_DEF_RECTANGLE:
            return tuple(
                (leftedge[0] + x, leftedge[1] + y)
                for x in range(width[0])
                for y in range(width[1])
            )        
        else:
            raise ValueError("Unrecognized window definition (self.config.window_def). "
                             "Use one of BumpHunterConfig.WINDOW_DEF_*")

    def sideband_width(self, window_width):
        """
        Return (widthx, widthy) of the sideband to be used for a particular window size.
        This is a single-valued function of central window width. TODO: is this ok?
        Do we need to use different sideband widths for different center locations
        with the same window_width?

        window_width - width of central window
        """
        if self.config.sideband_def == BumpHunterConfig.SIDEBAND_DEF_NONE:
            return (0, 0)
        return tuple(max(1, int(math.floor(x/2))) for x in window_width)

    def sideband(self, leftedge, window_width, sideband_width, window):
        """
        Return a collection of (binx, biny) tuples corresponding
        to the sideband of a particular central window. Exactly how the sideband is
        defined depends on self.config.sideband_def

        leftedge - coordinates of the left edge of the window
        window_width - x/y widths of the central window
        sideband_width - x/y widths of the sideband
        """
        if self.config.sideband_def != BumpHunterConfig.SIDEBAND_DEF_NONE and self.config.window_def != BumpHunterConfig.WINDOW_DEF_RECTANGLE:
            raise NotImplementedError("Sidebands not yet implemented for "
                                      "non-rectangular central window")
        if self.config.sideband_def == BumpHunterConfig.SIDEBAND_DEF_RECTANGLE:
            tot_width = tuple(sum(x) for x in zip(window_width, sideband_width, sideband_width))
            return tuple(
                cell for cell in set(
                    (leftedge[0] - sideband_width[0] + x, leftedge[1] - sideband_width[1] + y)
                    for x in range(-tot_width[0], tot_width[0])
                    for y in range(-tot_width[1], tot_width[1])
                ) if cell not in window
            )  ## TODO: kinda inefficient^
        if self.config.sideband_def == BumpHunterConfig.SIDEBAND_DEF_NONE:
            return []
        else:
            raise ValueError("Unrecognized sideband definition (self.config.sideband_def). "
                             "Use one of BumpHunterConfig.SIDEBAND_DEF_*")

    def step_size(self, window_width):
        """
        Return the step size along each dimension to use for the given window width
        """
        return tuple(max(1, int(math.floor(x/2))) for x in window_width)

    def windows_and_sidebands_for_width(self, window_width, sideband_width):
        min_leftedge = list(map(add, sideband_width, (1, 1)))
        max_leftedge = list((self.histo.GetNbinsX() - min_leftedge[0] - window_width[0],
                             self.histo.GetNbinsY() - min_leftedge[1] - window_width[1]))

        if self.config.allow_oob_x:
            min_leftedge[0] -= sideband_width[0]
            max_leftedge[0] += sideband_width[0]
        if self.config.allow_oob_y:
            min_leftedge[1] -= sideband_width[1]
            max_leftedge[1] += sideband_width[1]

        step_size = self.step_size(window_width)

        for leftedge in itertools.product(range(min_leftedge[0], max_leftedge[0], step_size[0]),
                                          range(min_leftedge[1], max_leftedge[1], step_size[1])):
            window = self.window(leftedge, window_width)
            sideband = self.sideband(leftedge, window_width, sideband_width, window)
            yield (leftedge, window_width, window, sideband)


class BumpHunter3D(BumpHunter):
    @property
    def _dim(self):
        return 3

    @property
    def _histo_cls(self):
        return TH3F

    @property
    def _fcn_cls(self):
        return TH3

    def make_bkg_histo(self, h, fit_fcn):
        """
        TF3.CreateHistogram() is a stub in ROOT -__-
        """
        fit = h.Fit(fit_fcn, 'q') # 'q' = quiet (no output)
        # fit = h.Fit(fit_fcn)
        bkg_histo = self.copy_histo_dimensions(h, "bkg", "Background Estimate")
        diff_histo = TH1F("diff", "(data - fit_fcn) for each bin", 100, -30, 30) # TODO: remove (debug)
        self.prettify_axes(bkg_histo)

        dx = (h.GetXaxis().GetXmax() - h.GetXaxis().GetXmin()) / h.GetNbinsX()
        dy = (h.GetYaxis().GetXmax() - h.GetYaxis().GetXmin()) / h.GetNbinsY()
        dz = (h.GetZaxis().GetXmax() - h.GetZaxis().GetXmin()) / h.GetNbinsZ()

        xmin = h.GetXaxis().GetXmin() + 0.5*dx
        ymin = h.GetYaxis().GetXmin() + 0.5*dy
        zmin = h.GetZaxis().GetXmin() + 0.5*dz

        for x_bin in range(h.GetNbinsX()):
            for y_bin in range(h.GetNbinsY()):
                for z_bin in range(h.GetNbinsZ()):
                    x, y, z = xmin + x_bin*dx, ymin + y_bin*dy, zmin + z_bin*dz
                    bkg_est = fit_fcn.Eval(x, y, z)
                    bkg_histo.SetBinContent(x_bin+1, y_bin+1, z_bin+1, bkg_est)
                    diff_histo.Fill(h.GetBinContent(x_bin+1, y_bin+1, z_bin+1) - bkg_est)

        return bkg_histo

    def copy_histo_dimensions(self, h, name, title):
        return TH3F(
            name, title,
            h.GetNbinsX(), h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax(),
            h.GetNbinsY(), h.GetYaxis().GetXmin(), h.GetYaxis().GetXmax(),
            h.GetNbinsZ(), h.GetZaxis().GetXmin(), h.GetZaxis().GetXmax(),
        )
    
    def windows_and_sidebands(self):
        """
        Yield tuples of (leftedge, window, sideband), where window
        and sideband are collections of (binx, biny) tuples. Each window
        must be defined so that the window+sideband is in the histogram, unless
        self.config.allow_oob_x/y is True, in which case the sideband may be partially or
        completely out of bounds.
        """
        return itertools.chain.from_iterable(
            self.windows_and_sidebands_for_width(window_width, self.sideband_width(window_width))
            for window_width in self.window_widths()
        )

    def window_widths(self):
        """
        Return an iterable containing (x, y) tuples of window widths to scan.
        Default behavior is to start with square windows, then
        skew every window by adding and subtracting each integer in
        [1, self.config.deviation_from_sq) from the xwidth of each window, leaving ywidth
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
        """
        max_x = int(math.floor(self.histo.GetNbinsX()/2))
        max_y = int(math.floor(self.histo.GetNbinsY()/2))
        max_z = int(math.floor(self.histo.GetNbinsZ()/2))
        for xwidth, ywidth, zwidth in zip(
            xrange(1, max_x),
            xrange(1, max_y),
            xrange(1, max_z),
        ):
            yield (xwidth, ywidth, zwidth)
            for i in range(1, self.config.deviation_from_sq):
                for j in range(1, self.config.deviation_from_sq):
                    _xwidth, xwidth_ = xwidth - i, xwidth + i
                    _ywidth, ywidth_ = ywidth - j, ywidth + j
                    if _xwidth >= 1 and _ywidth >= 1:
                        yield (_xwidth, _ywidth, zwidth)
                    if xwidth_ <= max_x and _ywidth >= 1:
                        yield (xwidth_, _ywidth, zwidth)
                    if _xwidth >= 1 and ywidth_ <= max_y:
                        yield (xwidth_, ywidth_, zwidth)
                    if xwidth_ <= max_x and ywidth_ <= max_y:
                        yield (xwidth_, ywidth_, zwidth)

    def window(self, leftedge, width):
        """
        Return a collection of (binx, biny) tuples corresponding to the window
        with the given left edge location and width. Exactly how the window is defined
        depends on self.config.window_def.

        leftedge - (x, y, z) tuple of window's left edge coordinate
        width - (x, y, z) tuple of window width
        """
        if self.config.window_def == BumpHunterConfig.WINDOW_DEF_RECTANGLE:
            return tuple(
                (leftedge[0] + x, leftedge[1] + y, leftedge[2] + z)
                for x in range(width[0])
                for y in range(width[1])
                for z in range(width[2])
            )        
        else:
            raise ValueError("Unrecognized window definition (self.config.window_def). "
                             "Use one of BumpHunterConfig.WINDOW_DEF_*")

    def sideband_width(self, window_width):
        # TODO: exactly same as 2D, refactor
        """
        Return (widthx, widthy, widthz) of the sideband to be used for a particular window size.
        This is a single-valued function of central window width. TODO: is this ok?
        Do we need to use different sideband widths for different center locations
        with the same window_width?

        window_width - width of central window
        """
        if self.config.sideband_def == BumpHunterConfig.SIDEBAND_DEF_NONE:
            return (0, 0, 0)
        return tuple(max(1, int(math.floor(x/2))) for x in window_width)

    def sideband(self, leftedge, window_width, sideband_width, window):
        """
        Return a collection of (binx, biny, binz) tuples corresponding
        to the sideband of a particular central window. Exactly how the sideband is
        defined depends on self.config.sideband_def

        leftedge - coordinates of the left edge of the window
        window_width - x/y/z widths of the central window
        sideband_width - x/y/z widths of the sideband
        """
        if self.config.sideband_def != BumpHunterConfig.SIDEBAND_DEF_NONE and self.config.window_def != BumpHunterConfig.WINDOW_DEF_RECTANGLE:
            raise NotImplementedError("Sidebands not yet implemented for "
                                      "non-rectangular central window")
        if self.config.sideband_def == BumpHunterConfig.SIDEBAND_DEF_RECTANGLE:
            tot_width = tuple(sum(x) for x in zip(window_width, sideband_width, sideband_width))
            return tuple(
                cell for cell in set(
                    (
                        leftedge[0] - sideband_width[0] + x,
                        leftedge[1] - sideband_width[1] + y,
                        leftedge[2] - sideband_width[2] + z,
                    )
                    for x in range(-tot_width[0], tot_width[0])
                    for y in range(-tot_width[1], tot_width[1])
                    for z in range(-tot_width[2], tot_width[2])
                ) if cell not in window
            )  ## TODO: kinda inefficient^
        if self.config.sideband_def == BumpHunterConfig.SIDEBAND_DEF_NONE:
            return []
        else:
            raise ValueError("Unrecognized sideband definition (self.config.sideband_def). "
                             "Use one of BumpHunterConfig.SIDEBAND_DEF_*")

    def step_size(self, window_width):
        # TODO exactly same as 2D, refactor
        """
        Return the step size along each dimension to use for the given window width
        """
        return tuple(max(1, int(math.floor(x/2))) for x in window_width)

    def windows_and_sidebands_for_width(self, window_width, sideband_width):
        min_leftedge = list(map(add, sideband_width, (1, 1, 1)))
        max_leftedge = list((self.histo.GetNbinsX() - min_leftedge[0] - window_width[0],
                             self.histo.GetNbinsY() - min_leftedge[1] - window_width[1],
                             self.histo.GetNbinsZ() - min_leftedge[2] - window_width[2]))

        if self.config.allow_oob_x:
            min_leftedge[0] -= sideband_width[0]
            max_leftedge[0] += sideband_width[0]
        if self.config.allow_oob_y:
            min_leftedge[1] -= sideband_width[1]
            max_leftedge[1] += sideband_width[1]
        # TODO: implement allow_oob_z (why would anyone need this?)

        step_size = self.step_size(window_width)

        for leftedge in itertools.product(range(min_leftedge[0], max_leftedge[0], step_size[0]),
                                          range(min_leftedge[1], max_leftedge[1], step_size[1]),
                                          range(min_leftedge[2], max_leftedge[2], step_size[2])):
            window = self.window(leftedge, window_width)
            sideband = self.sideband(leftedge, window_width, sideband_width, window)
            yield (leftedge, window_width, window, sideband)
