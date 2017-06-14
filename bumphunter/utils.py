from ROOT import Math, TMath

class MultiOutstream(object):
    def __init__(self, *streams):
        self.streams = streams

    def write(self, s):
        for stream in self.streams:
            stream.write(s)

def mean_poisson_pval(d, b, b_error, conv_width=3, step_size=1):
    """
    Mostly transcribed from
    https://svnweb.cern.ch/trac/atlasoff/browser/Trigger/TrigFTK/SuPlot/trunk/src/bumphunter/StatisticsAnalysis.C
    (including the docstring)

    Convolve a Gaussian (non-negative part) with a Poisson. Background is b+-b_error, and
    we need the mean Poisson probability to observe at least d (if d >=b, else at most d),
    given this PDF for b.
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
        this_pval = poisson_pval(d, bcenter)
        mean += this_pval*this_slice_weight
        total_weight += this_slice_weight
        l += step_size
    return mean / total_weight

def poisson_pval(d, b):
    """ Two sided only """
    return TMath.Gamma(d, b) if d >= b else 1 - TMath.Gamma(d+1, b)

