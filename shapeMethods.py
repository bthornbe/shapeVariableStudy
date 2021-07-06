import ROOT
import numpy as np


def from1dhisttoarray(hist, size):
    #!I probably don't need this method
    # returns a numpy array with the same information as hist
    # hist is a root histogram
    # size is the number of bins in the histogram
    array= np.array([])
    for i in range(size):
        array[i] = hist.GetBinContent(i)
    return array
