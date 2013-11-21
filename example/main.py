import ROOT
import math
import sys
import time
sys.path.append('/Users/jhetherly/Documents/Code_and_Scripts/python/CutOptimizer')
import CutOptimizer

def test_statistic(s, b):
    return s/math.sqrt(s+b)

def main():
    """n-sided cut test"""
    start_time = time.time()
    bkg = ROOT.TH1D()
    sig = ROOT.TH1D()
    f = ROOT.TFile()
    f = ROOT.TFile.Open("BkgDRll_processed.root")
    bkg = f.Get('DR ll')
    f = ROOT.TFile.Open("SigDRll240_processed.root")
    sig = f.Get('DR ll')
    n = 5

    #co = CutOptimizer.CutOptimizer(test_statistic)
    #co = CutOptimizer.CutOptimizer()
    co = CutOptimizer.CutOptimizer(CutOptimizer.ATLAS_test_statistic)
    #x_array, choice = co.n_sided_optimization(sig, bkg, 1)
    x_array, choice = co.n_sided_optimization(sig, bkg, n)

    print x_array, choice

    end_time = time.time()
    print "It took "+str(end_time-start_time)+"s to run the entire"\
        " program for order "+str(n)+"."

if __name__ == "__main__":
    main()
