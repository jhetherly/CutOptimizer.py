"""\package docstring
This module holds the CutOptimizer class used in optimizing cuts as well as
several test statistic functions.

This module relies on the parent module to import the pyROOT ("import ROOT")
library.

"""

#import ROOT
import math

def ATLAS_test_statistic(s, b):
    """Gives the standard ATLAS test statistic.
    
    \param[in] s the signal value
    \param[in] b the background value
    \return value of test statistic

    """
    return math.sqrt(2.0*((s+b)*math.log(1.0+s/b)-s))

class CutOptimizer:
    """Class that aids in optimizing cuts given signal and background 
    histograms.

    """
    def __init__(self, f = None):
        """Constructor requires a function passed in as the test statistic.
        
        \param[in] f test statistic funtion that should take two arguments
        (signal and background values)
        
        """
        if f == None:
            self.fitness_function = ATLAS_test_statistic
        else:
            self.fitness_function = f

    def format_range(self, array, side):
        formatted_string = ""
        for i in range(len(array)):
            if (i+1)%2 == 0:
                print "even"
        return formatted_string

    def auto_optimization(self, sig, bkg, rebin = 1, x_min = float("nan"), x_max = float("nan")):
        """Function for finding an optimized set of cuts.
        
        \param[in] sig the signal histogram (needs to have same axis as
        background)
        \param[in] bkg the background histogram
        \param[in] rebin number of bins to group together
        \param[in] x_min minimum x value to search from (not a bin number)
        \param[in] x_max maximum x value to search to (not a bin number)
        \return array of x cut values and even or oddness of cuts

        """
        max_statistic = 0.0
        max_s = 0.0
        max_b = 0.0
        data = []
        good_regions = []
        max_choice = ""
        bad_regions = []
        #optimization_order = 1
        max_x_array = []
        class region_info:
            pass
        if rebin >= 1:
            s = sig.Rebin(rebin, "s-john117")
            b = bkg.Rebin(rebin, "b-john117")
        bgr = round(.1*s.GetNbinsX())
        if math.isnan(x_min):
            x_min = s.GetXaxis().GetFirst()
        else:
            x_min = s.GetXaxis().FindBin(x_min)
        if math.isnan(x_max):
            x_max = s.GetXaxis().GetLast()
        else:
            x_max = s.GetXaxis().FindBin(x_max)
        #int_regions = [x_min, x_max]
        #compute test statistic for each bin
        for b_offset in range(x_max-x_min+1):
            data.append(region_info())
            data[b_offset].status = "NA"
            data[b_offset].s = s.Integral(b_offset+x_min, b_offset+x_min)
            data[b_offset].b = b.Integral(b_offset+x_min, b_offset+x_min)
            try:
                data[b_offset].stat = self.fitness_function(data[b_offset].s,\
                        data[b_offset].b)
            except ZeroDivisionError:
                data[b_offset].stat = 0.0
            except ValueError:
                data[b_offset].stat = 0.0
            data[b_offset].low_bin = b_offset+x_min
            data[b_offset].high_bin = b_offset+x_min
        #determine which bins that have a high test statistic
        data.sort(key=lambda x: x.stat)
        data.reverse()
        for b_offset in range(len(data)):
            if b_offset == 0:
                max_statistic = data[b_offset].stat
                max_s = data[b_offset].s
                max_b = data[b_offset].b
                data[b_offset].status = "good"
            else:
                max_s += data[b_offset].s
                max_b += data[b_offset].b
                try:
                    new_stat = self.fitness_function(max_s, max_b)
                except ZeroDivisionError:
                    new_stat = 0.0
                except ValueError:
                    new_stat = 0.0
                #if new_stat >= max_statistic:
                if new_stat > max_statistic:
                    max_statistic = new_stat
                    data[b_offset].status = "good"
                else:
                    max_s -= data[b_offset].s
                    max_b -= data[b_offset].b
                    data[b_offset].status = "bad"
        #merge "good" bins into good regions
        data.sort(key=lambda x: x.low_bin)
        for b_offset in range(len(data)):
            if data[b_offset].status == "good":
                if len(good_regions) == 0:
                    good_regions.append(data[b_offset])
                    good_regions[-1].stat = max_statistic
                else:
                    #if good_regions[-1].high_bin+1+bgr >= data[b_offset].low_bin:
                    if good_regions[-1].high_bin+1 >= data[b_offset].low_bin:
                        good_regions[-1].s += data[b_offset].s
                        good_regions[-1].b += data[b_offset].b
                        good_regions[-1].stat = max_statistic
                        good_regions[-1].high_bin = data[b_offset].high_bin
                    else:
                        good_regions.append(data[b_offset])
                        good_regions[-1].stat = max_statistic
            if data[b_offset].status == "bad":
                bad_regions.append(data[b_offset])
        
        #loop over "bad" regions to find optimal test statistic
        for i in range(len(bad_regions)):
            print i+1
#hybrid way
        #determine the integration ranges (in between the good regions)
        #if len(good_regions) > 1:
        #    if x_min+bgr >= good_regions[0].low_bin:
        #        good_regions[0].low_bin = x_min
        #    if good_regions[-1].high_bin+bgr >= x_max:
        #        good_regions[-1].high_bin = x_max
        #for i in range(len(good_regions)):
        #    int_regions.insert(-1, good_regions[i].low_bin-1)
        #    int_regions.insert(-1, good_regions[i].high_bin+1)
        #if int_regions[0] > int_regions[1]:
        #    int_regions.pop(0)
        #    int_regions.pop(0)
        #if int_regions[-1] < int_regions[-2]:
        #    int_regions.pop()
        #    int_regions.pop()
        #optimization_order = 2*len(good_regions)
        #if int_regions[0] != x_min:
        #    optimization_order -= 1
        #if int_regions[-1] != x_max:
        #    optimization_order -= 1
        #print int_regions, optimization_order
        #perform optimization search with int_regions of order 1 
        #(b/c first order rarely occurs in this algorithm naturally)
        #max_statistic, max_x_array, max_choice = self._compute_max_statistic(s, b, 1, "both", x_min, x_max)
        #perform restricted domain optimization search with int_regions of order optimization_order
        
#original way
        #if good_regions[0].low_bin == x_min:
        #    print "even"
        #else:
        #    print "odd"
        #for i in range(len(good_regions)):
        #    #compute actual integrals
        #    if bgr != 0:
        #        for b_offset in range(good_regions[i].high_bin-good_regions[i].low_bin+1):
        #            actual_s += data[b_offset+good_regions[i].low_bin-1].s
        #            actual_b += data[b_offset+good_regions[i].low_bin-1].b
        #    #fill max_x_array with actual x values
        #    lower_x = s.GetXaxis().GetBinLowEdge(good_regions[i].low_bin)
        #    upper_x = s.GetXaxis().GetBinUpEdge(good_regions[i].high_bin)
        #    if i == 0 and good_regions[i].low_bin > x_min:
        #        max_x_array.append(s.GetXaxis().GetBinLowEdge(x_min))
        #        max_choice = "even"
        #    elif i ==0:
        #        max_choice = "odd"
        #    max_x_array.append(lower_x)
        #    max_x_array.append(upper_x)
        #if max_x_array[-1] < s.GetXaxis().GetBinUpEdge(x_max):
        #    max_x_array.append(s.GetXaxis().GetBinUpEdge(x_max))
        #if bgr != 0:
        #    try:
        #        actual_statistic = self.fitness_function(actual_s, actual_b)
        #    except ZeroDivisionError:
        #        actual_statistic = 0.0
        #    except ValueError:
        #        actual_statistic = 0.0
        return max_x_array, max_choice


    def n_sided_optimization(self, sig, bkg, n, side = "both", rebin = 1, x_min = float("nan"), x_max = float("nan")):
        """Function for finding (up to) the n-sided optimized set of cuts.
        
        \param[in] sig the signal histogram (needs to have same axis as
        background)
        \param[in] bkg the background histogram
        \param[in] n order of cut (also checks all lower orders)
        \param[in] side select from "even", "odd", or "both"
        \param[in] rebin number of bins to group together
        \param[in] x_min minimum x value to search from (not a bin number)
        \param[in] x_max maximum x value to search to (not a bin number)
        \return array of x cut values and even or oddness of cuts

        """
        max_stat = 0.0
        max_x_array = []
        max_choice = ""
        side = side.lower()
        if rebin >= 1:
            s = sig.Rebin(rebin, "s-john117")
            b = bkg.Rebin(rebin, "b-john117")
        if math.isnan(x_min):
            x_min = s.GetXaxis().GetFirst()
        else:
            x_min = s.GetXaxis().FindBin(x_min)
        if math.isnan(x_max):
            x_max = s.GetXaxis().GetLast()
        else:
            x_max = s.GetXaxis().FindBin(x_max)
        for i in range(n):
            stat, x_array, choice = self._compute_max_statistic(s, b, i+1, side, x_min, x_max)
            if stat >= max_stat:
                max_stat = stat
                max_x_array = x_array
                max_choice = choice
        for i in range(len(max_x_array)):
            if i == 0:
                max_x_array[i] = s.GetXaxis().GetBinLowEdge(max_x_array[i])
            else:
                max_x_array[i] = s.GetXaxis().GetBinUpEdge(max_x_array[i])
        return max_x_array, max_choice

    
    def _move_x_array(self, x_array, x_array_initial):
        """Aids in finding the next set of bin numbers for the
        n_sided_optimization
        
        """
        x_array_status = ["NA","NA"]

        for i in range(len(x_array)-2):
            x_array_status.insert(-1, "OK")

        for i in range(len(x_array)):
            if x_array_status[i] == "NA":
                continue
            if x_array[i]+1 == x_array[i+1]:
                x_array_status[i] = "RESET"
            else:
                x_array[i] += 1
                break

        if x_array_status[-2] == "RESET":
            return False, x_array
        for i in range(len(x_array)):
            if x_array_status[i] == "NA":
                continue
            if x_array_status[i] == "RESET":
                x_array[i] = x_array_initial[i]
            else:
                break
        return True, x_array

    def _compute_max_statistic(self, sig, bkg, n, side, x_min, x_max):
        """Aids in integration and finding the max statistic in the 
        n_sided_optimization 
        
        note: upper edge is 'x', except x_min
        
        """
        max_statistic = 0.0
        max_choice = ""
        max_x_array = []
        x_array = []
        x_array.append(x_min)
        for i in range(n):
            x_array.append(x_min+i)
        x_array.append(x_max)
        x_array_initial = x_array[:]
        continue_to_move = True
        while continue_to_move == True:
            odd_stat = 0
            even_stat = 0
            sig_even_sum = 0
            sig_odd_sum = 0
            bkg_even_sum = 0
            bkg_odd_sum = 0
            sig_odd_sum += sig.Integral(x_min, x_array[1])
            bkg_odd_sum += bkg.Integral(x_min, x_array[1])
            for i, x in enumerate(x_array[1:-1]):
                if i % 2 == 0:
                    sig_even_sum += sig.Integral(x_array[i+1]+1, x_array[i+2])
                    bkg_even_sum += bkg.Integral(x_array[i+1]+1, x_array[i+2])
                else:
                    sig_odd_sum += sig.Integral(x_array[i+1]+1, x_array[i+2])
                    bkg_odd_sum += bkg.Integral(x_array[i+1]+1, x_array[i+2])
            try:
                odd_stat = self.fitness_function(sig_odd_sum, bkg_odd_sum)
            except ZeroDivisionError:
                pass
            except ValueError:
                pass
            try:
                even_stat = self.fitness_function(sig_even_sum, bkg_even_sum)
            except ZeroDivisionError:
                pass
            except ValueError:
                pass
            if side != "even":
                if odd_stat >= max_statistic:
                    max_statistic = odd_stat
                    max_x_array = x_array[:]
                    max_choice = "odd"
            if side != "odd":
                if even_stat >= max_statistic:
                    max_statistic = even_stat
                    max_x_array = x_array[:]
                    max_choice = "even"
            continue_to_move, x_array = self._move_x_array(x_array, x_array_initial)
        return max_statistic, max_x_array, max_choice
    
