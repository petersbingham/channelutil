import pynumwrap as nw
from units import *

signs_pos = "SignsPositive"
signs_specified = "SignsSpecified"

signs_bndandres = "SignsForBndAndRes"

signs_ana_over_axis = "SignsAnaOverAxis"
signs_ana_over_thres = "SignsAnaOverThres"

# Performs asymptotic calculations
class AsymCalc:
    def __init__(self, units, ls=None, thresholds=None, sign_sel=signs_pos,
                 signs=None):
        if ls is None:
            if thresholds is not None:
                self.ls = [0]*len(thresholds)
            else:
                self.ls = [0]
        else:
            self.ls = ls

        if thresholds is None:
            self.thresholds = [0.]*len(self.ls)
        else:
            self.thresholds = thresholds

        self.units = units
        if units == rydbergs:
            self.ene_conv = 1.0
        elif units == hartrees:
            self.ene_conv = 1./rydbergs_to_hartrees
        elif units == hartrees:
            self.ene_conv = 1./rydbergs_to_eVs
        self.sign_sel = sign_sel
        self.signs = signs

    def __str__(self):
        ret = "ls:" + str(self.ls) + ", thres:" + str(self.thresholds)
        ret += ", signs:"
        if self.signs is None:
            ret += self.sign_sel
        else:
            ret += str(map(lambda x:str(x),self.signs)).replace("'","")
        ret += ", " + self.units
        return ret

    # Converts electron wavenumber to kinetic energy
    def ke(self, k):
        ene = (1.0/self.get_ene_conv())*k**2
        return nw.complex(ene)

    # Converts election wavenumber to total energy
    def e(self, ch, k):
        return ke(k) + self.thresholds[ch]

    # Converts free electron energy to wavenumber
    def fk(self, ene): #free k
        k = nw.sqrt(self.get_ene_conv()*ene)
        return nw.complex(k)

    # Converts channel electron energy to wavenumber
    def k(self, ch, ene):
        if self.sign_sel == signs_pos:
            return self._kpos(ch, ene)
        elif self.sign_sel == signs_specified:
            return self._k_specified(ch, ene)
        elif self.sign_sel == signs_bndandres:
            return self._k_bnd_and_res(ch, ene)
        elif self.sign_sel == signs_ana_over_axis:
            return self._k_ana_over_axis(ch, ene)
        elif self.sign_sel == signs_ana_over_thres:
            return self._k_ana_over_thres(ch, ene)

    # Returns threshold potential for a channel
    def th(self, ch):
        return self.thresholds[ch]
        
    # Returns angular momentum for a channel
    def l(self, ch):
        return self.ls[ch]

    def get_number_channels(self):
        return len(self.thresholds)

    def get_units(self):
        return self.units

    def get_ene_conv(self):
        return self.ene_conv

    def is_elastic(self):
        return self.thresholds[1:] == self.thresholds[:-1]

    def _kpos(self, ch, ene):
        return nw.sqrt(self._get_value(ch, ene))

    def _get_value(self, ch, ene):
        #if ene.real < self.thresholds[ch]:
        #    print "WARNING!"
        #    print str(ene) + "   " + str(self.thresholds[ch])
        return self.get_ene_conv()*(ene - self.thresholds[ch])

    def _k_specified(self, ch, ene):
        mult = 1.0
        if self.signs is not None:
            mult = self.signs[ch]
        return mult * self._kpos(ch, ene)

    def _k_bnd_and_res(self, ch, ene):
        k = self._kpos(ch, ene)
        if ene.real <= self.thresholds[ch]: 
            # We want to be on the physical here
            if ene.imag >= 0.0:
                sign = 1.0
            else:
                sign = -1.0
        elif ene.real > self.thresholds[ch]: 
            # We want to be on the unphysical here
            if ene.imag >= 0.0:
                sign = -1.0
            else:
                sign = 1.0
        return sign*k

    # This keeps function analytical AS LONG as you don't cross thresholds.
    def _k_ana_over_axis(self, ch, ene):
        k = self._kpos(ch, ene)
        if ene.real <= self.thresholds[ch]: 
            #We want smooth transition over the real axis
            if ene.imag >= 0.0:
                sign = 1.0
            else:
                sign = -1.0
        elif ene.real > self.thresholds[ch]: 
            #We want smooth transition over the real axis
            sign = 1.0
        return sign*k

    # This keeps function analytical AS LONG as you don't cross axis.
    def _k_ana_over_thres(self, ch, ene):
        k = self._kpos(ch, ene)
        if ene.imag >= 0.0:
            sign = -1.0
        else:
            sign = 1.0
        return sign*k

def use_python_types(dps=nw.dps_default_python):
    nw.use_python_types(dps)

def use_mpmath_types(dps=nw.dps_default_mpmath):
    nw.use_mpmath_types(dps)
