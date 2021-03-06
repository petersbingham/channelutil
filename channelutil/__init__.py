import pynumwrap as nw
from units import *
from channelutil.release import __version__

signs_pos = "SignsPositive"
signs_specified = "SignsSpecified"

signs_bndandres = "SignsForBndAndRes"

signs_ana_over_axis = "SignsAnaOverAxis"
signs_ana_over_thres = "SignsAnaOverThres"

# Performs asymptotic calculations
class AsymCalc:
    def __init__(self, units, angmoms=None, tot_spin=None, targ_spins=None,
                 thresholds=None, sign_sel=signs_pos, signs=None):
        if angmoms is None:
            if thresholds is not None:
                self.angmoms = [0]*len(thresholds)
            else:
                self.angmoms = [0]
        else:
            self.angmoms = angmoms

        if tot_spin is None:
            self.tot_spin_ = 0.5

        if targ_spins is None:
            self.targ_spins_ = 0.

        if thresholds is None:
            self.thresholds = [0.]*len(self.angmoms)
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
        ret = "angmoms:" + str(self.angmoms) + ", tot_spin:" + str(self.tot_spin_)
        ret += ", targ_spins:" + str(self.targ_spins_)
        ret += ", thres:" + str(self.thresholds)
        ret += ", " + self.units
        return ret

    def sign_str(self):
        if self.signs is None:
            return self.sign_sel
        else:
            return str(map(lambda x:str(x),self.signs)).replace("'","")

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
    def thres(self, ch):
        return self.thresholds[ch]
        
    # Returns angular momentum for a channel
    def angmom(self, ch):
        return self.angmoms[ch]
    
    # Returns total spin for system
    def tot_spin(self):
        return self.tot_spin_

    # Returns spin for a channel
    def targ_spins(self, ch):
        try:
            return self.targ_spins_[ch]
        except TypeError:
            return self.targ_spins_

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

def use_python_types():
    nw.use_python_types()

def use_mpmath_types(dps=nw.dps_default_mpmath):
    nw.use_mpmath_types(dps)
