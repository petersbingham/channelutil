import pynumwrap as nw
from units import *

signs_pos = "SignsPositive"
signs_specified = "SignsSpecified"

signs_bndandres = "SignsForBndAndRes"

signs_ana_over_axis = "SignsAnaOverAxis"
signs_ana_over_thres = "SignsAnaOverThres"

# Performs asymptotic calculations
class asymCal:
    def __init__(self, units, ls=None, isPolar=False, thresholds=None, 
                 signSel=signs_pos, signs=None):
        if ls is None:
            if thresholds is not None:
                self.ls = [0]*len(thresholds)
            else:
                self.ls = [0]
        else:
            self.ls = ls

        self.polar = isPolar

        if thresholds is None:
            self.thresholds = [0.]*len(self.ls)
        else:
            self.thresholds = thresholds

        self.units = units
        if units == RYDs:
            self.eneConv = 1.0
        elif units == HARTs:
            self.eneConv = 1./RYD_to_HART
        elif units == HARTs:
            self.eneConv = 1./RYD_to_EV
        self.signSel = signSel
        self.signs = signs

    def __str__(self):
        if self.signs is None:
            return self.signSel
        else:
            return nw.floatList(self.signs)

    # Converts electron wavenumber to kinetic energy
    def ke(self, k):
        ene = (1.0/self.getEneConv())*k**2
        return nw.complex(ene)

    # Converts election wavenumber to total energy
    def e(self, ch, k):
        return ke(k) + self.thresholds[ch]

    # Converts free electron energy to wavenumber
    def fk(self, ene): #free k
        k = nw.sqrt(self.getEneConv()*ene)
        return nw.complex(k)

    # Converts channel electron energy to wavenumber
    def k(self, ch, ene):
        if self.signSel == signs_pos:
            return self._kpos(ch, ene)
        elif self.signSel == signs_specified:
            return self._kspecified(ch, ene)
        elif self.signSel == signs_bndandres:
            return self._kbndAndRes(ch, ene)
        elif self.signSel == signs_ana_over_axis:
            return self._kanaOverAxis(ch, ene)
        elif self.signSel == signs_ana_over_thres:
            return self._kanaOverThres(ch, ene)

    # Returns threshold potential for a channel
    def th(self, ch):
        return self.thresholds[ch]
        
    # Returns angular momentum for a channel
    def l(self, ch):
        return self.ls[ch]

    def getNumberChannels(self):
        return len(self.thresholds)

    def getUnits(self):
        return self.units

    def getEneConv(self):
        return self.eneConv

    def isElastic(self):
        return self.thresholds[1:] == self.thresholds[:-1]

    def isPolar(self):
        return self.polar

    def _kpos(self, ch, ene):
        return nw.sqrt(self._getValue(ch, ene))

    def _getValue(self, ch, ene):
        #if ene.real < self.thresholds[ch]:
        #    print "WARNING!"
        #    print str(ene) + "   " + str(self.thresholds[ch])
        return self.getEneConv()*(ene - self.thresholds[ch])

    def _kspecified(self, ch, ene):
        mult = 1.0
        if self.signs is not None:
            mult = self.signs[ch]
        return mult * self._kpos(ch, ene)

    def _kbndAndRes(self, ch, ene):
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
    def _kanaOverAxis(self, ch, ene):
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
    def _kanaOverThres(self, ch, ene):
        k = self._kpos(ch, ene)
        if ene.imag >= 0.0:
            sign = -1.0
        else:
            sign = 1.0
        return sign*k

def usePythonTypes(dps=nw.dps_default_python):
    nw.usePythonTypes(dps)

def useMpmathTypes(dps=nw.dps_default_mpmath):
    nw.useMpmathTypes(dps)
