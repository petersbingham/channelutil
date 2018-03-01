import pynumwrap as nw

REDUCED_MASS = 1.0  
   
SIGNS_POS = "SignsPositive"
SIGNS_SPECIFIED = "SignsSpecified"

SIGNS_BNDANDRES = "SignsForBndAndRes"

SIGNS_ANA_OVER_AXIS = "SignsAnaOverAxis"
SIGNS_ANA_OVER_THRES = "SignsAnaOverThres"

MASSMULT_RYDBERGS = 1.0
MASSMULT_HARTREES = 2.0

class calculator:
    def __init__(self, thresholds=None, ls=None, ktype=SIGNS_POS, ksigns=None, massMult=MASSMULT_RYDBERGS):
        if thresholds is None:
            self.thresholds = []
        else:
            self.thresholds = thresholds
        if ls is None:
            self.ls = [0]*len(thresholds)
        else:
            self.ls = ls
        self.massMult = massMult
        self.ktype = ktype
        self.ksigns = ksigns

    def __str__(self):
        if self.ksigns is None:
            return self.ktype
        else:
            return nw.floatList(self.ksigns)

    # Converts free electron wavenumber to energy
    def fe(self, k):
        ene = (1.0/self.getMult())*k**2
        return nw.complex(ene)

    # Converts free electron energy to wavenumber
    def fk(self, ene): #free k
        k = nw.sqrt(self.getMult()*ene)
        return nw.complex(k)

    # Converts channel electron energy to wavenumber
    def k(self, ch, ene):
        if self.ktype == SIGNS_POS:
            return self._kpos(ch, ene)
        elif self.ktype == SIGNS_SPECIFIED:
            return self._kspecified(ch, ene)
        elif self.ktype == SIGNS_BNDANDRES:
            return self._kbndAndRes(ch, ene)
        elif self.ktype == SIGNS_ANA_OVER_AXIS:
            return self._kanaOverAxis(ch, ene)
        elif self.ktype == SIGNS_ANA_OVER_THRES:
            return self._kanaOverThres(ch, ene)

    # Returns angular momentum for a channel
    def l(self, ch):
        return self.ls[ch]

    def getMult(self):
        return self.massMult*REDUCED_MASS

    def isElastic(self):
        return self.thresholds[1:] == self.thresholds[:-1]

    def _kpos(self, ch, ene):
        return nw.sqrt(self._getValue(ch, ene))

    def _getValue(self, ch, ene):
        #if ene.real < self.thresholds[ch]:
        #    print "WARNING!"
        #    print str(ene) + "   " + str(self.thresholds[ch])
        return self.getMult()*(ene - self.thresholds[ch])

    def _kspecified(self, ch, ene):
        mult = 1.0
        if self.ksigns is not None:
            mult = self.ksigns[ch]
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
