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

    def e(self, k, primType=False):
        ene = (1.0/self.getMult())*k**2
        if primType:
            return complex(ene)
        else:
            return nw.complex(ene)

    def fk(self, ene, primType=False): #free k
        k = nw.sqrt(self.getMult()*ene)
        if primType:
            return complex(k)
        else:
            return nw.complex(k)

    def kl(self, ch, ene, add=0.0):
        k = self.k(ch, ene)
        return nw.pow(k, self.ls[ch]+add)

    def k(self, ch, ene):
        if self.ktype == SIGNS_POS:
            return self.kpos(ch, ene)
        elif self.ktype == SIGNS_SPECIFIED:
            return self.kspecified(ch, ene)
        elif self.ktype == SIGNS_BNDANDRES:
            return self.kbndAndRes(ch, ene)
        elif self.ktype == SIGNS_ANA_OVER_AXIS:
            return self.kanaOverAxis(ch, ene)
        elif self.ktype == SIGNS_ANA_OVER_THRES:
            return self.kanaOverThres(ch, ene)

    def l(self, ch):
        return self.ls[ch]

    def kpos(self, ch, ene):
        return nw.sqrt(self._getValue(ch, ene))

    def kspecified(self, ch, ene):
        mult = 1.0
        if self.ksigns is not None:
            mult = self.ksigns[ch]
        return mult * self.kpos(ch, ene)

    def kbndAndRes(self, ch, ene):
        k = self.kpos(ch, ene)
        if ene.real <= self.thresholds[ch]: #We want to be on the physical here
            if ene.imag >= 0.0:
                sign = 1.0
            else:
                sign = -1.0
        elif ene.real > self.thresholds[ch]: #We want to be on the unphysical here
            if ene.imag >= 0.0:
                sign = -1.0
            else:
                sign = 1.0
        return sign*k

    def kanaOverAxis(self, ch, ene):
        #This keeps function analytical AS LONG as you don't cross thresholds.
        k = self.kpos(ch, ene)
        if ene.real <= self.thresholds[ch]: #We want smooth transition over the real axis
            if ene.imag >= 0.0:
                sign = 1.0
            else:
                sign = -1.0
        elif ene.real > self.thresholds[ch]: #We want smooth transition over the real axis
            sign = 1.0
        return sign*k

    def kanaOverThres(self, ch, ene):
        #This keeps function analytical AS LONG as you don't cross axis.
        k = self.kpos(ch, ene)
        if ene.imag >= 0.0:
            sign = -1.0
        else:
            sign = 1.0
        return sign*k

    def getPhase(self, ch, ene):
        if ene.real <= self.thresholds[ch]:
            return 0.0
        else:
            return nw.pi

    def getMult(self):
        return self.massMult*REDUCED_MASS

    def isElastic(self):
        return self.thresholds[1:] == self.thresholds[:-1]

    def _getValue(self, ch, ene):
        #if ene.real < self.thresholds[ch]:
        #    print "WARNING!"
        #    print str(ene) + "   " + str(self.thresholds[ch])
        return self.getMult()*(ene - self.thresholds[ch])
