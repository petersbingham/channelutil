import os
import sys
basedir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,basedir+'/../..')

import channelutil as cu

import unittest

class test_channelutil(unittest.TestCase):
    def runTest(self):
        self.assertEqual(True,True)

if __name__ == "__main__":
    #Just for debug
    b = test_channelutil()
    b.runTest()
