from __future__ import division
import time
import sys

toolbar_width = 40

class ProgressBar:
    def __init__(self, niter=40, ndash=40):
        self.ndash = ndash
        self.niter = niter
        self.dash  = 0
        self.iter  = 0

        # setup toolbar
        sys.stdout.write("[%s]" % (" " * self.ndash))
        sys.stdout.flush()
        # return to start of line, after '['
        sys.stdout.write("\b" * (self.ndash+1)) 

    def iterate(self):
        if self.iter > self.niter:
            print("Progressbar is full")
            return

        self.iter = self.iter + 1

        toadd = self.iter / self.niter * self.ndash - self.dash
        toadd = int(toadd)
        if toadd > 0:
            for i in xrange(toadd):
                sys.stdout.write("-")
                sys.stdout.flush()
                self.dash = self.dash + 1
