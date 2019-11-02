import re
import numpy as np
import os
import gzip
import scipy.linalg
import scipy.io
import time, sys
#import psutil
#from guppy import hpy

class Traj(object):

    def __init__(self, path, title):
        self.path = path
        self.natoms = 0
        self.mass = 28
        self.is_read = False
        self.title = title

    def read_traj_amber(self, crd = '' ):
        if self.is_read:
            print
            'Traj already read once'
        else:
            print "Start to read trajectory"
            if crd == '':
                crd = '%s.crd' % self.title
            f = scipy.io.netcdf.netcdf_file('%s/%s' % (self.path, crd) )
            print f.dimensions['atom']
