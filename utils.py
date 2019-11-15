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

    def read_traj_amber(self, crd = '', top = '' ):
        if self.is_read:
            print
            'Traj already read once'
        else:
            print "Start to read trajectory"
            if top == '':
                top = '%s.top' % self.title
            with open('%s/%s' % (self.path, top)) as f:
                lines = f.readlines()
                do = False
                self.atoms = []
                for i,line in enumerate(lines):
                    if 'ATOM_NAME' in line:
                        do = True
                    elif 'FLAG CHARGE' in line:
                        do = False
                    else:
                        if do and 'FORMAT' not in line:
                            self.atoms += line.split()
            self.set_atoms = set(self.atoms)
            if crd == '':
                crd = '%s.crd' % self.title
            f = scipy.io.netcdf.netcdf_file('%s/%s' % (self.path, crd) )
            self.natoms = f.dimensions['atom']
            if self.natoms != len(self.atoms):
                print "Error in number of atoms"
                print "self.natoms", self.natoms
                print "len(self.atoms)", len(self.atoms)
            else:
                print 'Number of atoms: %s' % self.natoms
            self.steps = f.dimensions['frame']
            if self.steps is None:
                self.steps = f.variables['time'].shape[0]
            print 'Steps: %s' % self.steps
            self.positions = np.zeros((self.steps, self.natoms, 3))
            self.velocities = np.zeros((self.steps, self.natoms, 3))
            self.forces = np.zeros((self.steps, self.natoms, 3))
            for step in range(self.steps):
                self.positions[step] = f.variables['coordinates'][step]

    def determine_z_density(self, dist, binwidth):
        dens_zs = {}
        nbins = int(dist/binwidth)
        self.bins_z = binwidth*np.array(range(nbins))
        for atom in self.set_atoms:
            dens_zs[atom] = np.zeros(nbins)
        for iatom in range(self.natoms):
            pos = self.positions[:,iatom, 2]
            pos = pos - dist * np.rint(pos/dist)
            pos = np.rint(pos / binwidth).astype(int)
            for pp in pos:
                dens_zs[self.atoms[iatom]][pp] += 1
        self.dens_zs = dens_zs

