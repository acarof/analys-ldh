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
        self.is_unitcell = False
        self.title = title
        self.bins = {}

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

    def wrap_unitcell(self, unitcell):
        shift = np.rint( self.positions.dot(unitcell) / np.sum(np.power(unitcell, 2), 0))
        self.positions_unitcell = self.positions - shift.dot(unitcell)

    def determine_xyz_density(self, unitcell, binwidth):
        self.wrap_unitcell(unitcell)
        self.dens = {}
        self.dens['z'] = {}
        self.dens['xy'] = {}
        distz = np.dot(unitcell[2], np.array([0,0,1]))
        nbins = int(distz/binwidth)
        self.bins['z'] = binwidth*np.array(range(nbins))
        for atom in self.set_atoms:
            self.dens['z'][atom] = np.zeros(nbins)
        vect = np.diag(unitcell[0:2,0:2])
        print vect
        nbins = [int(x/binwidth) for x in vect]
        self.bins['xy'] = np.meshgrid(*(binwidth*np.array(range(n)) for n in nbins))
        for atom in self.set_atoms:
            self.dens['xy'][atom] = np.zeros([n for n in nbins])
        for iatom in range(self.natoms):
            pos = self.positions_unitcell[:,iatom,:]
            pos = np.rint(pos/binwidth).astype(int)
            for pp in pos:
                self.dens['xy'][self.atoms[iatom]][pp[0], pp[1]]  += 1
                self.dens['z'][self.atoms[iatom]][pp[2]] += 1


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


    def determine_xy_density(self, vect, binwidth):
        vect = np.array(vect)
        dens_xys = {}
        nbins = [int(x/binwidth) for x in vect]
        self.bins_xy = np.meshgrid(*(binwidth*np.array(range(n)) for n in nbins))
        for atom in self.set_atoms:
            dens_xys[atom] = np.zeros([n for n in nbins])
        for iatom in range(self.natoms):
            pos = self.positions[:,iatom, 0:2]
            pos = pos - vect * np.rint(pos/vect)
            #print pos[0]
            #print pos[0]/binwidth
            #print np.rint(pos[0]/binwidth)
            pos = np.rint(pos / binwidth).astype(int)
            #print pos[0]
            for pp in pos:
                #print 'pp',pp
                dens_xys[self.atoms[iatom]][pp[0], pp[1]] += 1
                #print dens_xys[self.atoms[iatom]]
                #stop
        self.dens_xys = dens_xys

def write_map(path, map_):
    with open(path, 'w') as f:
        for line in map_:
            f.write('%s\n' % ' '.join(map(str, line)))

def write_function_z(path, bins, values):
    with open(path, 'w') as f:
         f.write('Z   F\n')
         for r, value in zip(bins, values):
             f.write('%s    %s\n' % (r, value))
