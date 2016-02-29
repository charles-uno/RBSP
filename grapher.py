#!/usr/bin/env python

# Charles McEachern

# Spring 2016

# #############################################################################
# #################################################################### Synopsis
# #############################################################################

# This routine reads in pickles created by grabber.py, rotates the data into
# the coordinates we want, and plots it. 

# #############################################################################
# ######################################################### Load Python Modules
# #############################################################################

try:
  import cPickle as pickle
except ImportError:
  import pickle
import matplotlib.pyplot as plt
import numpy as np
import os
from plotmod import *
from random import choice

# #############################################################################
# ######################################################################## Main
# #############################################################################

def main():

  # Location of the output directories. 
  outdir = '/media/My Passport/RBSP/pickles/'

  # Location where images should be saved. 
  savedir = '/home/user1/mceachern/Desktop/plots/' + now() + '/'
  if '-i' in argv and not os.path.exists(savedir):
    os.mkdir(savedir)

  # Each event has its own directory full of pickles. Go through all of them
  # if we're making images. For debugging, just look at a random one. 
  if '-i' in argv:
    pkldirs = sorted( os.listdir(outdir) )
  else:
    pkldirs = [ choice( os.listdir(outdir) ) ]
#    pkldirs = ['a_20121023_221000']

  for pkldir in pkldirs:

    # From a directory, construct an event object. 
    ev = event(outdir + pkldir + '/')

    # Make sure the event actually holds data before we try to access it. 
    if not ev.isok():
      print 'SKIPPING ' + pkldir + ' DUE TO BAD DATA. '
      continue

    # Initialize plot window object. Four double-wide rows. Label them. 
    PW = plotWindow(nrows=3, ncols=-2)
    ylabels = ( notex('Poloidal (\\frac{mV}{m}) (nT)'), 
#                notex('Poloidal FFT (\\frac{mV}{m}) (nT)'), 
                notex('Toroidal (\\frac{mV}{m}) (nT)'), 
                notex('Field-Aligned (\\frac{mV}{m}) (nT)') )
    [ PW[i].setParams(ylabel=ylbl) for i, ylbl in enumerate(ylabels) ]

    # Assemble plot title and labels based on event parameters. 
    title = notex('In Situ Electric (Blue) and Magnetic (Red) Fields')
    xlabel = notex( 'Time (hh:mm) on ' + calendar(ev.date) )
    collabel = ( notex( 'RBSP--' + ev.probe.upper() ) + ' \\qquad L \\!=\\! ' +
                 format(np.average( ev.get('lshell') ), '.1f') + ' \\qquad ' +
                 notex( clock( 3600*np.average( ev.get('MLT') ) ) +
                 ' MLT') + ' \\qquad ' +
                 format(np.average( ev.get('mlat') ), '.0f') + '^\\circ' +
                 notex(' Magnetic Latitude') )
    PW.setParams(title=title, collabels=[collabel], xlabel=xlabel)

    # Manually set the ticks and tick labels on the x axis. 
    xtks, xtls = range(11), ['']*11
    for i in range(1, 11, 2):
      xtls[i] = '$' + notex( clock(ev.get('t')[0] + 60*i) ) + '$'
    PW.setParams(x=ev.get('tcoord'), xticks=xtks, xticklabels=xtls)

    # Plot poloidal, toroidal, and field-aligned components. 
    [ PW[0].setLine(ev.get(v), c) for v, c in ( ('Ey', 'b'), ('Bx', 'r') ) ]
    [ PW[1].setLine(ev.get(v), c) for v, c in ( ('Ex', 'b'), ('By', 'r') ) ]
    [ PW[2].setLine(ev.get(v), c) for v, c in ( ('Ez', 'b'), ('Bz', 'r') ) ]

    # Fourier transform the poloidal components, to eyeball phase offset. 

    # Harmonic at data resolution. 
    t = ev.get('tcoord')
    def harm(n):
      return np.exp( 2j*np.pi*n*t / ( t[-1] - t[0] ) )

    # Harmonic at fine resolution, for a smooth plot. 
    nfine = 1000
    tfine = np.linspace(t[0], t[-1], nfine)
    def harmfine(n):
      return np.exp( 2j*np.pi*n*tfine/ ( tfine[-1] - tfine[0] ) )

    # The domain is finite, so we use a discrete set of Fourier weights. 
    nmax = 20
    wB, wE = np.zeros(nmax, dtype=np.complex), np.zeros(nmax, dtype=np.complex)

    # Compute Fourier series weights. 
    dt = t[1] - t[0]
    for n in range(nmax):
      wB[n] = np.sum(ev.get('Bx')*harm(-n)*dt)/np.sum(harm(n)*harm(-n)*dt)
      wE[n] = np.sum(ev.get('Ey')*harm(-n)*dt)/np.sum(harm(n)*harm(-n)*dt)

    # Reconstitute a waveform. 
    Bf = np.zeros(nfine, dtype=np.complex)
    Ef = np.zeros(nfine, dtype=np.complex)
    for n in range(nmax):
      Bf = Bf + harmfine(n)*wB[n]
      Ef = Ef + harmfine(n)*wE[n]

    # Compute the complex phase between the electric and magnetic fields. 
    Bangle, Eangle = np.angle(Bf), np.angle(Ef)
    phase = Eangle - Bangle
    # Make sure we stay between -180 and 180 degrees. 
    phase = np.where(phase<-np.pi, phase + 2*np.pi, phase)
    phase = np.where(phase>np.pi, phase - 2*np.pi, phase)

    # Plot the phase, in radians. It's weird to put it on the same scale as the
    # poloidal fields, but that's fine for the moment. 
    PW[0].setLine(tfine, phase, 'g')

    # Save the plot as an image. 
    if '-i' in argv:
      PW.render(savedir + pkldir + '.png')
    else:
      PW.render()

  return


# #############################################################################
# ################################################################ Event Object
# #############################################################################

# Given the name of a directory holding an event, this event grabs the data out
# of that directory and assembles it into a useful form. 
class event:

  # ---------------------------------------------------------------------------
  # --------------------------------------------------- Initialize Event Object
  # ---------------------------------------------------------------------------

  def __init__(self, datadir):

    # From the name of the directory, get the probe name and the timestamp. 
    self.probe, self.date, self.t0, self.t1 = self.parse(datadir)

    # Load the pickles for time, fields, and position. 
    for var in ('time', 'bgse', 'egse', 'xgse', 'lshell', 'mlt', 'mlat'):
      self.__dict__[var] = loadpickle(datadir + var + '.pkl')

    # Shift the time coordinate to be zero at midnight, rather than in 1970. 
    self.t = self.time - self.time[0]

    # Figure out the indeces that correspond to the event start and end. 
    self.i0, self.i1 = np.argmax(self.t>self.t0), np.argmax(self.t>self.t1)

    # Compute the background magnetic field using a rolling average. 
    self.b0gse = self.getbg(self.t, self.bgse)

    # Get the parallel, azimuthal, and crosswise unit vectors. 
    self.xhat, self.yhat, self.zhat = self.uvecs(self.xgse, self.b0gse)

    # Rotate the magnetic field into dipole coordinates. 
    self.bx, self.by, self.bz = self.rotate(self.bgse - self.b0gse)

    # Do the same for the magnetic fields. 
    self.ex, self.ey, self.ez = self.rotate(self.egse)

    return

  # ---------------------------------------------------------------------------
  # ------------------------------------------------------ Parse Directory Name
  # ---------------------------------------------------------------------------

  # Split the directory into the probe name, the date stamp, and the time. Note
  # that all events are ten minutes long by construction. 
  def parse(self, pkldir):
    probe, date, time = pkldir.rstrip('/').split('/')[-1].split('_')
    hh, mm, ss = int( time[0:2] ), int( time[2:4] ), int( time[4:6] )
    t0 = ss + 60*mm + 3600*hh
    return probe, date, t0, t0 + 600

  # ---------------------------------------------------------------------------
  # ----------------------------------------- Compute Background Magnetic Field
  # ---------------------------------------------------------------------------

  # Take a rolling average of the magnetic field to estimate the background
  # field. The average gets truncated at the edges -- notably, none of Lei's
  # events happen within ten minutes of midnight. 
  def getbg(self, t, B):
    B0 = np.empty(B.shape)
    for i in range( B.shape[1] ):
      # Find indeces five minutes in the past and five minutes in the future. 
      ipas, ifut = np.argmax(t > t[i] - 300), np.argmax(t > t[i] + 300)
      # Correct the upper bound in casae it overflows the array domain. 
      ifut = ifut if ifut>0 else B.shape[1]
      # Take the average. 
      B0[:, i] = np.average(B[:, ipas:ifut], axis=1)
    return B0

  # ---------------------------------------------------------------------------
  # ----------------------------------------------------- Compute Basis Vectors
  # ---------------------------------------------------------------------------

  # Compute the dipole coordinate directions. The zhat unit vector lines up
  # with the background magnetic field, yhat is azimuthally eastward, and xhat
  # completes the orthonormal coordinate system. 
  def uvecs(self, X, B0):
    zhat = unit(B0)
    yhat = unit( np.cross(zhat, X, axis=0) )
    return np.cross(yhat, zhat, axis=0), yhat, zhat

  # ---------------------------------------------------------------------------
  # ------------------------------------- Rotate from GSE to Dipole Coordinates
  # ---------------------------------------------------------------------------

  # Dot vectors in GSE coordinates with the basis vectors (also in GSE
  # coordinates) to end up with vector components in the dipole basis. 
  def rotate(self, vgse):
    return dot(vgse, self.xhat), dot(vgse, self.yhat), dot(vgse, self.zhat)

  # ---------------------------------------------------------------------------
  # --------------------------------------------------------- Access Event Data
  # ---------------------------------------------------------------------------

  # Access a quantity. Only the slice during the ten-minute event is returned. 
  def get(self, var):
    if var=='tcoord':
      return ( self.get('t') - self.get('t')[0] )/60
    else:
      return self.__dict__[ var.lower() ][..., self.i0:self.i1]

  # ---------------------------------------------------------------------------
  # -------------------------------------------------------- Check for Bad Data
  # ---------------------------------------------------------------------------

  # Data can be bad for a number of reasons. Sometimes a chunk of the day is
  # missing due to eclipse. Sometimes the spin axis is too close to a
  # measurement axis, so E dot B doesn't give good information. In any case, 
  # the bad data is pretty easy to identify, so we just skip it. 
  def isok(self):
    if self.get('t').size==0:
      return False
    return all( all( np.isfinite( self.get(var) ) ) for var in ('Ex', 'Ey') )

'''
# Make sure the given variable is a 3-by-N array.  
def is3d(x):
  return isinstance(x, np.ndarray) and len(x.shape)==2 and x[0]==3
'''

# #############################################################################
# ############################################################ Helper Functions
# #############################################################################

# Put the slashes back into a date string. 
def calendar(date):
  return date[0:4] + '--' + date[4:6] + '--' + date[6:8]

# Compute clock time from a count of seconds from midnight. 
def clock(sfm, seconds=False):
  if not np.isfinite(sfm):
    return '??:??'
  hh = sfm/3600
  if seconds:
    mm = (sfm%3600)/60
    ss = sfm%60
    return znt(hh, 2) + ':' + znt(mm, 2) + ':' + znt(ss, 2)
  else:
    mm = int( format( (sfm%3600)/60., '.0f') )
    return znt(hh, 2) + ':' + znt(mm, 2)

# Load all of the pickles in a directory into a dictionary. 
def load(datadir):
  data = {}
  for pklname in os.listdir(datadir):
    with open(datadir + pklname, 'rb') as handle:
      data[ pklname[:-4] ] = pickle.load(handle)
  return data

# Load a pickle file. 
def loadpickle(pklpath):
  with open(pklpath, 'rb') as handle:
    return pickle.load(handle)

# Based on the pickle directory name, return the probe name and the start and
# end time of the event, in seconds from midnight. Note that all events are ten
# minutes long by construction. 
def pdtt(pkldir):
  probe, date, time = pkldir.rstrip('/').split('/')[-1].split('_')
  hh, mm, ss = int( time[0:2] ), int( time[2:4] ), int( time[4:6] )
  t0 = ss + 60*mm + 3600*hh
  return probe, date, t0, t0 + 600

# Take a rolling average of the magnetic field to estimate the background
# field. The average gets truncated at the edges -- notably, none of Lei's
# events happen within ten minutes of midnight. 
def getbg(t, B, tavg=600):
  B0 = np.empty(B.shape)
  for i in range( B.shape[1] ):
    # Find indeces five minutes in the past and five minutes in the future. 
    ipas = np.argmax( t > t[i] - tavg/2 )
    ifut = np.argmax( t > t[i] + tavg/2 )
    # If we're looking for something past the edge of the array, argmax will
    # return 0. That a problem for the upper bound. 
    ifut = ifut if ifut>0 else B.shape[1]
    # Take the average. 
    B0[:, i] = np.average(B[:, ipas:ifut], axis=1)
  return B0

# Dot product of a pair of vectors or a pair of arrays of vectors. 
def dot(v0, v1, axis=0):
  return np.sum(v0*v1, axis=axis)

# Scale a vector, or an array of vectors, to unit vectors. 
def unit(v):
  return v/np.sqrt( dot(v, v) )

# #############################################################################
# ########################################################### For Importability
# #############################################################################

if __name__=='__main__':
  main()


