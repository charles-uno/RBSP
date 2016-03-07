#!/usr/bin/env python

# Charles McEachern

# Spring 2016

# Note: This document wraps at column 80. 

# #############################################################################
# #################################################################### Synopsis
# #############################################################################

# The event class reads in a bunch of pickles of RBSP data, and transforms them
# into a dipole-aligned coordinate system for easy plotting. It can also give
# spectral properties of the data, such as the coherence between poloidal
# electric and magnetic fields, or their cross-spectral density. Finally, the
# event class provides support for plotting -- it gives axis bounds, labels, 
# etc, nicely formatted for LaTeX. 

# #############################################################################
# ####################################################### Import Python Modules
# #############################################################################

from calendar import timegm
import numpy as np
from numpy import pi
import os
try:
  import cPickle as pickle
except ImportError:
  import pickle
from plotmod import *
from scipy import signal
from time import gmtime

# #############################################################################
# ################################################################ Event Object
# #############################################################################

# Given the name of a directory holding an event, this event grabs the data out
# of that directory and assembles it into a useful form. 
class event:

  # ---------------------------------------------------------------------------
  # --------------------------------------------------- Initialize Event Object
  # ---------------------------------------------------------------------------

  def __init__(self, probe, date, time='00:00:00', mins=10):

    # Where is the plasmapause at the time of this event? 
    self.lpp = lpp(date=date, time=time)

    # Store event identifiers, and index this event with a unique name. 
    self.probe, self.date = probe, date.replace('-', '')
    self.name = self.date + '_' + time.replace(':', '') + '_' + probe

    # Load the pickles for time, fields, and position. 
    datadir = '/media/My Passport/rbsp/pkls/' + self.date + '/' + probe + '/'
    for var in ('time', 'bgse', 'egse', 'xgse', 'lshell', 'mlt', 'mlat'):
      self.__dict__[var] = loadpickle(datadir + var + '.pkl')

    # Shift the time coordinate to be zero at midnight, rather than in 1970. 
    self.t, self.dt = self.time - self.time[0], self.time[1] - self.time[0]

    # Figure out the indeces that correspond to the event start and end. 
    self.t0 = timeint(date=None, time=time)
    self.t1 = self.t0 + 60*mins
    self.i0, self.i1 = np.argmax(self.t>self.t0), np.argmax(self.t>self.t1)

    # Compute the background magnetic field using a rolling average. 
    self.b0gse = self.getbg(self.t, self.bgse)

    # Get the parallel, azimuthal, and crosswise unit vectors. 
    self.xhat, self.yhat, self.zhat = self.uvecs(self.xgse, self.b0gse)

    # Rotate the magnetic field into dipole coordinates. 
    self.bx, self.by, self.bz = self.rotate(self.bgse - self.b0gse)

    # Do the same for the electric fields. 
    self.ex, self.ey, self.ez = self.rotate(self.egse)

    return

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------- Create Adjacent Events
  # ---------------------------------------------------------------------------

  # Note that this is not smart enough to handle day boundaries. 
  def prev(self):
    return event( probe=self.probe, date=self.date, 
                  time=timestr(2*self.t0 - self.t1)[1] )

  # Note that this is not smart enough to handle day boundaries. 
  def next(self):
    return event( probe=self.probe, date=self.date, time=timestr(self.t1)[1] )

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
  # --------------------------------------------------------- Event Data Access
  # ---------------------------------------------------------------------------

  # Access a quantity. Only the slice during the ten-minute event is returned. 
  def get(self, var):
    # Allow requests for poloidal and toroidal fields. 


    # The parallel electric field is negligible near apogee. But we care if Bz
    # is coherent with Ey. So just use Ey for Ez for the moment. 
    pt = {'ep':'ey', 'et':'ex', 'ez':'ey', 'bp':'bx', 'bt':'by'}
#    pt = {'ep':'ey', 'et':'ex', 'bp':'bx', 'bt':'by'}


    if var.lower() in pt:
      return self.get( pt[ var.lower() ] )
    if var=='tcoord':
      return ( self.get('t') - self.get('t')[0] )/60
    elif var=='tfine':
      return np.linspace(self.get('tcoord')[0], self.get('tcoord')[-1], 1000)
    # For fields, subtract off the average value. 
    elif var.lower() in ('bx', 'by', 'bz', 'ex', 'ey', 'ez'):
      arr = self.__dict__[ var.lower() ][..., self.i0:self.i1]
      return arr - np.mean(arr)
    else:
      return self.__dict__[ var.lower() ][..., self.i0:self.i1]

  # Average over the ten-minute event to get a single value.
  def avg(self, var):
    return np.average( self.get(var) )

  # Nicely-formatted label of parameter values. 
  def lbl(self, var):
    if var.lower()=='lshell':
      return format(self.avg('lshell'), '.1f')
    elif var.lower()=='mlt':
      return notex( timestr( 3600*self.avg('mlt') )[1][:5] )
    elif var.lower()=='mlat':
      return notex(format(self.avg('mlat'), '+.0f') + '^\\circ')
    else:
      return '???'

  # ---------------------------------------------------------------------------
  # ----------------------------------------------------------- Data Validation
  # ---------------------------------------------------------------------------

  # Data can be bad for a number of reasons. Sometimes a chunk of the day is
  # missing due to eclipse. Sometimes the spin axis is too close to a
  # measurement axis, so E dot B doesn't give good information. In any case, 
  # the bad data is pretty easy to identify, so we just skip it. 
  def isok(self):
    if self.get('t').size==0:
      return False
    return all( all( np.isfinite( self.get(var) ) ) for var in ('Ex', 'Ey') )

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------- Spectral Density
  # ---------------------------------------------------------------------------

  # Frequencies used by the coherence, etc. Does not depend on mode. Use mHz. 
  def frq(self, mode='p'):
    t, B, E = self.get('t'), self.get('B' + mode), self.get('E' + mode)
    return signal.coherence(B, E, fs=1/( t[1] - t[0] ), nperseg=t.size/2, 
                            noverlap=t.size/2 - 1)[0]*1e3



  # Fourier transform of a field component. Use the frequencies employed by the
  # coherence and spectral density routines. Note that this breaks from the
  # Numpy FFT normalization convention -- we introduce the factor of 1/N here,
  # so that the FFT weights can be combined directly with basis functions to
  # recover the original waveform. 
  def fft(self, var):
    t, V = self.get('t'), self.get(var)
    # Refall self.frq() gives frequencies in mHz. 
    temp = [ np.sum( V*np.exp(2j*pi*f*t) ) for f in 1e-3*self.frq() ]
    return np.array(temp)/len(temp)





  # Auto-spectral density. Not renormalized. Note that the fields are real, so
  # the autospectral density also must be real. 
  def asd(self, var):
    t, V = self.get('t'), self.get(var)
    return np.real( signal.csd(V, V, fs=1/( t[1] - t[0] ), 
                      nperseg=t.size/2, noverlap=t.size/2 - 1)[1] )





  # Mode coherence -- at which frequencies do the electric and magnetic field
  # line up for a given mode? Essentially, this is the normalized magnitude of
  # the cross-spectral density. 
  def coh(self, mode):
    t, B, E = self.get('t'), self.get('B' + mode), self.get('E' + mode)
    return signal.coherence(B, E, fs=1/( t[1] - t[0] ), nperseg=t.size/2, 
                            noverlap=t.size/2 - 1)[1]

  # Absolute value of cross-spectral density, normalized to the unit interval. 
  # The units are sorta arbitrary, so we're not really losing any information. 
  # Note that the angle can be recovered from the phase lag function below. 
  def csd(self, mode):
    t, B, E = self.get('t'), self.get('B' + mode), self.get('E' + mode)
    temp = signal.csd(B, E, fs=1/( t[1] - t[0] ), nperseg=t.size/2, 
                      noverlap=t.size/2 - 1)[1]
    return np.abs(temp)/np.max( np.abs(temp) )

  # Phase offset -- at each frequency, for each mode, what's the phase lag
  # between of the electric field relative to the magnetic field? This is
  # determined from the complex phase of the cross-spectral density. 
  def lag(self, mode, deg=True):
    t, B, E = self.get('t'), self.get('B' + mode), self.get('E' + mode)
    return np.angle(signal.csd(B, E, fs=1/( t[1] - t[0] ), nperseg=t.size/2, 
                      noverlap=t.size/2 - 1)[1], deg=True)

  # Spectral density ratio of the compressional and poloidal magnetic fields
  # (compared to the poloidal electric field). This is a way to estimate the
  # azimuthal modenumber; if it's large, Bz should decouple from Bx and Ey. 
  def comp(self, mode='p'):
    t = self.get('t')
    Bx, Bz, Ey = self.get('Bx'), self.get('Bz'), self.get('Ey')
    csdx = signal.csd(Bx, Ey, fs=1/( t[1] - t[0] ), nperseg=t.size/2, 
                      noverlap=t.size/2 - 1)[1]
    csdz = signal.csd(Bz, Ey, fs=1/( t[1] - t[0] ), nperseg=t.size/2, 
                      noverlap=t.size/2 - 1)[1]
    return np.abs(csdz/csdx)

  # Estimate the Alfven speed by looking at the ratio of the electric field to
  # the magnetic field for each frequency component. 
  def va(self, mode='p'):
    t, B, E = self.get('t'), self.get('B' + mode), self.get('E' + mode)
    bsd = signal.csd(B, B, fs=1/( t[1] - t[0] ), nperseg=t.size/2, 
                      noverlap=t.size/2 - 1)[1]
    esd = signal.csd(E, E, fs=1/( t[1] - t[0] ), nperseg=t.size/2, 
                      noverlap=t.size/2 - 1)[1]
    return np.sqrt( np.abs(esd/bsd) )

  # ---------------------------------------------------------------------------
  # -------------------------------------------------- Spectral Density Filters
  # ---------------------------------------------------------------------------

  # Array indicating which frequencies are coherent. By default, results are
  # significant to 1.5 standard deviations. 
  def iscoh(self, mode, n=1):
    return self.coh(mode) > 0.5
#    coh = self.coh(mode)
#    return np.array( [ c > np.mean(coh) + n*np.std(coh) for c in coh ] )

  # Array indicating which frequencies have B (respectively, E) fields above 0.25 nT (mV/m). 
  def has(self, var):
    return np.abs( self.fft(var) ) > 0.25

  # Array indicating which frequencies are spectrally dense. By default, 
  # results are significant to 1.5 standard deviations. 
  def iscsd(self, mode, n=1):
    csd = self.csd(mode)
    return np.array( [ c > np.mean(csd) + n*np.std(csd) for c in csd ] )

  # Array of integers guessing a harmonic for each frequency. The phase lag
  # between electric and magnetic fields gives the parity. For odd modes, the
  # first and third harmonic can be distinguished by frequency. If mlat is
  # within a few degrees of the equator, this function returns all zeros. 
  def harm(self, mode):
    return 2*self.iseven(mode) + 1*self.isodd(mode)

  # Array indicating which frequencies have a phase lag consistent with an odd
  # mode -- north of the magnetic equator, the electric field should lag, and
  # south it should lead. 
  def isodd(self, mode):
    lag = self.lag(mode)
    mlat = self.avg('mlat')
    return (np.abs(mlat) > 3)*( np.sign(mlat) != np.sign(lag) )

  # Array of booleans indicating which frequencies have a phase lag consistent
  # with an even harmonic. 
  def iseven(self, mode):
    lag = self.lag(mode)
    mlat = self.avg('mlat')
    return (np.abs(mlat) > 3)*( np.sign(mlat) == np.sign(lag) )

  # Array of booleans indicating if this frequency is in the pc4 band. Doesn't
  # actually depend on the mode. 
  def ispc4(self, mode='p'):
    return np.array( [ 6 < f < 23 for f in self.frq() ] )

  # ---------------------------------------------------------------------------
  # ---------------------------------------- Labels Indicating Event Properties
  # ---------------------------------------------------------------------------

  # Long label, for the side margin. 
  def label(self):
    return ( notex( 'RBSP--' + self.probe.upper() ) + ' \\qquad L \\!=\\! ' +
             self.lbl('lshell') + ' \\qquad ' + self.lbl('mlt') +
             notex(' MLT') + ' \\qquad ' + self.lbl('mlat') + 
             notex(' MLAT') + ' \\qquad L_{PP} \\! = \\! ' +
             format(self.lpp, '.1f') )

  # Short, multi-line label for a column label. 
  def lab(self):
    return ( notex( self.probe.upper() ) + '$\n$' + 'L \\! = \\! ' +
             self.lbl('lshell') + '$\n$' + self.lbl('mlt') + '$\n$' +
             self.lbl('mlat') )

  # ---------------------------------------------------------------------------
  # -------------------------------------------- Axis Limits, Ticks, and Labels
  # ---------------------------------------------------------------------------

  def coords(self, x='t', y='b', cramped=False):
    # Assemble a keyword dictionary to be plugged right into the Plot Window. 
    kargs = {}
    # Horizontal axis, time coordinate. 
    if x.lower()=='t':
      kargs['x'] = self.get('t')
      kargs['xlims'] = (self.t0, self.t1)
      kargs['xlabel'] = notex('Time (hh:mm)')
      # If the axis doesn't get the full width of the window, use fewer ticks.
      nxticks = 5 if cramped else 11
      kargs['xticks'] = np.linspace(self.t0, self.t1, nxticks)
      kargs['xticklabels'] = [ '$' + notex( timestr(t)[1] ) + '$' for t in kargs['xticks'] ]
      kargs['xticklabels'][::2] = ['']*len( kargs['xticklabels'][::2] )
    # Horizontal axis, frequency coordinate. 
    elif x.lower()=='f':
      kargs['x'] = self.frq()
      kargs['xlims'] = (0, 40)
      kargs['xlabel'] = notex('Frequency (mHz)')
      # If the axis doesn't get the full width of the window, use fewer ticks.
      nxticks = 5 if cramped else 9
      kargs['xticks'] = np.linspace(0, 40, nxticks)
      kargs['xticklabels'] = [ '$' + znt(t) + '$' for t in kargs['xticks'] ]
      kargs['xticklabels'][1::2] = ['']*len( kargs['xticklabels'][1::2] )
    else:
      print 'UNKNOWN X COORD ', x
    # Vertical axis, plotting electric and/or magnetic fields. 
    if y.lower() in ('e', 'b', 'fields', 'waveform'):
#      kargs['ylabel'] = 'B' + notex(' (nT)  ') + 'E' + notex('(\\frac{mV}{m})')
      kargs['ylabel'] = notex('\\cdots (nT ; \\frac{mV}{m})')
      kargs['ylabelpad'] = -2
      kargs['ylims'] = (-5, 5)
      kargs['yticks'] = range(-5, 6)
      kargs['yticklabels'] = ('', '$-4$', '', '$-2$', '', '$0$', '', '$+2$',
                              '', '$+4$', '')
    # Vertical axis, plotting Fourier magnitude and phase. 
    elif y.lower() in ('p', 's', 'phase', 'spectra'):
      kargs['ylabel'] = '\\cdots' + notex(' (^\\circ)')
      kargs['ylabelpad'] = -2
      kargs['ylims'] = (0, 1)
      kargs['yticks'] = (0, 0.25, 0.5, 0.75, 1)
      kargs['yticklabels'] = ('$-180$', '', '$0$', '', '$+180$')
    elif y.lower() in ('c', 'u', 'coherence', 'unit'):
      kargs['ylabel'] = notex('Coherence')
      kargs['ylims'] = (0, 1)
      kargs['yticks'] = (0, 0.25, 0.5, 0.75, 1)
      kargs['yticklabels'] = ('$0$', '', '', '', '$1$')
    else:
      print 'UNKNOWN Y COORD ', y
    # Return the keyword dictionary, ready to be plugged right into setParams. 
    return kargs

# #############################################################################
# ################################################# High-Level Helper Functions
# #############################################################################

# Estimate the location of the plasmapause at a given timestamp. This is done
# by looking at RBSP inbound and outbound crossing times. 
def lpp(date, time):
  # Read in Scott's list of plasmapause crossings. Make sure they're sorted by
  # time. Each entry is of the form (date, time, L, probe, in/out)
  crossings = sorted( x.split() for x in read('lpp.txt') )
  # Compare the event date and time to the crossing date and times to find the
  # crossings just before and just after this event. 
  evtime = timeint(date, time)
  crosstimes = np.array( [ timeint( *x[0:2] ) for x in crossings ] )
  iafter = np.argmax(crosstimes > evtime)
  ibefore = iafter - 1
  # Get the plasmapause sizes just before and just after the event. 
  lbefore = float( crossings[ibefore][2] )
  lafter = float( crossings[iafter][2] )
  # Return a weighted average. 
  wafter = evtime - timeint( *crossings[ibefore][0:2] )
  wbefore = timeint( *crossings[iafter][0:2] ) - evtime
  return (wbefore*lbefore + wafter*lafter)/(wbefore + wafter)

# Come up with a nice listing of plasmapause crossings. Tab-delimit everything,
# and put the inbound and outbound crossings in the same file. 
def lppcleanup():
  for probe in ('a',):
    for io in ('in', 'out'):
      readpath = ( '/home/user1/mceachern/Desktop/rbsp/lpp/' + probe + '_' +
                   io + 'bound.txt' )
      if not os.path.exists(readpath):
        print 'WARNING: No file found for ' + probe + ' ' + io
        continue
      for line in read(readpath):
        datetime, loc = line.split()
        date, time = datetime.split('/')
        append('\t'.join( (date, time, loc[:5], probe, io) ), 'lpp.txt')
  return

# #############################################################################
# ################################################## Low-Level Helper Functions
# #############################################################################

# Append text to a file. 
def append(text, filename=None):
  if filename is not None:
    with open(filename, 'a') as fileobj:
      fileobj.write(text + '\n')
  return text

# Dot product of a pair of vectors or a pair of arrays of vectors. 
def dot(v0, v1, axis=0):
  return np.sum(v0*v1, axis=axis)

# Load all of the pickles in a directory into a dictionary. 
def load(datadir):
  data = {}
  for pklname in os.listdir(datadir):
    with open(datadir + pklname, 'rb') as handle:
      data[ pklname[:-4] ] = pickle.load(handle)
  return data

# Load a pickle file. 
def loadpickle(pklpath):
  if not os.path.exists(pklpath):
    return None
  with open(pklpath, 'rb') as handle:
    return pickle.load(handle)

# Read in a file as a list of lines. 
def read(filename):
  with open(filename, 'r') as fileobj:
    return [ x.strip() for x in fileobj.readlines() ]

# Returns the time, in seconds, from 1970-01-01. 
def timeint(date=None, time=None):
  # Parse a string of the form hh:mm:ss. 
  if time is None:
    hh, mm, ss = 0, 0, 0
  else:
    # If seconds aren't included, add them. 
    hh, mm, ss = [ int(x) for x in (time + ':00')[:8].split(':') ]
  # Parse a string of the form yyyy-mm-dd. If no date is given, use 1970-01-01
  # so that the returned value is just in seconds from midnight. 
  if date is None:
    year, mon, day = 1970, 1, 1
  else:
    year, mon, day = [ int(x) for x in date.split('-') ]
  return timegm( (year, mon, day, hh, mm, ss) )

# Returns strings indicating the date and time. 
def timestr(ti):
  year, mon, day = gmtime(ti).tm_year, gmtime(ti).tm_mon, gmtime(ti).tm_mday
  hh, mm, ss = gmtime(ti).tm_hour, gmtime(ti).tm_min, gmtime(ti).tm_sec
  date = znt(year, 4) + '-' + znt(mon, 2) + '-' + znt(day, 2)
  time = znt(hh, 2) + ':' + znt(mm, 2) + ':' + znt(ss, 2)
  return date, time

# Scale a vector, or an array of vectors, to unit vectors. 
def unit(v):
  return v/np.sqrt( dot(v, v) )


