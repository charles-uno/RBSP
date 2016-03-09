#!/usr/bin/env python

# Charles McEachern

# Spring 2016

# #############################################################################
# #################################################################### Synopsis
# #############################################################################

# Lei's got some filters in his data which are designed to screen out
# fundamental mode Pc4 pulsations. Can we come up with our own list of events? 
# Specifically, a list of fundamental mode Pc4 pulsations? 

# #############################################################################
# ######################################################### Load Python Modules
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
from plotmod import notex

# #############################################################################
# ################################################################## Day Object
# #############################################################################

# This class holds one day of RBSP data. Sometimes we want to look at a bunch
# of events in a row, which makes it silly to read in the whole day each time. 

class day:

  garbage = False

  # ---------------------------------------------------------------------------
  # ----------------------------------------------------- Initialize Day Object
  # ---------------------------------------------------------------------------

  def __init__(self, probe, date):

    # Identify this day. 
    self.probe = probe
    self.date = timestr( timeint(date=date) )[0]

    # Load the pickles, and make sure they contain actual data. 
    pickles = self.loadpickles()
    if not all( x.size > 1 for x in pickles.values() ):
      self.garbage = True
      print 'WARNING: No data for ' + self.date + ' RBSP-' + self.probe.upper()
      return

    # Grab the position data out of the pickle dictionary. 
    self.lshell = pickles['lshell']
    self.mlt = pickles['mlt']
    self.mlat = pickles['mlat']

    # Offset the time to start at midnight, rather than in 1970. 
    self.t = pickles['time'] - timeint(date=self.date)

    # Compute the background magnetic field using a rolling average. 
    b0gse = self.rollavg( pickles['bgse'] )

    # Compute dipole coordinate basis vectors. 
    self.xhat, self.yhat, self.zhat = self.uvecs( pickles['xgse'], b0gse )

    # Rotate the fields from GSE coordinates to dipole coordinates. 
    self.bx, self.by, self.bz = self.rotate( pickles['bgse'] - b0gse )
    self.ex, self.ey, self.ez = self.rotate( pickles['egse'] )

    return

  # ---------------------------------------------------------------------------
  # ----------------------------------------------------------------- Load Data
  # ---------------------------------------------------------------------------

  # Path to the pickles corresponding to this day. 
  def path(self):
    return ( '/media/My Passport/rbsp/pkls/' + self.date.replace('-', '') +
             '/' + self.probe + '/' )

  # Load a pickle file. 
  def loadpickle(self, name):
    if not os.path.exists(self.path() + name + '.pkl'):
      return None
    with open(self.path() + name + '.pkl', 'rb') as handle:
      return pickle.load(handle)

  # Load all of the pickles for this day. 
  def loadpickles(self):
    pkls = [ x for x in os.listdir( self.path() ) if x.endswith('.pkl') ]
    return dict( ( p[:-4], self.loadpickle( p[:-4] ) ) for p in pkls )

  # ---------------------------------------------------------------------------
  # ----------------------------------------- Compute Background Magnetic Field
  # ---------------------------------------------------------------------------

  # A (by default) ten-minute rolling average is used to estimate the
  # background magnetic field. For now, the range is truncated at midnight. 
  def rollavg(self, b, tlen=600):
    b0 = np.zeros(b.shape)
    di = int( round( tlen*0.5/( self.t[1] - self.t[0] ) ) )
    for i in range( b.shape[1] ):
      # Find indeces corresponding to tlen/2 in the past and in the future. 
      ipast, ifut = max( 0, i - di ), min( i + di, b.shape[1] )
      # Take the average over that range. 
      b0[:, i] = np.average(b[:, ipast:ifut], axis=1)
    return b0

  # ---------------------------------------------------------------------------
  # ----------------------------------------------------- Compute Basis Vectors
  # ---------------------------------------------------------------------------

  # Scale a vector, or an array of vectors, to unit vectors. 
  def unit(self, v):
    return v/np.sqrt( self.dot(v, v) )

  # Dot product of a pair of vectors or a pair of arrays of vectors. 
  def dot(self, v, w, axis=0):
    return np.sum(v*w, axis=axis)

  # Compute the dipole coordinate directions. The zhat unit vector lines up
  # with the background magnetic field, yhat is azimuthally eastward, and xhat
  # completes the orthonormal coordinate system. 
  def uvecs(self, x, b0):
    zhat = self.unit(b0)
    yhat = self.unit( np.cross(zhat, x, axis=0) )
    return np.cross(yhat, zhat, axis=0), yhat, zhat

  # ---------------------------------------------------------------------------
  # ------------------------------------- Rotate from GSE to Dipole Coordinates
  # ---------------------------------------------------------------------------

  # Dot vectors in GSE coordinates with the basis vectors (also in GSE
  # coordinates) to end up with vector components in the dipole basis. 
  def rotate(self, vgse):
    return [ self.dot(vgse, w) for w in (self.xhat, self.yhat, self.zhat) ]

  # ---------------------------------------------------------------------------
  # ----------------------------------------------- Slice a Chunk from This Day
  # ---------------------------------------------------------------------------

  # Given a start time and either an end time or a duration, create and return
  # an event object holding that slice of this day. Start and end times can be
  # strings ('hh:mm' or 'hh:mm:ss') or integers (in seconds from midnight). If
  # neither a finish nor a duration is given, the event is 10 minutes long. 
  def getslice(self, start, finish=None, duration=600):
    t0 = timeint(time=start) if isinstance(start, str) else start
    # If a finish time is given, use it. 
    if finish is not None:
      t1 = timeint(time=finish) if isinstance(finish, str) else finish
    # Otherwise, figure out the finish time from the duration. 
    else:
      t1 = t0 + duration
    # If the data for this day is garbage, return an event that knows just
    # enough to tell you it's broken. 
    if self.garbage:
      return event( evdict={'probe':self.probe, 'date':self.date, 
                            'time':timestr(t0)[1], 't':np.array(None) } )
    # Find the indeces that correspond to the desired times. Note that the
    # timestamps don't quite start and end at midnight! 
    i0, i1 = np.argmax(self.t > t0), np.argmax(self.t > t1)
    # Pack up that slice of all the data arrays into a dictionary and use that
    # to create an event object. 
    return event( evdict={ 'probe':self.probe, 'date':self.date, 
                           'lshell':self.lshell[i0:i1], 'mlt':self.mlt[i0:i1], 
                           'mlat':self.mlat[i0:i1], 'bx':self.bx[i0:i1], 
                           'by':self.by[i0:i1], 'bz':self.bz[i0:i1], 't0':t0, 
                           'ex':self.ex[i0:i1], 'ey':self.ey[i0:i1], 't1':t1,
                           'ez':self.ez[i0:i1], 't':self.t[i0:i1], 
                           'time':timestr(t0)[1] } )


# #############################################################################
# ################################################################ Event Object
# #############################################################################

# Unlike the past incarnation, the event object can no longer read in data to
# create itself. Inctead, it must be created using the day.getevent() method. 
# The old convention had us reading in a pickle and throwing away all but ten
# minutes of it. That's quite inefficient when the intention is to iterate over
# an entire day! 
class event:

  # ---------------------------------------------------------------------------
  # --------------------------------------------------- Initialize Event Object
  # ---------------------------------------------------------------------------

  def __init__(self, evdict):

    # Initialize the event parameters and arrays from the day's data. 
    self.__dict__ = evdict

    # Uniquely identify this event. 
    self.name = ( self.date.replace('-', '') + '_' + 
                  self.time.replace(':', '') + '_' + self.probe )

    '''
    print self.name, '\t', self.t.size
    '''

    # Estimate the location of the plasmapause during this event. 
    self.lpp = lpp(date=self.date, time=self.time)

  # ---------------------------------------------------------------------------
  # ----------------------------------------------------------- Data Validation
  # ---------------------------------------------------------------------------

  # Typically, bad data means that the spin axis was too close to the magnetic
  # field direction, causing problems in the E dot B = 0 assumption. 
  def isok(self):
    return self.t.size > 1 and np.all( np.isfinite(self.ex + self.ey) )

  # ---------------------------------------------------------------------------
  # --------------------------------------------------------------- Data Access
  # ---------------------------------------------------------------------------

  # We allow 'Ep' for the poloidal electric field, etc. Also, for fields,
  # subtract off the DC offset. 
  def get(self, name):
    keys = {'bp':'bx', 'bt':'by', 'bz':'bz', 'ep':'ey', 'et':'ex', 'ez':'ez'}
    if name.lower() in keys:
      arr = self.__dict__[ keys[ name.lower() ] ]
      return arr - np.mean(arr)
    else:
      return self.__dict__[ name.lower() ]

  # Get the probe's average position during this event. 
  def avg(self, name):
    return np.mean( self.get(name) )

  # ---------------------------------------------------------------------------
  # ---------------------------------------------- Frequency-Domain Data Access
  # ---------------------------------------------------------------------------

  # Frequencies used by the coherence, etc, in mHz. Does not depend on mode. 
  def frq(self):
    dt = ( self.t[1] - self.t[0] )*1e-3
    return np.arange(0, self.t.size/4 + 1)/(0.5*self.t.size*dt)

  # Fourier transform coefficients. Note that half of the components get
  # tossed, since we want the frequency domain to match with the coherence and
  # cross-spectral density functions. Return only magnitude. 
  def fft(self, name):
    v = self.get(name)
    # Recall self.frq() gives frequencies in mHz. 
    temp = [ np.sum( v*np.exp(2j*pi*f*self.t) ) for f in 1e-3*self.frq() ]
    return np.abs( np.array(temp)/len(temp) )

  # Mode coherence -- at which frequencies do the electric and magnetic field
  # line up for a given mode? Essentially, this is the normalized magnitude of
  # the cross-spectral density. 
  def coh(self, mode):
    t, b, e = self.get('t'), self.get('b' + mode), self.get('e' + mode)
    return signal.coherence(b, e, fs=1/( t[1] - t[0] ), nperseg=t.size/2, 
                            noverlap=t.size/2 - 1)[1]

  # Absolute value of cross-spectral density, normalized to the unit interval. 
  # The units are sorta arbitrary, so we're not really losing any information. 
  # Note that the angle can be recovered from the phase lag function below. 
  def csd(self, mode):
    t, b, e = self.get('t'), self.get('b' + mode), self.get('e' + mode)
    temp = signal.csd(b, e, fs=1/( t[1] - t[0] ), nperseg=t.size/2, 
                      noverlap=t.size/2 - 1)[1]
    return np.abs(temp)/np.max( np.abs(temp) )

  # Phase offset -- at each frequency, for each mode, what's the phase lag
  # between of the electric field relative to the magnetic field? This is
  # determined from the complex phase of the cross-spectral density. 
  def lag(self, mode, deg=True):
    t, b, e = self.get('t'), self.get('b' + mode), self.get('e' + mode)
    return np.angle(signal.csd(b, e, fs=1/( t[1] - t[0] ), nperseg=t.size/2, 
                      noverlap=t.size/2 - 1)[1], deg=True)

  # ---------------------------------------------------------------------------
  # --------------------------------------------------- Frequency-Domain Ratios
  # ---------------------------------------------------------------------------

  # The ratio of the poloidal and compressional magnetic fields distinguishes
  # between high and low azimuthal modenumber. 
  def comp(self):
    return np.abs( self.fft('bz')/self.fft('bx') )

  # The ratio between the electric and magnetic field is a kludgey estimate of
  # the Alfven speed. 
  def va(self, mode='p'):
    return np.abs( self.fft('e' + mode)/self.fft('b' + mode) )

  # ---------------------------------------------------------------------------
  # --------------------------------------------- Frequency-Domain Data Filters
  # ---------------------------------------------------------------------------

  # Array indicating which frequencies are coherent. 
  def iscoh(self, mode, n=1):
    return self.coh(mode) > 0.75

  # Gives which Fourier components have B (E) above 0.25 nT (mV/m). 
  def has(self, name):
    return np.abs( self.fft(name) ) > 0.25

  # Array indicating which frequencies are spectrally dense. By default, 
  # results are significant to 1.5 standard deviations. 
  def iscsd(self, mode, n=1):
    csd = self.csd(mode)
    return np.array( [ c > np.mean(csd) + n*np.std(csd) for c in csd ] )

  # Array of integers guessing a harmonic for each frequency. The phase lag
  # between electric and magnetic fields gives the parity. If mlat is within a
  # few degrees of the equator, this function returns all zeros. 
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
  def ispc4(self):
    # Lei says 25mHz, but 22mHz might be more typical. 
    return np.array( [ 6 < f < 26 for f in self.frq() ] )

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------- Label Describing Event
  # ---------------------------------------------------------------------------

  # Nicely format a property of this event. 
  def lbl(self, name):
    if name.lower()=='lshell':
      return format(self.avg('lshell'), '.1f')
    elif name.lower()=='mlt':
      return notex( timestr( 3600*self.avg('mlt') )[1][:5] )
    elif name.lower()=='mlat':
      return notex(format(self.avg('mlat'), '+.0f') + '^\\circ')
    else:
      return '???'

  # Assemble event properties into a label. 
  def label(self):
    return notex( self.date + '\\quad' + self.time + '\\quad ' + 'RBSP-' +
                  self.probe.upper() + '\\quad{}L\\!=\\!' +
                  self.lbl('lshell') + '\\quad' + self.lbl('mlt') + 'MLT' +
                  '\\quad' + self.lbl('mlat') + 'MLAT\\quad{}L_{PP}\\!=\\!' +
                  format(self.lpp, '.1f') )

  # Label describing the best wave in the given mode, selected based on electric
  # field FFT magnitude. 
  def descr(self, mode):
    ipeak = np.argmax( self.fft('e' + mode) )
    return notex( format(self.fft('ey')[ipeak], '.1f') + '\\frac{mV}{m}' + 
                  '\\qquad' + format(self.frq()[ipeak], '.0f') + 'mHz\\qquad' +
                  format(self.lag('p')[ipeak], '+.0f') + '^\\circ\\qquad' +
                 '|\\overset{\\sim}{B}_z/\\overset{\\sim}{B}_x|\\!=\\!' +
                 format(100*self.comp()[ipeak], '.0f') + '\\%\\qquad' + 
                 '|\\overset{\\sim}{E}_y/\\overset{\\sim}{B}_x|\\!=\\!' +
                 format(1e3*self.va('p')[ipeak], '.0f')  + '\\frac{km}{s}' )

  # ---------------------------------------------------------------------------
  # -------------------------------------------- Axis Limits, Ticks, and Labels
  # ---------------------------------------------------------------------------

  def coords(self, style, cramped=False):
    # Assemble a keyword dictionary to be plugged right into the Plot Window. 
    kargs = {}
    # Plot electric and magnetic fields as a function of time. 
    if style.lower().startswith('wave'):
      # Horizontal axis. 
      kargs['x'] = self.get('t')
      kargs['xlims'] = (self.t0, self.t1)
      kargs['xlabel'] = notex('Time (hh:mm)')
      nxticks = 5 if cramped else 11
      kargs['xticks'] = np.linspace(self.t0, self.t1, nxticks)
      kargs['xticklabels'] = [ '$' + notex( timestr(t)[1] ) + '$' for t in kargs['xticks'] ]
      kargs['xticklabels'][::2] = ['']*len( kargs['xticklabels'][::2] )
      # Vertical axis. 
      kargs['ylabel'] = notex('\\cdots (nT ; \\frac{mV}{m})')
      kargs['ylabelpad'] = -2
      kargs['ylims'] = (-5, 5)
      kargs['yticks'] = range(-5, 6)
      kargs['yticklabels'] = ('', '$-4$', '', '$-2$', '', '$0$', '', '$+2$',
                              '', '$+4$', '')
    # Plot coherence, cross-spectral density, etc in the frequency domain. 
    else:
      # Horizontal axis. 
      kargs['x'] = self.frq()
      kargs['xlims'] = (0, 40)
      kargs['xlabel'] = notex('Frequency (mHz)')
      nxticks = 5 if cramped else 9
      kargs['xticks'] = np.linspace(0, 40, nxticks)
      kargs['xticklabels'] = [ '$' + znt(t) + '$' for t in kargs['xticks'] ]
      kargs['xticklabels'][1::2] = ['']*len( kargs['xticklabels'][1::2] )
      # Vertical axis. 
      kargs['ylabel'] = notex('Coherence')
      kargs['ylims'] = (0, 1)
      kargs['yticks'] = (0, 0.25, 0.5, 0.75, 1)
      kargs['yticklabels'] = ('$0$', '', '', '', '$1$')
    return kargs

# #############################################################################
# ############################################################ Helper Functions
# #############################################################################

# Append text to a file. 
def append(text, filename=None):
  if filename is not None:
    with open(filename, 'a') as fileobj:
      fileobj.write(text + '\n')
  return text

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
    # Account for missing colons and/or missing seconds. 
    hms = ( time.replace(':', '') + '00' )[:6]
    hh, mm, ss = int( hms[0:2] ), int( hms[2:4] ), int( hms[4:6] )
  # Parse a string of the form yyyy-mm-dd. If no date is given, use 
  # 1970-01-01 so that the returned value is just in seconds from midnight. 
  if date is None:
    year, mon, day = 1970, 1, 1
  else:
    # Allow missing dashes in the date. 
    ymd = date.replace('-', '')
    year, mon, day = int( ymd[0:4] ), int( ymd[4:6] ), int( ymd[6:8] )
  return timegm( (year, mon, day, hh, mm, ss) )

# Returns strings indicating the date and time. 
def timestr(ti):
  year, mon, day = gmtime(ti).tm_year, gmtime(ti).tm_mon, gmtime(ti).tm_mday
  hh, mm, ss = gmtime(ti).tm_hour, gmtime(ti).tm_min, gmtime(ti).tm_sec
  date = znt(year, 4) + '-' + znt(mon, 2) + '-' + znt(day, 2)
  time = znt(hh, 2) + ':' + znt(mm, 2) + ':' + znt(ss, 2)
  return date, time



