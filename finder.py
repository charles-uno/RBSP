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

# #############################################################################
# ######################################################################## Main
# #############################################################################

def main():
  # A date for which we have data, and neighboring day data. 
  return checkdate('2012-10-02')





# Break a given date into chunks and look for fundamental mode Pc4 events in
# that chunk. Each chunk is mpc minutes long -- minutes per chunk. 
def checkdate(d, mpc=20):

  print d

  for p in ('a', 'b'):

    today = day(probe=p, date=d)

    # Scroll through each 20 minute chunk. 
    for t in range(0, 86400, 60*mpc)[:3]:

      print '\t' + p + '\t' + today.timestr(t)[1], 

      ev = today.getslice(t, duration=60*mpc)

      if ev.isok():
        print '\tOK'
      else:
        print '\tX'



# #############################################################################
# ################################################################## Day Object
# #############################################################################

# This class holds one day of RBSP data. Sometimes we want to look at a bunch
# of events in a row, which makes it silly to read in the whole day each time. 

class day:

  # ---------------------------------------------------------------------------
  # ----------------------------------------------------- Initialize Day Object
  # ---------------------------------------------------------------------------

  def __init__(self, probe, date):

    # Identify this day. 
    self.probe, self.date = probe, date

    # Load the position, electric field, and magnetic field data. 
    pickles = self.loadpickles()

    # Grab the position data out of the pickle dictionary. 
    self.lshell = pickles['lshell']
    self.mlt = pickles['mlt']
    self.mlat = pickles['mlat']

    # Offset the time to start at midnight, rather than in 1970. 
    self.t = pickles['time'] - self.timeint(date=self.date)

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
  def loadpickle(self, pklpath):
    if not os.path.exists(pklpath):
      return None
    with open(pklpath, 'rb') as handle:
      return pickle.load(handle)

  # Load all of the pickles for this day. 
  def loadpickles(self):
    pkls = [ x for x in os.listdir( self.path() ) if x.endswith('.pkl') ]
    return dict( ( p[:-4], self.loadpickle(self.path() + p) ) for p in pkls )

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
  # ----------------------------------- Move Between Epoch Time and String Time
  # ---------------------------------------------------------------------------

  # Returns the time, in seconds, from 1970-01-01. 
  def timeint(self, date=None, time=None):
    # Parse a string of the form hh:mm:ss. 
    if time is None:
      hh, mm, ss = 0, 0, 0
    else:
      # If seconds aren't included, add them. 
      hh, mm, ss = [ int(x) for x in (time + ':00')[:8].split(':') ]
    # Parse a string of the form yyyy-mm-dd. If no date is given, use 
    # 1970-01-01 so that the returned value is just in seconds from midnight. 
    if date is None:
      year, mon, day = 1970, 1, 1
    else:
      year, mon, day = [ int(x) for x in date.split('-') ]
    return timegm( (year, mon, day, hh, mm, ss) )

  # Returns strings indicating the date and time. 
  def timestr(self, ti):
    year, mon, day = gmtime(ti).tm_year, gmtime(ti).tm_mon, gmtime(ti).tm_mday
    hh, mm, ss = gmtime(ti).tm_hour, gmtime(ti).tm_min, gmtime(ti).tm_sec
    date = znt(year, 4) + '-' + znt(mon, 2) + '-' + znt(day, 2)
    time = znt(hh, 2) + ':' + znt(mm, 2) + ':' + znt(ss, 2)
    return date, time

  # ---------------------------------------------------------------------------
  # ----------------------------------------------- Slice a Chunk from This Day
  # ---------------------------------------------------------------------------

  # Given a start time and either an end time or a duration, create and return
  # an event object holding that slice of this day. Start and end times can be
  # strings ('hh:mm' or 'hh:mm:ss') or integers (in seconds from midnight). If
  # neither a finish nor a duration is given, the event is 10 minutes long. 
  def getslice(self, start, finish=None, duration=600):
    t0 = self.timeint(start)[1] if isinstance(start, str) else start
    # If a finish time is given, use it. 
    if finish is not None:
      t1 = self.timeint(finish)[1] if isinstance(finish, str) else finish
    # Otherwise, figure out the finish time from the duration. 
    else:
      t1 = t0 + duration
    # Find the indeces that correspond to the desired times. Note that the
    # timestamps don't quite start and end at midnight! 
    i0, i1 = np.argmax(self.t > t0), np.argmax(self.t > t1)
    # Pack up that slice of all the data arrays into a dictionary and use that
    # to create an event object. 
    return event( evdict={ 'probe':self.probe, 'date':self.date, 
                           'lshell':self.lshell[i0:i1], 'mlt':self.mlt[i0:i1], 
                           'mlat':self.mlat[i0:i1], 'bx':self.bx[i0:i1], 
                           'by':self.by[i0:i1], 'bz':self.bz[i0:i1], 
                           'ex':self.bx[i0:i1], 'ey':self.by[i0:i1], 
                           'ez':self.bz[i0:i1], 't':self.t[i0:i1] } )


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

    self.__dict__ = evdict

    return

  # ---------------------------------------------------------------------------
  # ----------------------------------------------------------- Data Validation
  # ---------------------------------------------------------------------------

  # Typically, bad data means that the spin axis was too close to the magnetic
  # field direction, causing problems in the E dot B = 0 assumption. 
  def isok(self):
    return np.all( np.isfinite(self.ex + self.ey) )

  # ---------------------------------------------------------------------------
  # --------------------------------------------------------------- Data Access
  # ---------------------------------------------------------------------------

  # We allow 'Ep' for the poloidal electric field, etc. 
  def get(self, name):
    d = {'bp':'bx', 'bt':'by', 'ep':'ey', 'et':'ex'}
    key = d[ name.lower() ] if name.lower() in d else name.lower()
    return self.__dict__[key]

  # Get the probe's average position during this event. 
  def avg(self, name):
    return np.mean( self.get(name) )

  # ---------------------------------------------------------------------------
  # ---------------------------------------------- Frequency-Domain Data Access
  # ---------------------------------------------------------------------------

  # Frequencies used by the coherence, etc, in mHz. Does not depend on mode. 
  def frq(self, mode='p'):
    t, b, e = self.get('t'), self.get('b' + mode), self.get('e' + mode)
    return signal.coherence(b, e, fs=1/( t[1] - t[0] ), nperseg=t.size/2, 
                            noverlap=t.size/2 - 1)[0]*1e3

  # Fourier transform coefficients. Note that half of the components get
  # tossed, since we want the frequency domain to match with the coherence and
  # cross-spectral density functions. 
  def fft(self, name):
    t, v = self.get('t'), self.get(name)
    # Recall self.frq() gives frequencies in mHz. 
    temp = [ np.sum( v*np.exp(2j*pi*f*t) ) for f in 1e-3*self.frq() ]
    return np.array(temp)/len(temp)

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
  # --------------------------------------------- Frequency-Domain Data Filters
  # ---------------------------------------------------------------------------

  # Array indicating which frequencies are coherent. 
  def iscoh(self, mode, n=1):
    return self.coh(mode) > 0.5

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
  def ispc4(self, mode='p'):
    return np.array( [ 6 < f < 23 for f in self.frq() ] )




  '''
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


  '''



# #############################################################################
# ########################################################### For Importability
# #############################################################################

if __name__=='__main__':
  main()


