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
from numpy import pi
import os
from plotmod import *
from scipy import signal

# #############################################################################
# ######################################################################## Main
# #############################################################################

# Global variable indicating where images should be saved. 
savedir = '/home/user1/mceachern/Desktop/plots/' + now() + '/'

def main():

  # Flip through the events in random order. 
  for evline in np.random.permutation( read('events.txt') ):

    if plotboth(evline):
      break

  return

# #############################################################################
# ########################################################## Plotting Functions
# #############################################################################

# =============================================================================
# ====================================================== One Event, Both Probes
# =============================================================================

def plotboth(evline):
  global savedir

#  # Force it to use a nice event. 
#  evline = 'b\t2012-10-23\t19:40:00'

  # ---------------------------------------------------------------------------
  # ------------------------------------------------- Grab and Check Event Data
  # ---------------------------------------------------------------------------

  # Based on the event record, create a pair of event objects. 
  probe, date, time = evline.split()
  ev = { 'a':event(probe='a', date=date, time=time),
         'b':event(probe='b', date=date, time=time) }

  # If the events have flawed data, start over. 
  if not ev['a'].isok() or not ev['b'].isok():
    print 'SKIPPING ' + ev[probe].name + ' DUE TO BAD DATA. '
    return False

  # ---------------------------------------------------------------------------
  # -------------------------------------------------------- Set Up Plot Window
  # ---------------------------------------------------------------------------

  # Create a plot window to look at both probes at the same time -- both the
  # waveforms and the spectra. 
  PW = plotWindow(nrows=4, ncols=3)

  # Set the title and labels. Label the rows with probe position data. 
  title = notex(probe.upper() + '  ' + date + '  ' + time)
  collabels = ( notex('Poloidal'), notex('Toroidal'), notex('Parallel') )
  rowlabels = ( ev['a'].lab(), ev['b'].lab(), ev['a'].lab(), ev['b'].lab() )
  PW.setParams(title=title, collabels=collabels, rowlabels=rowlabels)

  # ---------------------------------------------------------------------------
  # ----------------------------------------- Add Waveforms and Spectra to Plot
  # ---------------------------------------------------------------------------

  # Index all probe/mode combinations. 
  pm = [ (p, m) for p in ('a', 'b') for m in ('p', 't', 'z') ]

  # For each probe and mode, compute coherence and (complex phase) lag. 
  frq = [ 1e3*ev[p].coh(m)[0] for p, m in pm ]
  coh = [ ev[p].coh(m)[1] for p, m in pm ]
  lag = [ ev[p].lag(m)[1] for p, m in pm ]

  # Iterate over the probe/mode combinations. 
  for i, (p, m) in enumerate(pm):

    # Plot the fields. 
    PW[i/3, i%3].setParams( **ev[p].coords('t', 'b', cramped=True) )
    PW[i/3, i%3].setLine(ev[p].get('E' + m), 'b')
    PW[i/3, i%3].setLine(ev[p].get('B' + m), 'r')

    # Plot the spectra. Scale to the unit interval. 
    PW[i/3 + 2, i%3].setParams( **ev[p].coords('f', 'c', cramped=True) )
    PW[i/3 + 2, i%3].setLine( frq[i], coh[i], 'k')
#    PW[i/3 + 2, i%3].setLine( frq[i], lag[i]/360 + 0.5, 'g')

#    # Summarize the spectra. 
#    fr, ph, pc = ev[p].fpp(m)
#    PW[i/3, i%3].setParams(text=fr + ' \\qquad ' + ph)

  # ---------------------------------------------------------------------------
  # -------------------------------------- Look for Fundamental Mode Pc4 Events
  # ---------------------------------------------------------------------------

  fundlabels = []
  for i, f in enumerate( frq[0] ):

    # Ignore frequencies outside the Pc4 band. 
    if not 7 < f < 25:
      continue

    # Find places where the coherence is significantly above average. 
    iscoherent = np.array( [ c[i] > np.mean(c) + np.std(c) for c in coh ] )

    # In fundamental modes, mlat and phase lag have opposite signs. 
    mlatsign = np.array( [ np.sign( ev[p].avg('mlat') ) for p, m in pm ] )
    lagsign = np.array( [ np.sign( l[i] ) for l in lag ] )

    # Find the frequency/mode combinations -- if any -- which are coherent and
    # fundamental. 
    iscoherentfund = (mlatsign!=lagsign)*iscoherent

    # Actually, let's filter just by poloidal components. 
    iscoherentfund = iscoherentfund*np.array( (1, 0, 0, 1, 0, 0) )

    # Note the nontrivial fundamental modes, if any. 
    for x in np.nonzero(iscoherentfund)[0]:
      modenames = {'p':'Poloidal', 't':'Toroidal', 'z':'Parallel'}
      fundlabels.append( pm[x][0].upper() + '  ' + format(f, '.0f') + 'mHz  ' +
                         modenames[ pm[x][1] ] )

  # If no fundamental modes were found, bail. 
  if not fundlabels:
    print 'SKIPPING ' + ev[probe].name + ' DUE TO NO FUNDAMENTAL MODES. '
    return False

  # Otherwise, note the events on the plot. 
  PW.setParams( sidelabel=' $\n$ '.join( notex(x) for x in fundlabels ) )

  # ---------------------------------------------------------------------------
  # ----------------------------------------------------------- Create the Plot
  # ---------------------------------------------------------------------------

  if '-i' not in argv:
    print 'Plotting ' + ev[probe].name
    return PW.render()
  else:
    if not os.path.exists(savedir):
      os.mkdir(savedir)
    return not PW.render(savedir + ev[probe].name + '.png')

# =============================================================================
# ============================================================= Field Waveforms
# =============================================================================

# Plot the waveform. Estimate frequency and phase offset with a Fourier series.
def plotwaveforms(datadir):
  global savedir

  # Create an event object from the given path, and make sure the data is ok. 
  ev = event(datadir)
  if not ev.isok():
    print 'SKIPPING ' + ev.name + ' DUE TO BAD DATA. '
    return False

  # Initialize the Plot Window. Three double-wide rows for poloidal,
  # toroidal, and parallel. 
  ylabels = ( notex('Poloidal (\\frac{mV}{m}) (nT)'), 
              notex('Toroidal (\\frac{mV}{m}) (nT)'), 
              notex('Field-Aligned (\\frac{mV}{m}) (nT)') )
  PW = plotWindow(nrows=len(ylabels), ncols=-2)

  # Manually set the ticks and tick labels on the horizontal axis. 
  PW.setParams( **ev.tparams() )

  # Plot poloidal, toroidal, and field-aligned components. 
  [ PW[0].setLine(ev.get(v), c) for v, c in ( ('Ey', 'b'), ('Bx', 'r') ) ]
  [ PW[1].setLine(ev.get(v), c) for v, c in ( ('Ex', 'b'), ('By', 'r') ) ]
  [ PW[2].setLine(ev.get(v), c) for v, c in ( ('Ez', 'b'), ('Bz', 'r') ) ]

  # Find the peak frequency, in terms of energy density, and pack that
  # information into a label for the plot. 
  pwr = ev.power('Bx', 'Ey')
  pct = notex(format(100*np.max(pwr)/np.sum(pwr), '.0f') + '\\%')
  freq = notex(format( ev.get('f')[ np.argmax(pwr) ], '.0f') + ' mHz')
  flabel = freq + notex(' (') + pct + notex(' Power)')

  # Compute the difference in (complex) phase between the poloidal fields as a
  # function of time. Put the mean and standard deviation in a label. 
  phs = ev.phase('Bx', 'Ey')
  mean = ( '+' + format(np.mean(phs), '.0f') ).replace('+-', '-') + '^\\circ'
  stdev = format(np.std(phs), '.0f') + '^\\circ'
  phslabel = mean + ' \\pm ' + stdev + notex(' Phase Offset')

  # Title the plot and add axis labels. 
  title = notex('In Situ Electric (Blue) and Magnetic (Red) Fields')
  PW.setParams( title=title, collabels=(flabel + ' \\qquad ' + phslabel,), 
                sidelabel=ev.label() )
  [ PW[i].setParams(ylabel=ylbl) for i, ylbl in enumerate(ylabels) ]

#  # For debugging, it's useful to plot the phase (in radians) next to the data.
#  PW[0].setLine(ev.get('tfine'), phs*np.pi/180, 'g')

#  PW[0].setLine(ev.get('tfine'), ev.series('Bx'), 'r:')

  # If we're debugging, show the plot, even if it's bad. 
  if '-i' not in argv:
    print 'Plotting ' + ev.name
    return PW.render()
  # If the power is spread out evenly across the harmonics, or if the phase is
  # all over the place, this event is probably garbage. 
  elif int( stdev.strip('\\^cir') )>99 or int( pct.strip('%\\operatnm{}') )<30:
    print 'SKIPPING ' + ev.name + ' DUE TO AMBIGUOUS PHASE. '
    return False

  # At the moment, we want only the events where the phase lag and the magnetic
  # latitude have opposite signs. 
  elif mean.startswith('+')==( np.mean( ev.get('mlat') )>0 ):
    print 'SKIPPING ' + ev.name + ' DUE TO WRONG SIGN OF MLAT*PHASE. '
    return False

  # Save good plots (after making sure there's a place to put them). 
  else:
    if not os.path.exists(savedir):
      os.mkdir(savedir)
    return not PW.render(savedir + ev.name + '.png')

# #############################################################################
# ################################################################ Event Object
# #############################################################################

# Given the name of a directory holding an event, this event grabs the data out
# of that directory and assembles it into a useful form. 
class event:

  # ---------------------------------------------------------------------------
  # --------------------------------------------------- Initialize Event Object
  # ---------------------------------------------------------------------------

  def __init__(self, probe=None, date=None, time=None, mins=10):

    # Store event identifiers, and index this event with a unique name. 
    self.probe, self.date = probe, date.replace('-', '')
    hh, mm, ss = [ int(x) for x in time.split(':') ]
    self.t0 = 3600*hh + 60*mm + ss
    self.t1 = self.t0 + 60*mins
    self.name = self.date + '_' + time.replace(':', '') + '_' + probe

    # Load the pickles for time, fields, and position. 
    datadir = '/media/My Passport/rbsp/pkls/' + self.date + '/' + probe + '/'
    for var in ('time', 'bgse', 'egse', 'xgse', 'lshell', 'mlt', 'mlat'):
      self.__dict__[var] = loadpickle(datadir + var + '.pkl')

    # Shift the time coordinate to be zero at midnight, rather than in 1970. 
    self.t, self.dt = self.time - self.time[0], self.time[1] - self.time[0]

    # Figure out the indeces that correspond to the event start and end. 
    self.i0, self.i1 = np.argmax(self.t>self.t0), np.argmax(self.t>self.t1)

    # Set the Fourier modes to use. 
    self.modes = np.arange( (self.i1 - self.i0 + 1)/2 )

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
    pt = {'ep':'ey', 'et':'ex', 'bp':'bx', 'bt':'by'}
    if var.lower() in pt:
      return self.get( pt[ var.lower() ] )
    if var=='tcoord':
      return ( self.get('t') - self.get('t')[0] )/60
    elif var=='tfine':
      return np.linspace(self.get('tcoord')[0], self.get('tcoord')[-1], 1000)
    elif var=='f':
      return 1e3*self.modes/600
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
      return clock( 3600*self.avg('mlt') )
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

  # Mode coherence -- at which frequencies do the electric and magnetic field
  # line up for a given mode? Essentially, this is the normalized magnitude of
  # the cross-spectral density. 
  def coh(self, mode):
    t, B, E = self.get('t'), self.get('B' + mode), self.get('E' + mode)
    return signal.coherence(B, E, fs=1/( t[1] - t[0] ), nperseg=t.size/2, 
                            noverlap=t.size/2 - 1)

  # Phase offset -- at each frequency, for each mode, what's the phase lag
  # between of the electric field relative to the magnetic field? This is
  # determined from the complex phase of the cross-spectral density. 
  def lag(self, mode, deg=True):
    t, B, E = self.get('t'), self.get('B' + mode), self.get('E' + mode)
    f, csd = signal.csd(E, B, fs=1/( t[1] - t[0] ), nperseg=t.size/2, 
                            noverlap=t.size/2 - 1)
    return f, np.angle(csd, deg=deg)

  # ---------------------------------------------------------------------------
  # ---------------------------------------- Labels Indicating Event Properties
  # ---------------------------------------------------------------------------

  # Long label, for the side margin. 
  def label(self):
    return ( notex( 'RBSP--' + self.probe.upper() ) + ' \\qquad L \\!=\\! ' +
             self.lbl('lshell') + ' \\qquad ' + self.lbl('mlt') +
             notex(' MLT') + ' \\qquad ' + self.lbl('mlat') + 
             notex(' Magnetic Latitude') )

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
      kargs['xticklabels'] = [ '$' + clock(t) + '$' for t in kargs['xticks'] ]
      kargs['xticklabels'][::2] = ['']*len( kargs['xticklabels'][::2] )
    # Horizontal axis, frequency coordinate. 
    elif x.lower()=='f':
      kargs['x'] = self.get('f')
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
      kargs['ylabel'] = 'B' + notex(' (nT)  ') + 'E' + notex('(\\frac{mV}{m})')
      kargs['ylabelpad'] = -2
      kargs['ylims'] = (-3, 3)
      kargs['yticks'] = range(-3, 4)
      kargs['yticklabels'] = ('$-3$', '', '', '$0$', '', '', '$+3$')
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
# ############################################################ Helper Functions
# #############################################################################

# Put the slashes back into a date string. 
def calendar(date):
  return notex( date[0:4] + '--' + date[4:6] + '--' + date[6:8] )

# Compute clock time from a count of seconds from midnight. 
def clock(sfm, seconds=False):
  if not np.isfinite(sfm):
    return notex('??:??')
  hh = sfm/3600
  if seconds:
    mm = (sfm%3600)/60
    ss = sfm%60
    return notex( znt(hh, 2) + ':' + znt(mm, 2) + ':' + znt(ss, 2) )
  else:
    mm = int( format( (sfm%3600)/60., '.0f') )
    return notex( znt(hh, 2) + ':' + znt(mm, 2) )

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

# Read in a file as a list of lines. 
def read(filename):
  with open(filename, 'r') as fileobj:
    return [ x.strip() for x in fileobj.readlines() ]

# Scale a vector, or an array of vectors, to unit vectors. 
def unit(v):
  return v/np.sqrt( dot(v, v) )

# #############################################################################
# ########################################################### For Importability
# #############################################################################

if __name__=='__main__':
  main()


