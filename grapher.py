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

from calendar import timegm
from time import gmtime

# #############################################################################
# ######################################################################## Main
# #############################################################################

# Global variable indicating where images should be saved. 
savedir = '/home/user1/mceachern/Desktop/plots/' + now() + '/'

def main():

#  for evline in ('b\t2012-10-23\t19:15:00', 'a\t2012-10-23\t21:40:00'):
#    plot0(evline)
#  return

  # Flip through the events in random order. This makes it easier to debug -- a
  # different event pops up first every time. 
  for evline in np.random.permutation( read('events.txt') ):

    if plot1(evline):
      break

#    if plot2(evline):
#      break


  return

# #############################################################################
# ########################################################## Plotting Functions
# #############################################################################

# =============================================================================
# ============================================================== Plot One Event
# =============================================================================

# Plot the waveform. Estimate frequency and phase offset with a Fourier series.
def plot1(evline):
  global savedir

#  # Force it to use a nice event. 
#  evline = 'a\t2014-05-02\t10:10:00'

  # ---------------------------------------------------------------------------
  # ------------------------------------------------- Grab and Check Event Data
  # ---------------------------------------------------------------------------

  # Based on the event record, create an event object. Make sure it's ok. 
  probe, date, time = evline.split()
  ev = event(probe=probe, date=date, time=time)
  if not ev.isok():
    print 'SKIPPING ' + ev.name + ' DUE TO BAD DATA. '
    return False

  # ---------------------------------------------------------------------------
  # -------------------------------------------------------- Set Up Plot Window
  # ---------------------------------------------------------------------------

  # Create a plot window to look at both probes at the same time -- both the
  # waveforms and the spectra. 
  PW = plotWindow(nrows=3, ncols=2)

  # Set the title and labels. 
  title = notex( 'RBSP-' + probe.upper() + '  on  ' + date + '  from  ' +
                 notex( timestr(ev.t0)[1] )  + '  to  ' +
                 notex( timestr(ev.t0)[1] ) )

  title = ev.label()

  rowlabels = ( notex('Poloidal'), notex('Toroidal'), notex('Parallel') )
  collabels = ( notex('B (Red) ; E (Blue)'), notex('Coh. (Black) ; CSD (Green) ; Par. (Violet)') )
  PW.setParams(title=title, collabels=collabels, rowlabels=rowlabels)

  # ---------------------------------------------------------------------------
  # ------------------------------------------------- Add Waveforms and Spectra
  # ---------------------------------------------------------------------------

  # Index the poloidal, toroidal, and field-aligned modes. 
  modes = ('p', 't', 'z')

  # Plot waveforms. 
  PW[:, 0].setParams( **ev.coords('t', 'b', cramped=True) )
  [ PW[i, 0].setLine(ev.get('B' + m), 'r') for i, m in enumerate(modes) ]
  [ PW[i, 0].setLine(ev.get('E' + m), 'b') for i, m in enumerate(modes) ]

  # Plot coherence. 
  PW[:, 1].setParams( **ev.coords('f', 'c', cramped=True) )
  [ PW[i, 1].setLine( ev.csd(m), 'g' ) for i, m in enumerate(modes) ]
  [ PW[i, 1].setLine( ev.coh(m), 'k' ) for i, m in enumerate(modes) ]
  [ PW[i, 1].setLine( 0.05 + 0.9*ev.isodd(m), 'm' ) for i, m in enumerate(modes) ]

  # ---------------------------------------------------------------------------
  # -------------------------------------------------------------- Screen Event
  # ---------------------------------------------------------------------------

  # Let's keep track of all the strong waves we find. 
  waves = []

  # Keep track of MLT. 
  mlt = ev.avg('mlt')

  # Only look at the poloidal mode. Ignore any components of the spectra that
  # are not coherent or not spectrally dense. 
  for m in modes[:1]:
    mname = {'p':'Poloidal', 't':'Toroidal', 'z':'Parallel'}[m]

    # Find any harmonics worth talking about. 
    harm = ev.harm(m)*ev.iscoh(m)*ev.iscsd(m)*ev.ispc4()
    for i in np.nonzero(harm)[0]:
      waves.append( {'strength':ev.coh(m)[i]*ev.csd(m)[i], 
                     'frequency':ev.frq()[i],
                     'harmonic': harm[i], 
                     'mode':mname, 
                     'mlt':mlt} )

  if len(waves) == 0:
    print 'SKIPPING ' + ev.name + ' DUE TO NO SUITABLE MODES. '
    return False

  for w in waves:
    print w

  print 'Estimate of plasmapause location... ', lpp(evline)



  # ---------------------------------------------------------------------------
  # ----------------------------------------------------------- Create the Plot
  # ---------------------------------------------------------------------------

  if '-i' not in argv:
    print 'Plotting ' + ev.name
    return PW.render()
  else:
    if not os.path.exists(savedir):
      os.mkdir(savedir)
    return not PW.render(savedir + ev.name + '.png')

# =============================================================================
# =========================================================== Plot a Long Event
# =============================================================================

# Plot the waveform. Estimate frequency and phase offset with a Fourier series.
def plot0(evline):
  global savedir

  # ---------------------------------------------------------------------------
  # ------------------------------------------------- Grab and Check Event Data
  # ---------------------------------------------------------------------------

  # Based on the event record, create an event object. Make sure it's ok. 
  probe, date, time = evline.split()
  ev = event(probe=probe, date=date, time=time, mins=60)
  if not ev.isok():
    print 'SKIPPING ' + ev.name + ' DUE TO BAD DATA. '
    return False

  # ---------------------------------------------------------------------------
  # -------------------------------------------------------- Set Up Plot Window
  # ---------------------------------------------------------------------------

  PW = plotWindow(nrows=3, ncols=-2)
  title = ev.label()
  rowlabels = ( notex('Poloidal'), notex('Toroidal'), notex('Parallel') )
  PW.setParams(title=title, rowlabels=rowlabels)

  # ---------------------------------------------------------------------------
  # ------------------------------------------------- Add Waveforms and Spectra
  # ---------------------------------------------------------------------------

  # Index the poloidal, toroidal, and field-aligned modes. 
  modes = ('p', 't', 'z')

  # Plot waveforms. 
  PW.setParams( **ev.coords('t', 'b') )
  [ PW[i].setLine(ev.get('B' + m), 'r') for i, m in enumerate(modes) ]
  [ PW[i].setLine(ev.get('E' + m), 'b') for i, m in enumerate(modes) ]

  # ---------------------------------------------------------------------------
  # ----------------------------------------------------------- Create the Plot
  # ---------------------------------------------------------------------------

  if '-i' not in argv:
    print 'Plotting ' + ev.name
    return PW.render()
  else:
    if not os.path.exists(savedir):
      os.mkdir(savedir)
    return not PW.render(savedir + ev.name + '.png')

# =============================================================================
# =========================================== Plot Both Probes During One Event
# =============================================================================

def plot2(evline):
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
  frq = [ ev[p].coh(m)[0] for p, m in pm ]
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

    # Fundamental FLRs should be around 10mHz. 
    if not 6 < f < 20:
      continue

    # Find places where the coherence is significantly above average. 
    iscoherent = np.array( [ c[i] > np.mean(c) + 2*np.std(c) for c in coh ] )

    # In fundamental modes, mlat and phase lag have opposite signs. 
    mlatsign = np.array( [ np.sign( ev[p].avg('mlat') ) for p, m in pm ] )
    lagsign = np.array( [ np.sign( l[i] ) for l in lag ] )

    # Skip anything within one degree of the magnetic equator. 
    isoffeq = np.array( [ np.abs( ev[p].avg('mlat') ) > 1 for p, m in pm ] )

    # Find the frequency/mode combinations -- if any -- which are coherent and
    # fundamental. 
    iscoherentfund = (mlatsign!=lagsign)*iscoherent*isoffeq

#    # Only keep events which trigger in the poloidal mode. 
#    if not any( iscoherentfund[::3] ):
#      continue

    # Can we find any simultaneous observations by both probes? 
    if not any( iscoherentfund[:3] ) or not any( iscoherentfund[3:] ):
      continue

    # Note the nontrivial fundamental modes, if any. 
    for x in np.nonzero(iscoherentfund)[0]:
      modenames = {'p':'Poloidal', 't':'Toroidal', 'z':'Parallel'}
      fundlabels.append( pm[x][0].upper() + '  ' + format(f, '.0f') + 'mHz  ' +
                         modenames[ pm[x][1] ] )

  # If no fundamental modes were found, bail. 
  if not fundlabels:
    print 'SKIPPING ' + ev[probe].name + ' DUE TO NO SUITABLE MODES. '
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

# #############################################################################
# ################################################################ Event Object
# #############################################################################

# Given the name of a directory holding an event, this event grabs the data out
# of that directory and assembles it into a useful form. 
class event:

  # ---------------------------------------------------------------------------
  # --------------------------------------------------- Initialize Event Object
  # ---------------------------------------------------------------------------

  def __init__(self, probe=None, date=None, time=None, t0=None, mins=10):

    # Store event identifiers, and index this event with a unique name. 
    self.probe, self.date = probe, date.replace('-', '')

    if t0 is not None:
      self.t0, self.t1 = t0, t0 + 60*mins
    else:
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
  # ---------------------------------------------------- Create Adjacent Events
  # ---------------------------------------------------------------------------

  # Note that this is not smart enough to handle day boundaries. 
  def prev(self):
    return event(probe=self.probe, date=self.date, t0=self.t0 - 600)

  # Note that this is not smart enough to handle day boundaries. 
  def next(self):
    return event(probe=self.probe, date=self.date, t0=self.t1)

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
      return notex( timestr( 3600*self.avg('mlt') )[1] )
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
  def csd(self, mode, deg=True):
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

  # ---------------------------------------------------------------------------
  # -------------------------------------------------- Spectral Density Filters
  # ---------------------------------------------------------------------------

  # Array indicating which frequencies are coherent. By default, results are
  # significant to 1.5 standard deviations. 
  def iscoh(self, mode, n=1):
    coh = self.coh(mode)
    return np.array( [ c > np.mean(coh) + n*np.std(coh) for c in coh ] )

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
    odds = np.where( self.frq() < 18, 1, 3 )
    return 2*self.iseven(mode) + odds*self.isodd(mode)

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
    return np.array( [ 6 < f < 22 for f in self.frq() ] )

  # ---------------------------------------------------------------------------
  # ---------------------------------------- Labels Indicating Event Properties
  # ---------------------------------------------------------------------------

  # Long label, for the side margin. 
  def label(self):
    return ( notex( 'RBSP--' + self.probe.upper() ) + ' \\qquad L \\!=\\! ' +
             self.lbl('lshell') + ' \\qquad ' + self.lbl('mlt') +
             notex(' MLT') + ' \\qquad ' + self.lbl('mlat') + 
             notex(' MLAT') )

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
# ################################################# High-Level Helper Functions
# #############################################################################

# Estimate the location of the plasmapause at a given timestamp. This is done
# by looking at RBSP inbound and outbound crossing times. 
def lpp(evline):
  # Read in Scott's list of plasmapause crossings. Make sure they're sorted by
  # time. Each entry is of the form (date, time, L, probe, in/out)
  crossings = sorted( x.split() for x in read('lpp.txt') )
  # Compare the event date and time to the crossing date and times to find the
  # crossings just before and just after this event. 
  evtime = timeint( *evline.split()[1:3] )
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

# #############################################################################
# ########################################################### For Importability
# #############################################################################

if __name__=='__main__':
  main()


