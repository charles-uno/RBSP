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
from numpy import pi, where, conj, angle
import os
from plotmod import *

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

'''
  # Location of the output directories. 
  outdir = '/media/My Passport/RBSP/pickles/'
  # Each event is a directory full of pickles. Go through them in random order.
  # This makes us more likely to see bugs when just looking at the first plot.
  for pkldir in np.random.permutation( os.listdir(outdir) ):
#    # If we want to force it to look at a nice one... 
#    pkldir = 'a_20121023_220000'
    # If we're debugging, we just want to stop after a single plot. Some events
    # are bad, so we loop until the plotter tells us it's succeeded. 
    if plotwaveforms(outdir + pkldir + '/'):
      break
  return
'''

# #############################################################################
# ########################################################## Plotting Functions
# #############################################################################

# =============================================================================
# ====================================================== One Event, Both Probes
# =============================================================================

def plotboth(evline):
  global savedir

  # Based on the event record, create a pair of event objects. 
  probe, date, time = evline.split()
  ev = { 'a':event(probe='a', date=date, time=time),
         'b':event(probe='b', date=date, time=time) }

  # If the events have flawed data, start over. 
  if not ev['a'].isok() or not ev['b'].isok():
    print 'SKIPPING ' + ev[probe].name + ' DUE TO BAD DATA. '
    return False

  # Create a plot window to look at both probes at the same time -- both the
  # waveforms and the spectra. 
  PW = plotWindow(nrows=4, ncols=3)

  # Set the title and labels. Label the rows with probe position data. 
  title = notex(probe.upper() + '  ' + date + '  ' + time)
  collabels = ( notex('Poloidal'), notex('Toroidal'), notex('Parallel') )
  rowlabels = ( ev['a'].lab(), ev['b'].lab(), ev['a'].lab(), ev['b'].lab() )
  PW.setParams(title=title, collabels=collabels, rowlabels=rowlabels)

  # Index all probe/mode combinations. 
  pm = [ (p, m) for p in ('a', 'b') for m in ('p', 't', 'z') ]

  # Grab the cross-spectral weights and compute a shared normalization. 
  cw = [ ev[p].crossweights(m) for p, m in pm ]
  norm = max( np.max( np.abs(w) ) for w in cw )

  # Get the spectral magnitudes and phases. 
  mag = [ 2*pi*np.abs(w)/norm for w in cw ]
  ang = [ where( angle(w)<0, angle(w) + 2*pi, angle(w) ) for w in cw ]

  # Iterate over the probe/mode combinations. 
  for i, (p, m) in enumerate(pm):

    # Plot the fields. 
    PW[i/3, i%3].setParams( **ev[p].coords('t', 'b', cramped=True) )
    PW[i/3, i%3].setLine(ev[p].get('E' + m), 'b')
    PW[i/3, i%3].setLine(ev[p].get('B' + m), 'r')

    # Plot the spectra. 
    PW[i/3 + 2, i%3].setParams( **ev[p].coords('f', 's', cramped=True) )
    PW[i/3 + 2, i%3].setLine(mag[i], 'k')
    PW[i/3 + 2, i%3].setLine(ang[i], 'g')

    # Summarize the spectra. 
    frq, phs, pct = ev[p].fpp(m)
    PW[i/3, i%3].setParams(text=frq + ' \\quad ' + phs + ' \\quad ' + pct)

  print 'Plotting ' + ev[probe].name
  return PW.render()

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

# =============================================================================
# ============================================================== Coherence, etc
# =============================================================================

# Plot the waveform. Estimate frequency and phase offset with a Fourier series.
def plotcoherence(datadir):
  global savedir

  # Create an event object from the given path, and make sure the data is ok. 
  ev = event(datadir)
  if not ev.isok():
    print 'SKIPPING ' + ev.name + ' DUE TO BAD DATA. '
    return False

  PW = plotWindow(nrows=2)

  B, E = ev.get('Bx'), ev.get('Ey')
  N = B.size

  t = ev.get('t')
  dt = t[1] - t[0]

#  f = ev.get('f')/1e3

  def cc(x, y, n):
    if n==0:
      return np.sum(x*y)
    else:
      return np.sum( x[n:]*y[:-n] ) if n>0 else np.sum( x[:n]*y[-n:] )

  def csd(x, y, f):
    N = x.size
    return sum( cc(x, y, n)*np.exp(-1j*f*n*dt) for n in range(-N, N) )/(2*pi)

  GBB = np.zeros(ev.get('f').size, dtype=np.complex)
  GEE = np.zeros(ev.get('f').size, dtype=np.complex)
  GEB = np.zeros(ev.get('f').size, dtype=np.complex)

  for i, f in enumerate(ev.get('f')/1e3):
    GBB[i] = csd(B, B, f)
    GEE[i] = csd(E, E, f)
    GEB[i] = csd(E, B, f)
    print '\t', i, '\t', format(f*1e3, '3.0f'), '\t', format(np.abs(GBB[i]), '6.2f'), '\t', format(np.abs(GEE[i]), '6.2f'), '\t', format(np.abs(GEB[i]), '6.2f'), '\t', (np.abs(GEB)**2 / np.abs(GEE*GBB))[i]


  coherence = np.abs(GEB)**2 / np.abs(GEE*GBB)
  PW[0].setParams( ylims=(0, 1), **ev.fparams() )
  PW[0].setLine(coherence)




  # FFT validation. 
  '''
  PW.setParams( **ev.fparams() )
  wB = ev.weights('Bx')
  wE = ev.weights('Ey')
  f = ev.get('f')
  nB = np.fft.rfft( ev.get('Bx') )
  nE = np.fft.rfft( ev.get('Ey') )
  PW[0].setLine(f, np.abs(wB) )
  PW[1].setLine(f, np.abs(wE) )
  PW[0].setLine(f, np.abs(nB) )
  PW[1].setLine(f, np.abs(nE) )
  '''

#  collabels = ( 'E_y' + notex('(Blue) and ') + 'B_x' + notex('(Red)'), notex('Spectral Density') )
#  title = notex('title')

#  PW[1].setParams( ylims=(-5, 5), **ev.tparams(cramped=True) )
#  [ PW[1].setLine(ev.get(v), c) for v, c in ( ('Ey', 'b'), ('Bx', 'r') ) ]

#  PW[0].setParams( ylims=(0, 1), **ev.fparams(cramped=True) )
#  PW[0].setLine(coherence)

#  PW[0].setLine( ev.get('f'), BB )
#  PW[1].setLine( ev.get('f'), EE )
#  PW[2].setLine( ev.get('f'), np.real(EB) )
#  PW[3].setLine( ev.get('f'), np.imag(EB) )

  '''
  wB, wE = ev.weights('Bx'), ev.weights('Ey')
  EE = np.abs(wE)**2
  BB = np.abs(wB)**2
  EB = wE*np.conj(wB)
  coherence = np.abs(EB)**2 / (EE*BB)
  PW[0].setLine(ev.get('f'), EE)
  PW[1].setLine(ev.get('f'), BB)
  PW[2].setLine( ev.get('f'), np.real(EB) )
  PW[3].setLine( ev.get('f'), np.imag(EB) )
  PW[4].setLine(ev.get('f'), coherence)
  PW.setParams( ylims=(0, 2) )
  '''

  # If we're debugging, show the plot, even if it's bad. 
  if '-i' not in argv:
    print 'Plotting ' + ev.name
    return PW.render()
  # If the power is spread out evenly across the harmonics, or if the phase is
  # all over the place, this event is probably garbage. 
  elif int(stdev)>75 or int( powerpct.strip('%\\operatorname{}') )<35:
    print 'SKIPPING ' + ev.name + ' DUE TO AMBIGUOUS PHASE. '
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

  def __init__(self, datadir=None, probe=None, date=None, time=None):


    if datadir is not None:

      # Get the event's name. 
      self.name = datadir.rstrip('/').split('/')[-1]
      # From the name of the directory, get the probe name and the timestamp. 
      self.probe, self.date, self.t0, self.t1 = self.parse(datadir)

    else:

      self.date = date.replace('-', '')
      self.probe = probe

      hh, mm, ss = [ int(x) for x in time.split(':') ]
      self.t0 = 3600*hh + 60*mm + ss
      self.t1 = self.t0 + 600

      self.name = self.date + '_' + time.replace(':', '') + '_' + probe

      datadir = '/media/My Passport/rbsp/pkls/' + self.date + '/' + probe + '/'

    # Load the pickles for time, fields, and position. 
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
  # ----------------------------------------- Label Indicating Event Properties
  # ---------------------------------------------------------------------------

  # Nicely-formatted average probe position. 
  def getpos(self):
    lshell = format(self.avg('lshell'), '.1f')
    mlt = notex( clock( 3600*self.avg('mlt') ) )
    mlat = notex(format(self.avg('mlat'), '+.0f') + '^\\circ')
    return lshell, mlt, mlat

  # Long label, for the side margin. 
  def label(self):
    return ( notex( 'RBSP--' + self.probe.upper() ) + ' \\qquad L \\!=\\! ' +
             format(np.average( self.get('lshell') ), '.1f') + ' \\qquad ' +
             notex( clock( 3600*np.average( self.get('MLT') ) ) +
             ' MLT') + ' \\qquad ' +
             format(np.average( self.get('mlat') ), '+.0f') + '^\\circ' +
             notex(' Magnetic Latitude') )

  # Condensed label, for putting over a cramped frame. 
  def lbl(self):
    return ( notex( self.probe.upper() ) + ' \\quad ' +
             notex(format(self.avg('lshell'), '.1f') + 'R_E') + ' \\quad ' +
             notex(format(self.avg('mlat'), '+.0f') + '^\\circ') + ' \\quad ' +
             notex( clock( 3600*self.avg('mlt') ) + ' MLT' ) )

  # Short, multi-line label for a column label. 
  def lab(self):
    lshell, mlt, mlat = self.getpos()
    return ( notex( self.probe.upper() ) + '$\n$' + 'L \\! = \\! ' + lshell + '$\n$' + mlt + '$\n$' + mlat )

  # Summarize the Fourier transform: the frequency with the maximum
  # cross-spectral weight, the phase of that weight, and the percent of the
  # wave made up by that weight. Format nicely. 
  def fpp(self, mode):
    cw, acw = self.crossweights(mode), np.abs( self.crossweights(mode) )
    frq = notex(format(self.get('f')[ np.argmax(acw) ], '.0f') + 'mHz')
    ang = angle(cw[ np.argmax(acw) ], deg=True)
    phs = notex(format( ang if ang>0 else ang + 360 , '.0f') + '^\\circ')
    pct = notex(format(100*np.max(acw)/np.sum(acw), '.0f') + '\\%')
    return frq, phs, pct

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

  # ---------------------------------------------------------------------------
  # -------------------------------------------- Axis Limits, Ticks, and Labels
  # ---------------------------------------------------------------------------

  def coords(self, x='t', y='b', cramped=False):
    # Assemble a keyword dictionary to be plugged right into the Plot Window. 
    kargs = {}
    # Horizontal axis, time coordinate. 
    if x.lower()=='t':
      kargs['x'] = self.get('tcoord')
      kargs['xlims'] = ( kargs['x'][0], kargs['x'][0] + 10)
      kargs['xlabel'] = notex('Time (hh:mm)')
      # If the axis doesn't get the full width of the window, use fewer ticks.
      if cramped:
        kargs['xticks'] = (0, 2.5, 5, 7.5, 10)
        kargs['xticklabels'] = ['']*len( kargs['xticks'] )
        for i in (0, 2, 4):
          kargs['xticklabels'][i] = '$' + notex( clock(self.t0 + 150*i) ) + '$'
      else:
        kargs['xticks'] = range(11)
        kargs['xticklabels'] = ['']*len( kargs['xticks'] )
        for i in range(1, 11, 2):
          kargs['xticklabels'][i] = '$' + notex( clock(self.t0 + 60*i) ) + '$'
    # Horizontal axis, frequency coordinate. 
    elif x.lower()=='f':
      kargs['x'] = self.get('f')
      kargs['xlims'] = (0, 40)
      kargs['xlabel'] = notex('Frequency (mHz)')
      # If the axis doesn't get the full width of the window, use fewer ticks.
      if cramped:
        kargs['xticks'] = (0, 10, 20, 30, 40)
        kargs['xticklabels'] = ['']*len( kargs['xticks'] )
        kargs['xticklabels'][::2] = ('$0$', '$20$', '$40$')
      else:
        kargs['xticks'] = range(0, 41, 5)
        kargs['xticklabels'] = ['']*len( kargs['xticks'] )
        kargs['xticklabels'][::2] = ('$0$', '$10$', '$20$', '$30$', '$40$')
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
      kargs['ylabel'] = '\\cdots' + notex(' (rad)')
      kargs['ylabelpad'] = -2
      kargs['ylims'] = (0, 2*pi)
      kargs['yticks'] = (0, pi/2, pi, 3*pi/2, 2*pi)
      kargs['yticklabels'] = ('$0$', '', '${\\displaystyle \\pi}$', '',
                              '${\\displaystyle 2 \\pi}$')
    else:
      print 'UNKNOWN Y COORD ', y
    # Return the keyword dictionary, ready to be plugged right into setParams. 
    return kargs

  # ---------------------------------------------------------------------------
  # ------------------------------------------------ Fourier Series of the Data
  # ---------------------------------------------------------------------------

  # Harmonic basis function. Note, awkwardly, that the time steps are not quite
  # uniform. This means that spelling out harmonics explicitly gives slightly
  # better results than using an FFT. 
  def harm(self, m, coord='tcoord'):
    t = self.get(coord)
    return np.exp( m*t*2j*np.pi / (t[-1] - t[0]) )

  # (Complex) Fourier weights for a given field. 
  def weights(self, var):
    wts = np.zeros(self.modes.size, dtype=np.complex)
    # A Fourier transform can be normalized per a handful of different
    # conventions. We choose to match the Numpy FFT convention, where the
    # weights are computed directly then a factor of 1/N is applied when they
    # are reconstituted. 
    for m in self.modes:
      wts[m] = np.sum( self.get(var)*self.harm(-m) )
    return wts

  def crossweights(self, mode):
    return self.weights('E' + mode)*conj( self.weights('B' + mode) )

  # Evaluate the Fourier series to recover a waveform. 
  def series(self, var, real=True):
    wts = self.weights(var)
    srs = np.zeros(self.get('tfine').shape, dtype=np.complex)
    # Sum over all modes. 
    for m in self.modes:
      srs = srs + self.harm(m, coord='tfine')*wts[m]/(self.i1 - self.i0 + 1)
    return np.real(srs) if real is True else srs 

  # Get the power in a given mode based on its Fourier weights. 
  def power(self, *args):
    # Alphabetize the field names to determine magnetic from electric. 
    bvar, evar = sorted(args)
    return np.abs( self.weights(bvar) )**2 + np.abs( self.weights(evar) )**2

  # Compute the complex phase between the electric and magnetic fields. 
  def phase(self, *args):
    # Alphabetize the field names to determine magnetic from electric. 
    bvar, evar = sorted(args)
    bang = np.angle(self.series(bvar, real=False), deg=True)
    eang = np.angle(self.series(evar, real=False), deg=True)
    # Keep it in the domain of -180 to 180. 
    phs = eang - bang
    return np.where(phs<-180, phs + 360, np.where(phs>180, phs - 360, phs) )

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


