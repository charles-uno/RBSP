#!/usr/bin/env python

# Charles McEachern

# Spring 2016

# #############################################################################
# #################################################################### Synopsis
# #############################################################################

# This script makes nice-looking plots of RBSP data. 

# #############################################################################
# ######################################################### Load Python Modules
# #############################################################################

from day import *
from plotmod import *

from numpy.ma import masked_where

# #############################################################################
# ######################################################################## Main
# #############################################################################

# A timestamped directory, in case we want to save any plots.  
plotdir = '/home/user1/mceachern/Desktop/plots/' + now() + '/'

# What should the minimum amplitude be? Anything below there is noise. 
thresh = 0.

def main():

#  return dungey()

#  # Location of the usable data. 
#  return posplot(save='-i' in argv)

#  # Location of all events, by parity and polarization. 
#  return eventplot(save='-i' in argv)

#  # Location of simultaneous poloidal-toroidal events. 
#  return doubleplot(save='-i' in argv)

#  # Location of poloidal events by compressional coupling. 
#  return azmplot(save='-i' in argv)

#  # Location of poloidal events by spectral width. 
#  return [ fwhmplot(mode, split=2, save='-i' in argv) for mode in ('p', 't') ]

  return paramplot(name='fwhm', save='-i' in argv)


  return

# #############################################################################
# #################################### Outer Magnetosphere and the Dungey Cycle
# #############################################################################

'''
def gauss(x, amp=0, avg=np.pi/2, std=1):
  return amp*np.exp( -(x - avg)**2/(2.*std**2) )

def lshell(L, stretch=1., qstretch=0., *args, **kargs):
  q0 = np.arcsin( np.sqrt( 1./np.abs(L) ) )
  q = np.linspace(q0, np.pi - q0, 100)
  r = L*np.sin(q)**2 + gauss(q, *args, **kargs)
  x, z = r*np.sin(q), r*np.cos(q)
  dx = x - x[0]
  xnew = x[0] + stretch*dx + qstretch*dx**2
  return xnew, z

def straighten(x, y):
  n = x.size
  xp, yp = x, y  
  x0, y0 = x[n/2-1], y[n/2-1]
  xp[n/2:] = np.linspace(x0, 30*np.sign(x0), n/2)
  yp[n/2:] = y0
  return xp, yp

def q(L):
  q0 = np.arcsin( np.sqrt( 1./np.abs(L) ) )
  return np.linspace(q0, np.pi - q0, 100)
'''

# L-shell, with a quadratic stretch. 
def lshell(L, qs=0.):
  x = np.linspace(-20, 20, 1000)
  x0 = np.sqrt( np.abs(1./L) ) if L!=0 else 0
  if L>0:
    xp = np.where( x<x0, x, x + qs*(x - x0)**2  )
  else:
    xp = np.where( x>x0, x, x - qs*(x - x0)**2  )
  zz = np.where( np.sign(L)==np.sign(xp), np.abs(L*xp**2)**(2./3) - xp**2, 0)
  z = np.sqrt( np.maximum(zz, 0 ) )
  return x, np.where( x**2 + z**2 > 1, z, 0)

def openline(x0):
  x = np.linspace(-20, 20, 1000)
  c = 0.1*(20 - x0)**2 / np.sqrt(20 - x0)
  if x0 < 0: 
    z = 100*c*np.sqrt( np.maximum(0, x0 - x) )
  else: 
    z = c*np.sqrt( np.maximum(x - x0, 0) )
  return x, z


def dungey():
  PW = plotWindow(ncols=-2, colorbar=None)
  PW.setParams( xlims=(-20, 20), ylims=(-8, 8) )
  PW.setParams(earth='left')
  dstretch = 1.2e-3
  nstretch = 1.2e-3
  [ PW.setLine( *lshell(L, qs=-nstretch*L), color='b' ) for L in range(2, 12, 4) ]
  [ PW.setLine( *lshell(L, qs=-dstretch*L), color='b' ) for L in range(-10, 0, 4) ]
  [ PW.setLine( *openline(x0), color='r' ) for x0 in (-18, -14, 16, 18) ]
  for L, x0 in ( (-12, -10), (-20, -6), (-32, -1.5), (12, 14), (18, 10), (32, 2) ):
    x, zc = lshell(L, qs=-dstretch*L)
    x, zo = openline(x0)
    p = 3.
    PW.setLine(x, (zo**p + zc**p)**(1/p), 'm')
  return PW.render()

  '''
  for i, L in enumerate( (2, 4, 6, 8, 10) ):
    PW.setLine( *lshell(-L), color='b' )
    PW.setLine( *lshell(-L, qstretch=i/200.), color='r' )
    PW.setLine( *lshell(L), color='b' )
    x, z = straighten( *lshell(L, qstretch=i/100., std=0.1, amp=i/10.) )
    PW.setLine( x, z, color='r' )
#  for i, x0 in enumerate( range(0, 20, 2) ):
#    z = (10 - i)*np.sqrt( np.maximum(x - x0, 0)/(20. - x0) )
#    PW.setLine(x, z, 'g')
  z = np.linspace(-8, 8, 100)
  for x0 in range(-20, 20, 2):
    x = x0 + 10*(20 - x0)**-2 * z**2
    PW.setLine(x, z, 'g')
  return PW.render()
#  PW.setLine( *lshell(4, amp=0.5, std=1), color='b' )
#  PW.setLine( *lshell(6, amp=1, std=0.5), color='b' )
#  PW.setLine( *lshell(8,  amp=1, std=0.1), color='b' )
#  PW.setLine( *lshell(10,  amp=2, std=0.1), color='b' )
  x, y0 = straighten( *lshell(12, amp=2, std=0.1) )
  PW.setLine( x, y0, color='m' )
#  y1 = np.sqrt( np.maximum(x - 10, 0) )
#  PW.setLine( x, y1, color='g' )
#  PW.setLine( x, np.sqrt(y0**2 + y1**2), color='orange' )
  y1 = 2*np.sqrt( np.maximum(x - 10, 0)/10. )
  PW.setLine(x, y1, 'r')
  PW.setLine(x, np.sqrt(y0**2 + y1**2), 'g')
  y2 = 1*np.sqrt( np.maximum(x - 14, 0)/6. )
  PW.setLine(x, y2, 'r')
#  for i, x0 in enumerate( (16, 14, 12, 10, 8) ):
#    y = (i+1)*np.sqrt( np.maximum(x - x0, 0)/(20. - x0) )
#    PW.setLine(x, y, 'r')
  PW.render()
  '''
  return

# #############################################################################
# ########################################################## Position Histogram
# #############################################################################

def posplot(save=False, unit='days'):
  global plotdir
  # Grab the position data. 
  unit = 'days'
  pos = getpos(unit=unit, dl=0.5, dm=1)
  x, y, z = [ pos[key] for key in ('x', 'y', 'z') ]
  date0, date1 = pos['dates']
  dt = np.sum(z)
  # Create the plot window using the bullseye params helper function. 
  PW = plotWindow( colorbar='pos', **bep() )
  title = notex('Usable Data ' + date0 + ' to ' + date1 + ' (' + 
                format(dt, '.0f') + ' ' + unit + ' total)')
  PW.setParams( title=title, unitlabel=notex(unit) )
  # Add the data to the plot. 
  PW.setMesh(x, y, z)
  # Show the plot, or save it as an image. 
  if save is True:
    return PW.render(plotdir + 'pos.pdf')
  else:
    return PW.render()

# #############################################################################
# ############################################################ Event Histograms
# #############################################################################

# =============================================================================
# =============================================================== Double Events
# =============================================================================

def doubleplot(save=False):
  global plotdir
  # Set up the grid, and get the probe position for normalization. 
  pos = getpos(dl=3, dm=2, lmin=4)
  x, y, z, hargs = [ pos[key] for key in ('x', 'y', 'z', 'hargs') ]
  d0, d1 = pos['dates']
  # Create a plot window to show different subsets of the events. 
  PW = plotWindow( ncols=2, nrows=2, colorbar='log', ncolors=6, **bep() )
  # Title and labels. 
  title = notex('Double Pc4 Occurrence Rate by Modes, ' + d0 + ' to ' + d1)
  rowlabels = ( notex('Odd\nPoloidal'), notex('Even\nPoloidal') )
  collabels = ( notex('Odd Toroidal'), notex('Even Toroidal') )
  PW.setParams(collabels=collabels, rowlabels=rowlabels, title=title, unitlabel='\\%')
  # Grab double events. We don't actually care which is big and which is small. 
  events = [ line for line in read('events.txt') if 'BIG' in line or 'SMALL' in line ]
  pairs = [ e0 + '\n' + e1 for e0, e1 in zip( events[0::2], events[1::2] ) ]
  # Flip through the possible poloidal-toroidal comparisons. 
  rownames, colnames = ('P1', 'P2'), ('T1', 'T2')
  for row, rn in enumerate(rownames):
    for col, cn in enumerate(colnames):
      # Grab just the pairs that have the poloidal/toroidal combo we're looking for. For each, grab the position. 
      filtered = [ p for p in pairs if rn in p and cn in p ]
      lmm = np.array( [ [ float(g) for g in f.split()[3:6] ] for f in filtered ] )
      # Can't call a histogram with no data. 
      if lmm.size==0:
        h = 0.*z
      else:
        h = 1800*np.histogram2d(lmm[:, 0], lmm[:, 1], **hargs)[0]/86400.
      # Build it into a histogram. Normalize based on sampling. 
      rate = 100*zmask(h)/z
      PW[row, col].setMesh(x, y, rate)
  # Show or save the plot. 
  if save is True:
    return PW.render(plotdir + 'double.pdf')
  else:
    return PW.render()

# =============================================================================
# ========================================================== All Events by Mode
# =============================================================================

def eventplot(save=False):
  global plotdir
  # Set up the grid and 2D histogram based on probe position. 
  pos = getpos(dl=3, dm=2, lmin=4)
  x, y, z, hargs = [ pos[key] for key in ('x', 'y', 'z', 'hargs') ]
  d0, d1 = pos['dates']
  # Create a plot window to show different subsets of the events. 
  PW = plotWindow( ncols=2, nrows=2, colorbar='log', ncolors=6, **bep() )
  # Title and labels. 
  title = notex('Pc4 Occurrence Rate by Mode, ' + d0 + ' to ' + d1)
  collabels = ( notex('Poloidal'), notex('Toroidal') )
  rowlabels = ( notex('Odd\nHarmonic'), notex('Even\nHarmonic') )
  PW.setParams(collabels=collabels, rowlabels=rowlabels, title=title, unitlabel='\\%')
  # Iterate over the filters. 
  mfilt, hfilt = ('P', 'T'), ('1', '2')
  for row, hf in enumerate(hfilt):
    for col, mf in enumerate(mfilt):
      # Grab a histogram of the events, filtered by mode and harmonic. 
      events = getevents(hargs, filt=mf+hf)
      # Normalize by how long each region was sampled. 
      rate = 100*zmask(events)/z # zmask( z, thresh=0.05*np.max(z) )

      # Add the mesh to the plot. 
      PW[row, col].setMesh(x, y, rate)
  # Show or save the plot. 
  if save is True:
    return PW.render(plotdir + 'rate.pdf')
  else:
    return PW.render()

# =============================================================================
# ========================= Poloidal Compressional and Non-Compressional Events
# =============================================================================

# Let's take a look at a plot of how FWHM depends on mode. 
def azmplot(save=False):
  global plotdir
  # Set up the grid and 2D histogram based on probe position. 
  pos = getpos(dl=3, dm=2, lmin=4)
  x, y, z, hargs = [ pos[key] for key in ('x', 'y', 'z', 'hargs') ]
  d0, d1 = pos['dates']
  # Set up the window. 
  PW = plotWindow( ncols=2, nrows=2, colorbar='log', ncolors=6, zmax=10, **bep() )
  title = notex('Poloidal Pc4 by Compressional Coupling, ' + d0 + ' to ' + d1)
  rowlabels = ( notex('Odd\nHarmonic'), notex('Even\nHarmonic') )
  collabels = ( notex('Small ') + 'm', notex('Large ') + 'm' )
  PW.setParams(collabels=collabels, rowlabels=rowlabels, title=title)
  # Iterate over the filters. 
  for row, filt in enumerate( ('P1', 'P2') ):
    # Grab a histogram of the events, filtered by mode and harmonic. 
    events = getevents(hargs, filt=filt, splitcomp=0.2)
    # Normalize by how long each region was sampled. 
    rates = [ 100*zmask(e)/z for e in events ]
    # Add the mesh to the plot. 
    [ PW[row, i].setMesh(x, y, r) for i, r in enumerate(rates) ]
  # Show or save the plot. 
  if save is True:
    return PW.render(plotdir + 'azm_rate.pdf')
  else:
    return PW.render()

# =============================================================================
# ========================================= Poloidal or Toroidal Events by FWHM
# =============================================================================

# Let's take a look at a plot of how FWHM depends on mode. 
def fwhmplot(mode, split=1., save=False):
  global plotdir
  # Set up the grid and 2D histogram based on probe position. 
  pos = getpos(dl=3, dm=2, lmin=4)
  x, y, z, hargs = [ pos[key] for key in ('x', 'y', 'z', 'hargs') ]
  d0, d1 = pos['dates']
  # Set up the window. 
  PW = plotWindow( ncols=2, nrows=2, colorbar='log', ncolors=6, zmax=10, **bep() )
  modename = 'Poloidal' if mode=='p' else 'Toroidal'
  title = notex(modename + ' Pc4 by Spectral Width, ' + d0 + ' to ' + d1)
  rowlabels = ( notex('Odd\nHarmonic'), notex('Even\nHarmonic') )
  collabels = ( notex('FWHM < ' + str(split) + 'mHz'), notex('FWHM > ' + str(split) + 'mHz') )
  PW.setParams(collabels=collabels, rowlabels=rowlabels, title=title)
  # Iterate over the filters. 
  for row, filt in enumerate( ('P1', 'P2') if mode=='p' else ('T1', 'T2') ):
    # Grab a histogram of the events, filtered by mode and harmonic. 
    events = getevents(hargs, filt=filt, splitfwhm=split)
    # Normalize by how long each region was sampled. 
    rates = [ 100*zmask(e)/z for e in events ]
    # Add the mesh to the plot. 
    [ PW[row, i].setMesh(x, y, r) for i, r in enumerate(rates) ]
  # Show or save the plot. 
  if save is True:
    return PW.render(plotdir + 'fwhm_rate.pdf')
  else:
    return PW.render()

# =============================================================================
# ===================================================== Parameter Distributions
# =============================================================================

def peek(arr, name=''):
  print '\t' + name + ' min    = ', np.min(arr)
  print '\t' + name + ' max    = ', np.max(arr)
  print '\t' + name + ' median = ', np.median(arr)
  print '\t' + name + ' mean   = ', np.mean(arr)
  print '\t' + name + ' stdev  = ', np.std(arr)
  return

def getparam(name, filt=''):
  global thresh
  # Grab lines that are the correct mode. 
  evlines = np.array( [ line for line in read('events.txt') if filt in line ] )
  # Filter out anything below an amplitude threshold. 
  amp = np.array( [ float( l.split()[10] ) for l in evlines ] )
  bigenough = np.nonzero(amp - thresh > 0)[0]
  evlines = evlines[bigenough]
  # Match name to column. 
  col = {'f':8, 'fwhm':9, 'amp':10, 'comp':11}[name]
  # Figure out an appropriate range for the histogram. 
  rng = {'f':(7, 25), 'fwhm':(0, 5), 'amp':(0, 1), 'comp':(0, 1)}[name]
  bins = {'f':18, 'fwhm':20, 'amp':10, 'comp':10}[name]
  # Grab the list of values. 
  arr = np.array( [ float( l.split()[col] ) for l in evlines ] )

  peek(arr)

  # Compute the histogram. 
  vals, edges = np.histogram(arr, range=rng, bins=bins)
  # Put points at bin centers. Normalize to the sum for rate. 
  return 0.5*( edges[1:] + edges[:-1] ), vals
#  return 0.5*( edges[1:] + edges[:-1] ), vals*1./np.sum(vals)

# =============================================================================
# ============================================= Look at Parameter Distributions
# =============================================================================

# Let's take a look at a plot of how FWHM depends on mode. 
def paramplot(name, save=False):
  global plotdir
  # Set up the window. 
  PW = plotWindow(ncols=2, nrows=2, colorbar=None)
  ttl = {'f':'Frequency', 'fwhm':'FWHM', 'amp':'Amplitude', 'comp':'Compressional Coupling'}[name]
  title = notex('Pc4 ' + ttl)
  rowlabels = ( notex('Odd\nHarmonic'), notex('Even\nHarmonic') )
  collabels = ( notex('Poloidal'), notex('Toroidal') )
  PW.setParams(collabels=collabels, rowlabels=rowlabels, title=title)
  xlims = { 'f':(7, 25), 'fwhm':(0, 10), 'amp':(0, 1), 'comp':(0, 1) }[name]
  PW.setParams( xlims=xlims, yticklabels=(), ylabel=notex('Rate') )
  # Create a plot window to show different subsets of the events. 
  mfilt, hfilt = ('P', 'T'), ('1', '2')
  # Iterate over the filters. 
  for row, hf in enumerate(hfilt):
    for col, mf in enumerate(mfilt):
      x, y = getparam(name, filt=mf+hf)
      dx = x[1] - x[0]
      ax = PW.cells[row, col].ax
      ax.bar( x - dx/2, y, width=dx )
  # Show or save the plot. 
  if save is True:
    return PW.render(plotdir + name + '.pdf')
  else:
    return PW.render()

# #############################################################################
# ################################################################## Data Input
# #############################################################################

# =============================================================================
# ==================================================================== Position
# =============================================================================

def getpos(dl=0.5, dm=1, lmin=None, lmax=None, unit='days'):
  # The options for units are hours and days. 
  secondsper = 86400. if unit=='days' else 1440.
  # The orbit of both RBSP paths has been broken into five-minute chunks. Grab
  # the position of each chunk which gives good data for E dot B = 0. 
  poslines = [ line for line in read('pos.txt') if line.endswith('ok') ]
  # Get the date range. 
  dates = ( poslines[0].split()[1], poslines[-1].split()[1] )
  # Arrange the positions as an array of floats. 
  pos = np.array( [ [ float(x) for x in p.split()[3:6] ] for p in poslines ] )
  # Figure out the histogram bounds. 
  if lmin is None:
    lmin = np.floor( np.min( pos[:, 0] ) )
  if lmax is None:
    lmax = np.ceil( np.max( pos[:, 0] ) )
  # Center MLT bins on the hour, at least at midnight. 
  mmin, mmax = -dm/2., 24 - dm/2.
  # We want a bin to be centered at zero. That means anything between (24-dm/2)
  # and 24 should be mapped to the range (-dm/2) to 0. 
  posm = np.where( pos[:, 1] > mmax, pos[:, 1] - 24, pos[:, 1] )
  # Number of bins in each direction. 
  lbins, mbins = int( (lmax - lmin)/dl ) + 1, int( (mmax - mmin)/dm ) + 1
  # Keyword arguments for the histogram2d call. 
  hargs = { 'range':( (lmin, lmax), (mmin, mmax) ), 'bins':(lbins-1, mbins-1) }
  # Bin bounds in terms of L and MLT. 
  l, m = np.mgrid[lmin:lmax:lbins*1j, mmin:mmax:mbins*1j]
  # Map to GSE coordinates. Put midnight at the bottom. 
  x, y = l*np.sin(2*pi*m/24.), -l*np.cos(2*pi*m/24.)
  # Bin the position data into a 2D histogram, then scale it to days. 
  h = np.histogram2d(pos[:, 0], posm, **hargs)[0]
  z = 300*h/secondsper
  # Return the position data. Total amount of usable time too. 
  return {'dates':dates, 'l':l, 'm':m, 'x':x, 'y':y, 'z':z, 'hargs':hargs}

# =============================================================================
# ====================================================================== Events
# =============================================================================

# Find bigger-than-average array values. 
def wherebig(arr):
  med = np.median(arr)
  return np.nonzero(arr - med > 0)[0]

# Find smaller-than-average array values. 
def wheresmall(arr):
  med = np.median(arr)
  return np.nonzero(arr - med < 0)[0]

# Get a histogram of event locations, perhaps filtered. 
def getevents(hargs, unit='days', filt='', filt2='', splitcomp=None, splitfwhm=None):
  global thresh
 
  secondsper = 86400. if unit=='days' else 1440.
  # Grab lines that are the correct mode. 
  evlines = np.array( [ line for line in read('events.txt') if filt in line and filt2 in line ] )
  # Filter out anything below an amplitude threshold. 
  amp = np.array( [ float( l.split()[10] ) for l in evlines ] )
  bigenough = np.nonzero(amp - thresh > 0)[0]
  evlines = evlines[bigenough]
  # Grab the position of each event. 
  pos = np.array( [ [ float(x) for x in l.split()[3:6] ] for l in evlines ] )

  # If we're splitting between high and low compresison, do that. 
  if splitcomp is not None:
    comp = np.array( [ float( l.split()[11] ) for l in evlines ] )
    ibig = np.nonzero(comp - splitcomp >= 0)[0]
    ismall = np.nonzero(comp - splitcomp < 0)[0]
    posbig, possmall = pos[ibig], pos[ismall]
    return [ 1800*np.histogram2d(p[:, 0], p[:, 1], **hargs)[0]/secondsper for p in (posbig, possmall) ]

  if splitfwhm is not None:
    fwhm = np.array( [ float( l.split()[9] ) for l in evlines ] )
    ibig = np.nonzero(fwhm - splitfwhm >= 0)[0]
    ismall = np.nonzero(fwhm - splitfwhm < 0)[0]
    posbig, possmall = pos[ibig], pos[ismall]
    return [ 1800*np.histogram2d(p[:, 0], p[:, 1], **hargs)[0]/secondsper for p in (posbig, possmall) ]

  '''
  fwhm = np.array( [ float( l.split()[9] ) for l in evlines ] )
  freq = np.array( [ float( l.split()[8] ) for l in evlines ] )
  amp = np.array( [ float( l.split()[10] ) for l in evlines ] )
  comp = np.array( [ float( l.split()[11] ) for l in evlines ] )
  if only=='bigfwhm':
    pos = pos[ wherebig(fwhm) ]
  elif only=='smallfwhm':
    pos = pos[ wheresmall(fwhm) ]
  elif only=='bigamp':
    pos = pos[ wherebig(amp) ]
  elif only=='smallamp':
    pos = pos[ wheresmall(amp) ]
  '''

  return 1800*np.histogram2d(pos[:, 0], pos[:, 1], **hargs)[0]/secondsper

# #############################################################################
# ############################################################ Helper Functions
# #############################################################################

# Make sure we don't throw anything infinite on the plot. 
def fmask(x):
  return masked_where(np.logical_not( np.isfinite(x) ), x)

# Also leave out any zeros. 
def zmask(x, thr=0):
  return masked_where(np.abs(x)<=thr, x)

# Bullseye params for the plotter. 
def bep():
  tls = ('$-8$', '', '$-4$', '', '$0$', '', '$+4$', '', '$+8$')
  return {'earth':'top', 'flipx':True, 'grid':True, 'square':True, 
          'xlabel': 'Y' + notex(' (R_E)'), 'xlims':(-8, 8),
          'xticks':np.mgrid[-8:8:9j], 'xticklabels':tls, 
          'ylabel': 'X' + notex(' (R_E)'), 'ylims':(-8, 8), 
          'yticks':np.mgrid[-8:8:9j], 'yticklabels':tls, 'ylabelpad':-2 }

# Space out a nice right-justified column. 
def col(x, width=12, unit='', digs=2):
  d = str(digs)
  if isinstance(x, float):
    return format(x, str(width - len(unit) - 1) + '.' + d + 'f') + unit + ' '
  else:
    return str(x).rjust(width - len(unit) - 1) + unit + ' '

# #############################################################################
# ########################################################### For Importability
# #############################################################################

if __name__=='__main__':
  main()


