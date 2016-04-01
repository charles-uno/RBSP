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
thresh = 0.01 if 'thresh' in argv else 0.001

# Bin all position plots the same way. 
pargs = {'dl':0.5, 'dm':1} if 'sharp' in argv else {'dl':3, 'dm':2, 'lmin':4}

def main():

#  return dungey()

#  # Location of the usable data. 
#  return posplot(save='-i' in argv)

#  # Location of all events, regardless of mode or harmonic. 
#  return allplot(save='-i' in argv)

  # Location of events inside or outside the plasmapause. 
  return [ llppplot(mode, save='-i' in argv) for mode in ('p', 't') ]

#  # Location of all events, by parity and polarization. 
#  return modeplot(save='-i' in argv)

#  # Location of simultaneous poloidal-toroidal events. 
#  return doubleplot(save='-i' in argv)

#  # Location of poloidal events by compressional coupling. 
#  return azmplot(save='-i' in argv)

  # Location of poloidal events by spectral width. 
  return [ fwhmplot(mode, split=1.5, save='-i' in argv) for mode in ('p', 't') ]

  return paramplot(name='fwhm', save='-i' in argv)

  return

# #############################################################################
# ########################################################## Position Histogram
# #############################################################################

def posplot(save=False):
  global plotdir
  # Grab the position data. 
  pos = getpos(dl=0.5, dm=1)
  x, y, z = [ pos[key] for key in ('x', 'y', 'z') ]
  date0, date1 = pos['dates']
  dt = np.sum(z)
  # Create the plot window using the bullseye params helper function. 
  PW = plotWindow( **bep() )
  title = notex('Usable Data ' + date0 + ' to ' + date1 + ' (' + 
                format(dt, '.0f') + ' days total)')
  PW.setParams( title=title, unitlabel=notex('days') )
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
# ========================================================== All Events by Mode
# =============================================================================

def allplot(save=False):
  global pargs, plotdir
  # Set up the grid and 2D histogram based on probe position. 
  pos = getpos(**pargs)
  x, y, z, hargs = [ pos[key] for key in ('x', 'y', 'z', 'hargs') ]
  # Create a plot window to show different subsets of the events. 
  PW = plotWindow( **bep() )
  # Title and labels. 
  title = notex( 'Pc4 Occurrence Rate: ' + tex('ImS') +
                 ' \\geq ' + str(thresh) + tex('mW/m^2') )
  PW.setParams(title=title, unitlabel='\\%')
  # Grab the events histogram. 
  eh = eventhist(hargs)
  # Normalize it by the sampling rate and plot it. 
  PW.setMesh(x, y, 100*zmask(eh)/z)

  print 'total percent: ' + format(100*np.sum(eh)/np.sum(z), '.1f') + '%'

  # Show or save the plot. 
  if save is True:
    return PW.render(plotdir + 'rate_all.pdf')
  else:
    return PW.render()

# =============================================================================
# =============================================================== Double Events
# =============================================================================

def doubleplot(save=False):
  global pargs, plotdir, thresh
  # Set up the grid, and get the probe position for normalization. 
  pos = getpos(**pargs)
  x, y, z, hargs = [ pos[key] for key in ('x', 'y', 'z', 'hargs') ]
  # Create a plot window to show different subsets of the events. 
  PW = plotWindow( ncols=2, nrows=2, **bep() )
  # Title and labels. 
  title = notex( 'Double Pc4 Occurrence Rate by Modes: ' + tex('ImS') +
                 ' \\geq ' + str(thresh) + tex('mW/m^2') )
  rlabs = ( notex('Odd\nPoloidal'), notex('Even\nPoloidal') )
  clabs = ( notex('Odd Toroidal'), notex('Even Toroidal') )
  PW.setParams(collabels=clabs, rowlabels=rlabs, title=title, unitlabel='\\%')
  # Flip through the possible poloidal-toroidal comparisons. 
  rownames, colnames = ('P1', 'P2'), ('T1', 'T2')
  for row, rn in enumerate(rownames):
    for col, cn in enumerate(colnames):
      # Grab a histogram of appropriately-filtered double events. 
      dh = doublehist(hargs, pmode=rn, tmode=cn)
      # Build it into a histogram. Normalize based on sampling. 
      rate = 100*zmask(dh)/z
      PW[row, col].setMesh(x, y, rate)
  # Show or save the plot. 
  if save is True:
    return PW.render(plotdir + 'double_rate.pdf')
  else:
    return PW.render()

# =============================================================================
# ========================================================== All Events by Mode
# =============================================================================

def modeplot(save=False):
  global pargs, plotdir
  # Set up the grid and 2D histogram based on probe position. 
  pos = getpos(**pargs)
  x, y, z, hargs = [ pos[key] for key in ('x', 'y', 'z', 'hargs') ]
  # Create a plot window to show different subsets of the events. 
  PW = plotWindow( ncols=2, nrows=2, **bep() )
  # Title and labels. 
  title = notex( 'Pc4 Occurrence Rate by Mode: ' + tex('ImS') +
                 ' \\geq ' + str(thresh) + tex('mW/m^2') )
  collabels = ( notex('Poloidal'), notex('Toroidal') )
  rowlabels = ( notex('Odd\nHarmonic'), notex('Even\nHarmonic') )
  PW.setParams(collabels=collabels, rowlabels=rowlabels, title=title, unitlabel='\\%')
  # Iterate over the filters. 
  mfilt, hfilt = ('P', 'T'), ('1', '2')
  for row, hf in enumerate(hfilt):
    for col, mf in enumerate(mfilt):
      # Grab a histogram of the events, filtered by mode and harmonic. 
      eh = eventhist(hargs, mode=mf+hf)
      # Normalize by how long each region was sampled. 
      rate = 100*zmask(eh)/z

      print 'overall percent for ' + mf + hf + ' = ' + format(100*np.sum(eh)/np.sum(z), '.1f') + '%'

      # Add the mesh to the plot. 
      PW[row, col].setMesh(x, y, rate)
  # Show or save the plot. 
  if save is True:
    return PW.render(plotdir + 'mode_rate.pdf')
  else:
    return PW.render()

# =============================================================================
# ========================= Poloidal Compressional and Non-Compressional Events
# =============================================================================

# Let's take a look at a plot of how FWHM depends on mode. 
def azmplot(save=False):
  global pargs, plotdir
  # Set up the grid and 2D histogram based on probe position. 
  pos = getpos(**pargs)
  x, y, z, hargs = [ pos[key] for key in ('x', 'y', 'z', 'hargs') ]
  # Set up the window. 
  PW = plotWindow( ncols=2, nrows=2, **bep() )
  title = notex('Poloidal Pc4 by Compressional Coupling: ' + tex('ImS') +
                 ' \\geq ' + str(thresh) + tex('mW/m^2') )
  rowlabels = ( notex('Odd\nHarmonic'), notex('Even\nHarmonic') )
  collabels = ( notex('Small ') + 'm', notex('Large ') + 'm' )
  PW.setParams(collabels=collabels, rowlabels=rowlabels, title=title)
  # Iterate over the filters. 
  for row, mode in enumerate( ('P1', 'P2') ):
    # Grab a histogram of the events, filtered by spectral width. 
    smallm  = eventhist(hargs, mode=mode, comp_ge=0.2)
    bigm = eventhist(hargs, mode=mode, comp_lt=0.2)
    # Normalize by how long each region was sampled. 
    rates  = 100*zmask(smallm)/z, 100*zmask(bigm)/z

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
  global pargs, plotdir
  # Set up the grid and 2D histogram based on probe position. 
  pos = getpos(**pargs)
  x, y, z, hargs = [ pos[key] for key in ('x', 'y', 'z', 'hargs') ]
  # Set up the window. 
  PW = plotWindow( ncols=2, nrows=2, **bep() )
  modename = 'Poloidal' if mode=='p' else 'Toroidal'
  title = notex(modename + ' Pc4 by Spectral Width: ' + tex('ImS') +
                 ' \\geq ' + str(thresh) + tex('mW/m^2') )
  rowlabels = ( notex('Odd\nHarmonic'), notex('Even\nHarmonic') )
  collabels = ( notex('FWHM') + ' < ' + str(split) + notex('mHz'), notex('FWHM') + ' > ' + str(split) + notex('mHz') )
  PW.setParams(collabels=collabels, rowlabels=rowlabels, title=title, unitlabel='\\%')
  # Iterate over the filters. 
  for row, mname in enumerate( ('P1', 'P2') if mode=='p' else ('T1', 'T2') ):
    # Grab a histogram of the events, filtered by spectral width. 
    broad  = eventhist(hargs, mode=mname, fwhm_lt=split)
    narrow = eventhist(hargs, mode=mname, fwhm_ge=split)

    print mname + ' overall broad  rate: ' + format(100*np.sum(broad)/np.sum(z), '.1f') + '%'
    print mname + ' overall narrow rate: ' + format(100*np.sum(narrow)/np.sum(z), '.1f') + '%'

    # Normalize by how long each region was sampled. 
    rates  = 100*zmask(broad)/z, 100*zmask(narrow)/z
    # Add the mesh to the plot. 
    [ PW[row, i].setMesh(x, y, r) for i, r in enumerate(rates) ]
  # Show or save the plot. 
  if save is True:
    return PW.render(plotdir + 'fwhm_rate.pdf')
  else:
    return PW.render()

# =============================================================================
# ====================================== Poloidal or Toroidal Events by L - LPP
# =============================================================================

# Let's take a look at a plot of how FWHM depends on mode. 
def llppplot(mode, split=0., save=False):
  global pargs, plotdir
  # Set up the grid and 2D histogram based on probe position. 
  pos = getpos(**pargs)
  x, y, z, hargs = [ pos[key] for key in ('x', 'y', 'z', 'hargs') ]
  # Set up the window. 
  PW = plotWindow( ncols=2, nrows=2, **bep() )
  modename = 'Poloidal' if mode=='p' else 'Toroidal'
  title = notex(modename + ' Pc4 by Plasmapause Location: ' + tex('ImS') +
                 ' \\geq ' + str(thresh) + tex('mW/m^2') )
  rowlabels = ( notex('Odd\nHarmonic'), notex('Even\nHarmonic') )
  collabels = ( 'L - L_{PP} < ' + str(split), 'L - L_{PP} > ' + str(split) )
  PW.setParams(collabels=collabels, rowlabels=rowlabels, title=title, unitlabel='\\%')
  # Iterate over the filters. 
  for row, mname in enumerate( ('P1', 'P2') if mode=='p' else ('T1', 'T2') ):
    # Grab a histogram of the events, filtered by spectral width. 
    inside  = eventhist(hargs, mode=mname, llpp_lt=split)
    outside = eventhist(hargs, mode=mname, llpp_ge=split)

    print mname + ' overall inside  rate: ' + format(100*np.sum(inside)/np.sum(z), '.1f') + '%'
    print mname + ' overall outside rate: ' + format(100*np.sum(outside)/np.sum(z), '.1f') + '%'

    # Normalize by how long each region was sampled. 
    rates  = 100*zmask(inside)/z, 100*zmask(outside)/z
    # Add the mesh to the plot. 
    [ PW[row, i].setMesh(x, y, r) for i, r in enumerate(rates) ]
  # Show or save the plot. 
  if save is True:
    return PW.render(plotdir + 'llpp_rate.pdf')
  else:
    return PW.render()

# =============================================================================
# ========================================== Poloidal or Toroidal Events by Dst
# =============================================================================

# Let's take a look at a plot of how FWHM depends on mode. 
def dstplot(mode, split=0., save=False):
  global pargs, plotdir
  # Set up the grid and 2D histogram based on probe position. 
  pos = getpos(**pargs)
  x, y, z, hargs = [ pos[key] for key in ('x', 'y', 'z', 'hargs') ]
  # Set up the window. 
  PW = plotWindow( ncols=2, nrows=2, **bep() )
  modename = 'Poloidal' if mode=='p' else 'Toroidal'
  title = notex(modename + ' Pc4 by Dst: ' + tex('ImS') +
                 ' \\geq ' + str(thresh) + tex('mW/m^2') )
  rowlabels = ( notex('Odd\nHarmonic'), notex('Even\nHarmonic') )
  collabels = ( 'Dst < ' + str(split), 'Dst > ' + str(split) )
  PW.setParams(collabels=collabels, rowlabels=rowlabels, title=title, unitlabel='\\%')
  # Iterate over the filters. 
  for row, mname in enumerate( ('P1', 'P2') if mode=='p' else ('T1', 'T2') ):
    # Grab a histogram of the events, filtered by spectral width. 
    calm  = eventhist(hargs, mode=mname, dst_lt=split)
    storm = eventhist(hargs, mode=mname, dst_ge=split)

    print mname + ' overall calm  rate: ' + format(100*np.sum(calm)/np.sum(z), '.1f') + '%'
    print mname + ' overall storm rate: ' + format(100*np.sum(storm)/np.sum(z), '.1f') + '%'

    # Normalize by how long each region was sampled. 
    rates  = 100*zmask(inside)/z, 100*zmask(outside)/z
    # Add the mesh to the plot. 
    [ PW[row, i].setMesh(x, y, r) for i, r in enumerate(rates) ]
  # Show or save the plot. 
  if save is True:
    return PW.render(plotdir + 'llpp_rate.pdf')
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
  bigenough = np.nonzero(amp - thresh >= 0)[0]
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

def getpos(dl=0.5, dm=1, lmin=None, lmax=None):
  # The options for units are hours and days. 
  secondsperday = 86400.
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
  x, y = -l*np.sin(2*pi*m/24.), -l*np.cos(2*pi*m/24.)
  # Bin the position data into a 2D histogram, then scale it to days. 
  h = np.histogram2d(pos[:, 0], posm, **hargs)[0]
  z = 300*h/secondsperday
  # Return the position data. Total amount of usable time too. 
  return {'dates':dates, 'l':l, 'm':m, 'x':x, 'y':y, 'z':z, 'hargs':hargs}

# =============================================================================
# ================================================== Loading and Binning Events
# =============================================================================

# Returns a filtered list of events. 
def loadevents(mode=None, fwhm_ge=None, fwhm_lt=None, comp_ge=None, comp_lt=None, double=False, llpp_lt=None, llpp_ge=None):
  global thresh
  # Grab the contents of the event file as an array of strings. 
  events = g2a( read('events.txt') )
  # Filter for simultaneous events. 
  if double is True:
    events = g2a( line for line in events if 'BIG' in line or 'SMALL' in line )
  # Filter based on mode. 
  if mode is not None:
    events = g2a( line for line in events if mode in line )
  # Filter on amplitude. 
  if thresh > 0:
    amp = g2a( float( line.split()[10] ) for line in events )
    inew = np.nonzero(amp >= thresh)[0]
    events = events[inew]
  # Filter on compressional coupling (lower bound). 
  if comp_ge is not None:
    comp = g2a( float( line.split()[11] ) for line in events )
    inew = np.nonzero(comp >= comp_ge)[0]
    events = events[inew]
  # Filter on compressional coupling (upper bound). 
  if comp_lt is not None:
    comp = g2a( float( line.split()[11] ) for line in events )
    inew = np.nonzero(comp < comp_lt)[0]
    events = events[inew]
  # Filter on spectral width (lower bound). 
  if fwhm_ge is not None:
    fwhm = g2a( float( line.split()[9] ) for line in events )
    inew = np.nonzero(fwhm >= fwhm_ge)[0]
    events = events[inew]
  # Filter on compressional coupling (upper bound). 
  if fwhm_lt is not None:
    fwhm = g2a( float( line.split()[9] ) for line in events )
    inew = np.nonzero(fwhm < fwhm_lt)[0]
    events = events[inew]
  # Filter on position relative to the plasmapause (lower bound). 
  if llpp_ge is not None:
    lpp = g2a( float( line.split()[6] ) for line in events )
    lshell = g2a( float( line.split()[3] ) for line in events )
    inew = np.nonzero(lshell - lpp >= llpp_ge)[0]
    events = events[inew]
  # Filter on position relative to the plasmapause (upper bound). 
  if llpp_lt is not None:
    lpp = g2a( float( line.split()[6] ) for line in events )
    lshell = g2a( float( line.split()[3] ) for line in events )
    inew = np.nonzero(lshell - lpp < llpp_lt)[0]
    events = events[inew]
  # Return the remaining events. 
  return events

# Returns a histogram of filtered events. 
def eventhist(hargs, **kargs):
  # Days per event. Scales to days from a count of half-hour events. 
  dpe = 1800./86400.
  # Grab the events. 
  events = loadevents(**kargs)

  # If there are no events, return an empty histogram. 
  if events.size == 0:
    return dpe*np.histogram2d([], [], **hargs)[0]

  # Get the position from each. 
  pos = g2a( [ float(x) for x in line.split()[3:6] ] for line in events )
  lshell, mlat, mlt = pos[:, 0], pos[:, 1], pos[:, 2]
  # Assemble the positions of these events into a histogram, using the bins
  # defined by hargs. Scale from a count of half-hour events to days. 
  return dpe*np.histogram2d(lshell, mlt, **hargs)[0]

# Returns a histogram of double events: those that trigger the poloidal and 
# toroidal channels simultaneously. 
def doublehist(hargs, pmode=None, tmode=None, **kargs):
  # Days per event. Scales to days from a count of half-hour events. 
  dpe = 1800./86400.
  # Grab the poloidal and toroidal events. 
  pev = loadevents(mode=pmode, double=True, **kargs)
  tev = loadevents(mode=tmode, double=True, **kargs)
  # The number of events should be pretty small, so just find matches by brute
  # force. The first 26 characters of the line give probe, date, time. 
  doubles = g2a( p + '\n' + t for p in pev for t in tev if p[:26]==t[:26] )
  # If there are no events, return an empty histogram. 
  if doubles.size == 0:
    return dpe*np.histogram2d([], [], **hargs)[0]

  dates = g2a( d[3:13] for d in doubles )

  # Tally up the days. 
  print pmode, tmode
  print '\tnumber of events: ', dates.size
  print '\tnumber of dates:  ', g2a( set(dates) ).size

  # Get the position for each event, then return a histogram of those events,
  # scaled to units of days. 
  pos = g2a( [ float(x) for x in line.split()[3:6] ] for line in doubles )
  lshell, mlat, mlt = pos[:, 0], pos[:, 1], pos[:, 2]
  return dpe*np.histogram2d(lshell, mlt, **hargs)[0]

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
# ############################################################ Helper Functions
# #############################################################################

# Turns a generator expression into a list, if necessary, then turns the list
# into an array. 
def g2a(expr):
  return np.array( list(expr) )

# Make sure we don't throw anything infinite on the plot. 
def fmask(x):
  return masked_where(np.logical_not( np.isfinite(x) ), x)

# Also leave out any zeros. 
def zmask(x, thr=0):
  return masked_where(np.abs(x) <= thr, x)

# Bullseye params for the plotter. 
def bep():
  tls = ('$-8$', '', '$-4$', '', '$0$', '', '$+4$', '', '$+8$')
  return {'earth':'top', 'flipx':True, 'grid':True, 'square':True, 
          'xlabel': 'Y' + notex(' (R_E)'), 'xlims':(-8, 8),
          'xticks':np.mgrid[-8:8:9j], 'xticklabels':tls, 
          'ylabel': 'X' + notex(' (R_E)'), 'ylims':(-8, 8), 
          'yticks':np.mgrid[-8:8:9j], 'yticklabels':tls, 'ylabelpad':-2, 
          'colorbar':'pos' }

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


