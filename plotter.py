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

def main():

#  # Plot the location of the usable data. 
#  return posplot(save='-i' in argv)

  # Plot the location of events... ?
  return eventplot(save='-i' in argv)

  return

# #############################################################################
# ########################################################## Position Histogram
# #############################################################################

def posplot(save=False, unit='days'):
  global plotdir
  # Grab the position data. 
  unit = 'days'
  pos = getpos(unit=unit)
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
# ############################################################ Events Histogram
# #############################################################################

# Make sure we don't throw anything infinite on the plot. 
def fmask(x):
  return masked_where(np.logical_not( np.isfinite(x) ), x)

# Also leave out any zeros. 
def zmask(x, thresh=0):
  return masked_where(np.abs(x)<=thresh, x)

def eventplot(save=False):
  global plotdir
  # Set up the grid and 2D histogram based on probe position. 
  pos = getpos(dl=1)
  x, y, z, hargs = [ pos[key] for key in ('x', 'y', 'z', 'hargs') ]
  d0, d1 = pos['dates']
  # Create a plot window to show different subsets of the events. 
  mfilt, hfilt = ('P', 'T'), ('1', '2')
  PW = plotWindow( ncols=2, nrows=2, colorbar='pos', **bep() )
  # Iterate over the filters. 
  for row, hf in enumerate(hfilt):
    for col, mf in enumerate(mfilt):
      # Grab a histogram of the events, filtered by mode and harmonic. 
      events = getevents(hargs, filt=mf+hf)
      # Normalize by how long each region was sampled. 
      rate = 100*zmask(events)/zmask( z, thresh=0.05*np.max(z) )
      # Add the mesh to the plot. 
      PW[row, col].setMesh(x, y, rate)
  # Title and labels. 
  title = notex('Pc4 Occurrence Rate by Mode, ' + d0 + ' to ' + d1)
  collabels = ( notex('Poloidal'), notex('Toroidal') )
  rowlabels = ( notex('Odd\nHarmonic'), notex('Even\nHarmonic') )
  PW.setParams(collabels=collabels, rowlabels=rowlabels, title=title, unitlabel='\\%')
  # Show or save the plot. 
  if save is True:
    return PW.render(plotdir + 'rate.pdf')
  else:
    return PW.render()

# #############################################################################
# ################################################################## Data Input
# #############################################################################

# =============================================================================
# ==================================================================== Position
# =============================================================================

def getpos(dl=0.5, dm=1, lmin=None, lmax=None, mmin=None, mmax=None, 
           unit='days'):
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
  if mmin is None:
    mmin = np.floor( np.min( pos[:, 1] ) )
  if mmax is None:
    mmax = np.ceil( np.max( pos[:, 1] ) )
  # Number of bins in each direction. 
  lbins, mbins = int( (lmax - lmin)/dl ) + 1, int( (mmax - mmin)/dm ) + 1
  # Keyword arguments for the histogram2d call. 
  hargs = { 'range':( (lmin, lmax), (mmin, mmax) ), 'bins':(lbins-1, mbins-1) }
  # Bin bounds in terms of L and MLT. 
  l, m = np.mgrid[lmin:lmax:lbins*1j, mmin - dm/2.:mmax - dm/2.:mbins*1j]
  # Map to GSE coordinates. Put midnight at the bottom. 
  x, y = l*np.sin(2*pi*m/24.), -l*np.cos(2*pi*m/24.)
  # Bin the position data into a 2D histogram, then scale it to days. 
  h = np.histogram2d(pos[:, 0], pos[:, 1], **hargs)[0]
  z = 300*h/secondsper
  # Return the position data. Total amount of usable time too. 
  return {'dates':dates, 'l':l, 'm':m, 'x':x, 'y':y, 'z':z, 'hargs':hargs}

# =============================================================================
# ====================================================================== Events
# =============================================================================

def getevents(hargs, unit='days', filt=''):
  secondsper = 86400. if unit=='days' else 1440.
  evlines = [ line for line in read('events.txt') if filt in line ]
  pos = np.array( [ [ float(x) for x in l.split()[3:6] ] for l in evlines ] )
  return 1800*np.histogram2d(pos[:, 0], pos[:, 1], **hargs)[0]/secondsper

# #############################################################################
# ############################################################ Helper Functions
# #############################################################################

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


