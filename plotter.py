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

# #############################################################################
# ######################################################################## Main
# #############################################################################

# A timestamped directory, in case we want to save any plots.  
plotdir = '/home/user1/mceachern/Desktop/plots/' + now() + '/'

def main():

  # Plot the location of the usable data. 
  return posplot(save='-i' in argv)

#  # Plot the location of events... ?
#  return rateplot(save='-i' in argv)

  return

# #############################################################################
# ##################################################### Probe Position Handling
# #############################################################################

# =============================================================================
# ========================================================= Parse Position Data
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
  # Bin bounds in terms of L and MLT. 
  l, m = np.mgrid[lmin:lmax:lbins*1j, mmin - dm/2.:mmax - dm/2.:mbins*1j]
  # Map to GSE coordinates. Put MLT=0 at the bottom. 
  x, y = l*np.sin(2*pi*m/24.), -l*np.cos(2*pi*m/24.)
  # Bin the position data into a 2D histogram, then scale it to days. 
  hist = np.histogram2d( pos[:, 0], pos[:, 1], bins=(lbins-1, mbins-1), 
                         range=[ (lmin, lmax), (mmin, mmax) ] )[0]
  z = 300*hist/secondsper
  # Return the position data. Total amount of usable time too. 
  return {'dates':dates, 'x':x, 'y':y, 'z':z}

# =============================================================================
# ===================================================== Plot Position Histogram
# =============================================================================

def posplot(save=False, unit='days'):
  global plotdir
  # Grab the position data. 
  unit = 'days'
  pos = getpos(unit=unit)
  x, y, z = [ pos[key] for key in ('x', 'y', 'z') ]
  date0, date1 = pos['dates']
  dt = np.sum(z)
  # Create the plot window and scale the axes properly. 
  PW = plotWindow(square=True, colorbar='pos')
  lms, tks = (-8, 8), np.mgrid[-8:8:9j]
  tls = [ '$' + format(t, '+.0f') + '$' if t%4==0 else '' for t in tks  ]
  tls[4] = '$0$'
  PW.setParams(xlims=lms, ylims=lms, xticks=tks, yticks=tks, 
               xticklabels=tls, yticklabels=tls)
  # Set the title and labels. 
  xlbl, ylbl = 'Y' + notex(' (R_E)'), 'X' + notex(' (R_E)')
  ulabel = notex(unit)
  title = notex('Usable Data ' + date0 + ' to ' + date1 + ' (' + 
                format(dt, '.0f') + ' ' + unit + ' total)')
  PW.setParams(title=title, unitlabel=ulabel, xlabel=xlbl, ylabel=ylbl)
  # Draw Earth. Mark the grid. 
  PW.setParams(earth='top', grid=True)
  # Add the data to the plot. 
  PW.setMesh(x, y, z)
  # Show the plot, or save it as an image. 
  if save is True:
    return PW.render(plotdir + 'pos.pdf')
  else:
    return PW.render()

# #############################################################################
# ################################################### Event Occurrence Handling
# #############################################################################

# =============================================================================
# ============================================================ Parse Event Data
# =============================================================================

def getevents():

  eventlines = [ line for line in read('events.txt') ]




# =============================================================================
# ===================================================== Plot Position Histogram
# =============================================================================

def rateplot(save=False):
  global plotdir
  unit = 'days'
  secondsper = 86400. if unit=='days' else 1440.
  # The orbit of both RBSP paths has been broken into five-minute chunks. Grab
  # the position of each chunk which gives good data for E dot B = 0. 
  poslines = [ line for line in read('pos.txt') if line.endswith('ok') ]
  # Tally up the total amount of usable data at five minutes per line, in days.
  utime = 300*len(poslines)/secondsper
  # Get the date range. 
  date0, date1 = poslines[0].split()[1], poslines[-1].split()[1]
  # Parse the positions into an array of floats. 
  pos = np.array( [ [ float(x) for x in p.split()[3:6] ] for p in poslines ] )
  # Figure out appropriate bin sizes for the position data. 
  dl, dm = 0.5, 1
  lmin, lmax = np.floor( np.min( pos[:, 0] ) ), np.ceil( np.max( pos[:, 0] ) )
  mmin, mmax = np.floor( np.min( pos[:, 1] ) ), np.ceil( np.max( pos[:, 1] ) )
  lbins, mbins = int( (lmax - lmin)/dl ) + 1, int( (mmax - mmin)/dm ) + 1
  # Bin the position data into a 2D histogram, then scale it to days. 
  hist = np.histogram2d( pos[:, 0], pos[:, 1], bins=(lbins-1, mbins-1), 
                         range=[ (lmin, lmax), (mmin, mmax) ] )[0]
  z = 300*hist/secondsper
  # Bin bounds in terms of L and MLT. 
  l, m = np.mgrid[lmin:lmax:lbins*1j, mmin - dm/2.:mmax - dm/2.:mbins*1j]
  # Map to GSE coordinates. Put MLT=0 at the bottom. 
  x, y = l*np.sin(2*pi*m/24.), -l*np.cos(2*pi*m/24.)
  # Create the plot window and scale the axes properly. 
  PW = plotWindow(square=True, colorbar='pos')
  lms, tks = (-8, 8), np.mgrid[-8:8:9j]
  tls = [ '$' + format(t, '+.0f') + '$' if t%4==0 else '' for t in tks  ]
  tls[4] = '$0$'
  PW.setParams(xlims=lms, ylims=lms, xticks=tks, yticks=tks, 
               xticklabels=tls, yticklabels=tls)
  # Set the title and labels. 
  xlbl, ylbl = 'Y' + notex(' (R_E)'), 'X' + notex(' (R_E)')
  ulabel = notex(unit)
  title = notex('Usable Data ' + date0 + ' to ' + date1 + ' (' + 
                format(utime, '.0f') + ' ' + unit + ' total)')
  PW.setParams(title=title, unitlabel=ulabel, xlabel=xlbl, ylabel=ylbl)
  # Draw the grid. 
  [ PW.setLine(x[i, :], y[i, :], 'k') for i in range( x.shape[0] ) ]
  [ PW.setLine(x[:, j], y[:, j], 'k') for j in range( x.shape[1] ) ]
  # Draw Earth. This is a bit kludgey. 
  ax = ax = PW.cells.flatten()[0].ax
  ax.add_artist( Wedge( (0, 0), 1, 0, 180, fc='w' ) )
  ax.add_artist( Wedge( (0, 0), 1, 180, 360, fc='k' ) )

  # Draw an asterisk in the middle of the largest bin. 
  i, j = np.unravel_index(np.argmax(z), z.shape)

  print 'loc = ', i, j

  xmid = 0.5*( x[:-1, :] + x[1:, :]  )
  ymid = 0.5*( y[:, :-1] + y[:, 1:]  )

  print 'x, y = ', xmid[i, j], ymid[i, j]

  kargs = {'x':xmid[i, j], 'y':ymid[i, j], 'horizontalalignment':'center', 'verticalalignment':'center', 'fontsize':15}

  ax.text(s='$*$', **kargs)





  # Add the data to the plot. 
  PW.setMesh(x, y, z)
  # Show the plot, or save it as an image. 
  if save is True:
    return PW.render(plotdir + 'rate.pdf')
  else:
    return PW.render()

# #############################################################################
# ############################################################ Helper Functions
# #############################################################################

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


