#!/usr/bin/env python

# Charles McEachern

# Spring 2016

# #############################################################################
# #################################################################### Synopsis
# #############################################################################

# This script makes plots of RBSP data. 

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

#  # What dates do we have data for? 
#  dates = sorted( os.listdir('/media/My Passport/rbsp/pkls/') )

#  # If we're saving our data, nuke the previous list to avoid double-counting. 
#  if '-i' in argv and os.path.exists('pos.txt'):
#    print 'Removing old position listing'
#    os.remove('pos.txt')

#  # Tally the probes' positions. 
#  for date in dates:
#    print date
#    [ trackpos(probe, date, mpc=5) for probe in ('a', 'b') ]

#  # If we're saving our data, nuke the previous list to avoid double-counting. 
#  if '-i' in argv and os.path.exists('events.txt'):
#    print 'Removing old event listing'
#    os.remove('events.txt')

#  # Search for events. Do the days in random order, for easier debugging. We
#  # can just look at the first event we find. 
#  for date in np.random.permutation(dates):
#    print date
#    # Check both probes. Thirty minute chunks. 
#    [ checkdate(probe, date, mpc=30) for probe in ('a', 'b') ]
  return


# #############################################################################
# ###################################################### Tabulate RBSP Position
# #############################################################################

def trackpos(probe, date, mpc=5):
  # This takes forever to run... 
  print 'THERE IS NO REASON YOU SHOULD BE CALLING THIS AGAIN. '
  exit()
  # Load the day's data into a day object. 
  today = day(probe=probe, date=date)
  # If the event is no good, bail. 
  if today.garbage:
    return
  # Scroll through the day a few minutes at a time. 
  for t in range(0, 86400, 60*mpc):
    # Grab a slice of the day. Print location and if the data is OK. 
    ev = today.getslice(t, duration=60*mpc)
    lshell, mlt, mlat = [ ev.avg(name) for name in ('lshell', 'mlt', 'mlat') ]
    evline = ( ev.probe + '\t' + ev.date + '\t' + ev.time + '\t' +
               format(ev.avg('lshell'), '.1f') + '\t' +
               format(ev.avg('mlt'), '.1f') + '\t' + 
               format(ev.avg('mlat'), '.1f') )
    append(evline + '\t' + ( 'ok' if ev.isok() else 'X' ), 'pos.txt')
  return

# #############################################################################
# ############################################################# Tabulate Events
# #############################################################################

# =============================================================================
# ===================================================== Search a Day for Events
# =============================================================================

# The day is broken into chunks (mpc is minutes per chunk). Each chunk is run
# through a few different filters to seek out odd-harmonic poloidal Pc4 waves. 
def checkdate(probe, date, mpc=30):
  # Load the day's data into a day object. 
  today = day(probe=probe, date=date)
  # Iterate over each chunk of the day. 
  for t in range(0, 86400, 60*mpc):
#    print '\t' + timestr(t)[1]
    # Check for poloidal and toroidal events independently. 
    ev = today.getslice(t, duration=60*mpc)
    evdicts = [ ev.standing(m, pc4=True, thresh=0.01) for m in ('p', 't') ]
    # If there's anything to save, do so, then plot it. 
    if keepevent(evdicts):
      plotevent(ev, save='-i' in argv)
      # If we're debugging, stop after a single plot. 
      if '-i' not in argv:
        exit()
  return

# =============================================================================
# ============================================================== Store an Event
# =============================================================================

# Assemble a dictionary about the event into a one-line summary. 
def evline(d):
  pdt = col( d['probe'] ) + col( d['date'] ) + col( d['time'] )
  pos = col( d['lshell'] ) + col( d['mlt'] ) + col( d['mlat'] )
  lpp = col( d['lpp'] )
  mh = col( d['mode'].upper() + str( d['harm'] ) )
  fdf = col( d['f'] ) + col( 2.355*d['df'] )
  mag = col( d['s'] )
  comp = col( d['comp'] )
  return pdt + pos + lpp + mh + fdf + mag + comp

# Write out the event to file. If there are two simultaneous events, indicate
# which is larger. 
def keepevent(evdicts):
  # Filter out non-events.
  evds = [ d for d in evdicts if d ]
  # If there are no events, bail. 
  if len(evds)==0:
    return 0
  # If there's one event, save it. 
  elif len(evds)==1:
    text = evline( evds[0] )
  # If there are two events, indicate which is larger. 
  elif len(evds)==2:
    ip = 0 if evds[0]['s'] > evds[1]['s'] else 1
    text = ( evline( evds[ip] ) + col('BIG') + '\n' + 
             evline( evds[1-ip] ) + col('SMALL') )
  # If we're storing the data, do so. In either case, print it. 
  if '-i' in argv:
    print append(text, 'events.txt')
  else:
    print text
  # Return an indication of how many events. 
  return len(evds)

# #############################################################################
# ############################################################### Plot an Event
# #############################################################################

def evtitle(d):
  mh = d['mode'].upper() + str( d['harm'] )
  frq = 'f\\!=\\!' + fmt(d['f'], digs=2) + tex('mHz')
  fwhm = '\\delta\\!f\\!=\\!' + fmt(2.355*d['df'], digs=2) + tex('mHz')
  mag = tex('imag') + tex('L3S') + '\\!=\\!' + fmt(d['s'], digs=2) + tex('mW/m^2')
  comp = tex('real') + tex( 'BB' + d['mode'] ) + '\\!=\\!' + fmt(d['comp'], digs=2)
  return '\\qquad{}'.join( (mh, frq, fwhm, mag, comp) )

def plotevent(ev, save=False):
  global plotdir
  # Create plot window to hold waveforms and spectra. 
  PW = plotWindow(nrows=3, ncols=2)
  # Index the poloidal, toroidal, and field-aligned modes. 
  modes = ('p', 't', 'z')
  # Plot waveforms as a function of time. 
  PW[:, 0].setParams( **ev.coords('waveform', cramped=True) )
  [ PW[i, 0].setLine(ev.get('B' + m), 'r') for i, m in enumerate(modes) ]
  [ PW[i, 0].setLine(ev.get('E' + m), 'b') for i, m in enumerate(modes) ]
  # Grab real and imaginary Fourier-domain Poynting flux, scaled by L^3. 
  sift = [ np.abs( np.imag( ev.sfft(m) ) ) for m in modes ]
  srft = [ np.abs( np.real( ev.sfft(m) ) ) for m in modes ]
  # Compute poloidal and/or toroidal standing waves. 
  stand = [ ev.standing(m, pc4=True, thresh=0.01) for m in ('p', 't') ]
  # Scale the Poynting flux by the largest value or the largest fit. 
  standmax = max( st['s'] for st in stand if st is not None )
  snorm = max(np.max(sift), np.max(srft), standmax)
  # Plot the real and imaginary spectral components. 
  PW[:, 1].setParams( **ev.coords('spectra', cramped=True) )
  [ PW[i, 1].setLine(s/snorm, 'm') for i, s in enumerate(sift) ]
  [ PW[i, 1].setLine(s/snorm, 'g') for i, s in enumerate(srft) ]
  # Plot the Gaussian fit of the imaginary spectral component. 
  f = np.linspace(0, 50, 1000)
  for i, st in enumerate(stand):
    args = (0, 0, 1) if st is None else ( st['s'], st['f'], st['df'] )
    PW[i, 1].setLine(f, ev.gauss(f, *args)/snorm, 'k--')
  # Plot title and labels. 
  rowlabels = ( notex('Poloidal'), notex('Toroidal'), notex('Parallel') )
  collabels = ( notex('B (Red) ; E (Blue)'), 
                tex('imag') + tex('L3S') + notex(' (Magenta) ; ') + 
                tex('real') + tex('L3S') + notex(' (Green)') )
  PW.setParams(collabels=collabels, title=ev.label(), rowlabels=rowlabels)
  # Information about the wave(s) goes in the side label. 
  tlist = [ evtitle(st) for st in stand if st is not None ]
  PW.setParams( sidelabel='$\n$'.join(tlist) )
  # Show the plot, or save it as an image. 
  if save is True:
    return PW.render(plotdir + ev.name + '.png')
  else:
    return PW.render()

# #############################################################################
# ######################################################## Plot RBSP's Position
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
# ============================================================ Parse Event Data
# =============================================================================

def getevents():

  eventlines = [ line for line in read('events.txt') ]










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




# =============================================================================
# ===================================================== Plot Position Histogram
# =============================================================================








# #############################################################################
# ############################################################# Plot Event Rate
# #############################################################################

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


