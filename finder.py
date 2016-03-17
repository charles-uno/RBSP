#!/usr/bin/env python

# Charles McEachern

# Spring 2016

# #############################################################################
# #################################################################### Synopsis
# #############################################################################

# Lei's got some filters in his data which are designed to screen out
# fundamental mode Pc4 pulsations. We run through the same date with different
# filters, specifically to find those fundamental mode poloidal Pc4 events. 

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



  # Let's plot where the data is! 
  return bullseye(save='-i' in argv)






  '''
  # Debug with a known super-nice event. 
  date = '2012-10-23'
  probe = 'a'
  t = timeint(time='22:00:00')
  today = day(probe=probe, date=date)
  mpc = 30
  ev = today.getslice(t, duration=60*mpc)
  evline = checkevent(ev)
  print evline
  return plot(ev)
  '''

#  # Nuke any previous event list to avoid double counting. 
#  if os.path.exists('oddevents.txt'):
#    os.remove('oddevents.txt')

  # There's one pickle directory for each date we're supposed to look at. 
  for date in sorted( os.listdir('/media/My Passport/rbsp/pkls/') ):

    print date

#    # Figure where the probes spent their time. Five minute chunks. 
#    [ trackpos(probe, date, mpc=5) for probe in ('a', 'b') ]

#    # Check each date for both probes. Thirty minute chunks. 
#    [ checkdate(probe, date, mpc=30) for probe in ('a', 'b') ]

  return

# #############################################################################
# #################################################### Tallying RBSP's Position
# #############################################################################

def trackpos(probe, date, mpc=5):

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

    evline = ( ev.probe + '\t' + ev.date + '\t' + ev.time + '\t' + format(ev.avg('lshell'), '.1f') + '\t' + format(ev.avg('mlt'), '.1f') + '\t' + format(ev.avg('mlat'), '.1f') )

    append(evline + '\t' + ( 'ok' if ev.isok() else 'X' ), 'pos.txt')

  return


#    evline = checkevent(ev)
#    if evline:
#      print append(evline, 'oddevents.txt')
#      plot(ev, save=True)




# =============================================================================
# ======================================================= Plotting the Position
# =============================================================================

def bullseye(save=False, unit='days'):
  global plotdir

#  unit = 'hours'

  secondsper = 86400. if unit=='days' else 1440.


  # The orbit of both RBSP paths has been broken into five-minute chunks. Grab
  # the position of each chunk which gives good data for E dot B = 0. 
  poslines = [ line for line in read('pos.txt') if line.endswith('ok') ]

  # Tally up the total amount of usable data at five minutes per line, in days.
  utime = 300*len(poslines)/secondsper

  # Get the date range. 
  date0, date1 = poslines[0].split()[1], poslines[-1].split()[1]
#  start = {'date':poslines[0].split()[1], 'time':poslines[0].split()[2]}
#  finish = { 'date':poslines[-1].split()[1], 'time':poslines[-1].split()[2] }
#  pdays = 2*( timeint(**finish) - timeint(**start) )/86400.
#  title = notex(format(udays, '.1f') + ' usable probe-days over ' + format(pdays, '.1f') + ' probe-days')

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
  title = notex('Usable Data ' + date0 + ' to ' + date1 + ' (' + format(utime, '.0f') + ' ' + unit + ' total)')
  PW.setParams(title=title, unitlabel=ulabel, xlabel=xlbl, ylabel=ylbl)

  # Draw the grid. 
  [ PW.setLine(x[i, :], y[i, :], 'k') for i in range( 1, x.shape[0] ) ]
  [ PW.setLine(x[:, j], y[:, j], 'k') for j in range( x.shape[1] ) ]

  # Draw Earth. This is a bit kludgey. 
  ax = ax = PW.cells.flatten()[0].ax
  ax.add_artist( Wedge( (0, 0), 1, 0, 180, fc='w' ) )
  ax.add_artist( Wedge( (0, 0), 1, 180, 360, fc='k' ) )

  # Add the data to the plot. 
  PW.setMesh(x, y, z)

  # Show the plot, or save it as an image. 
  if save is True:
    if not os.path.exists(plotdir):
      os.makedirs(plotdir)
    return PW.render(plotdir + 'pos.pdf')
  else:
    return PW.render()

# #############################################################################
# ######################################################## Searching for Events
# #############################################################################

# =============================================================================
# ====================================================== Check a Day for Events
# =============================================================================

# The day is broken into chunks (mpc is minutes per chunk). Each chunk is run
# through a few different filters to seek out odd-harmonic poloidal Pc4 waves. 
def checkdate(probe, date, mpc=30):

  # Load the day's data into a day object. 
  today = day(probe=probe, date=date)

  # Break the day into chunks and look for events in each chunk. 
  for t in range(0, 86400, 60*mpc):
    # Whenever an event is found, save it to file. 
    ev = today.getslice(t, duration=60*mpc)
    evline = checkevent(ev)
    if evline:
      print append(evline, 'oddevents.txt')
      plot(ev, save=True)

  return

# =============================================================================
# ====================================================== Check a Chunk of a Day
# =============================================================================

# Check if the given event -- which is a slice of a day's worth of data -- has
# any odd-mode poloidal Pc4 events. 
def checkevent(ev):

  # If the event's data is no good, obviously there can be no event. 
  if not ev.isok():
    return False

  # This filter is a bit trickier than Lei's, because we're going a level
  # deeper. We need not only the magnetic field, but also the electric field
  # and how they relate to one another. 

  # Get the (complex) poloidal Poynting flux in Fourier space. Find the point
  # at which its imaginary (standing wave) component is strongest. The value
  # is then scaled by L^3/mu0 (to map it to Poynting flux at the atmosphere). 
  sfft = ev.sfft('p')
  imax = np.argmax( np.abs( np.imag(sfft) ) )

  # Filter out anything that's not in the Pc4 frequency range. 
  if not ev.ispc4()[imax]:
    return False

  # Insist that all events be highly coherent so that phase offset is meaningful. 
  if not ev.iscoh('p')[imax]:
    return False

  # The magnetic latitude and the EB* phase must have opposite signs. This is
  # how we distinguish odd from even harmonics. 
  if not ev.isodd('p')[imax]:
    return False

  # Impose a cutoff on magnitude after mapping the Poynting flux to Earth. 
  if np.abs( np.imag( sfft[imax] ) ) < 0.02:
    return False

  # Notably, we are not filtering by spectral width. This should allow us to see giant pulsations in proportion. 

  # We're also not filtering by how much of the cross-spectral density is in the standing mode... this might have to change. 

  # If an event is found, return a line describing it. 
  return col(ev.probe) + col(ev.date) + col(ev.time)
#           col(ev.frq()[imax], unit='mHz')
#           col(np.max(csd), unit='nTmV/m')
#           col(ev.lag('p')[imax], unit='deg') 
#           col(ev.avg('lshell'), unit='L') 
#           col(ev.avg('mlt'), unit='MLT')
#           col(ev.avg('mlat'), unit='*MLAT')
#           col(ev.lpp, unit='LPP') 
#           col(ev.coh('p')[imax], unit='Coh')

# #############################################################################
# ############################################################# Plotting Events
# #############################################################################

def plot(ev, save=False):
  global plotdir
  # Create plot window to hold waveforms and spectra. 
  PW = plotWindow(nrows=3, ncols=2)
  # Index the poloidal, toroidal, and field-aligned modes. 
  modes = ('p', 't', 'z')
  # Plot waveforms as a function of time. 
  PW[:, 0].setParams( **ev.coords('waveform', cramped=True) )
  [ PW[i, 0].setLine(ev.get('B' + m), 'r') for i, m in enumerate(modes) ]
  [ PW[i, 0].setLine(ev.get('E' + m), 'b') for i, m in enumerate(modes) ]
  # Plot Fourier component magnitudes as a function of frequency. 
  PW[:, 1].setParams( **ev.coords('spectra', cramped=True) )
  # Dotted lines for the magnitudes of field transforms. 
  bfft = [ np.abs( ev.fft('b' + m) ) for m in modes ]
  efft = [ np.abs( ev.fft('e' + m) ) for m in modes ]
  [ PW[i, 1].setLine(b/np.max(b), 'r:') for i, b in enumerate(bfft) ]
  [ PW[i, 1].setLine(e/np.max(e), 'b:') for i, e in enumerate(efft) ]
  # Also -- more interestingly -- plot the real and imaginary components of the
  # Fourier-space Poynting flux. 
  sfft = [ ev.sfft(m)/np.max( np.abs( ev.sfft(m) ) ) for m in modes ]
  sfftr = [ np.abs( np.real(s) ) for s in sfft ]
  sffti = [ np.abs( np.imag(s) ) for s in sfft ]
  [ PW[i, 1].setLine(s, 'm') for i, s in enumerate(sffti) ]
  [ PW[i, 1].setLine(s, 'g') for i, s in enumerate(sfftr) ]
  # Put title and labels on the plot. 
  rowlabels = ( notex('Poloidal'), notex('Toroidal'), notex('Parallel') )
  collabels = ( notex('B (Red) ; E (Blue)'), 
                tex('imag') + tex('SFFT') + notex(' (Magenta) ; ') + 
                tex('real') + tex('SFFT') + notex(' (Green)') )
  PW.setParams(collabels=collabels, sidelabel=ev.label(), title=ev.descr('p'), 
               rowlabels=rowlabels)
  # Show the plot, or save it as an image. 
  if save is True:
    if not os.path.exists(plotdir):
      os.makedirs(plotdir)
    return PW.render(plotdir + ev.name + '.png')
  else:
    return PW.render()

# #############################################################################
# ############################################################ Helper Functions
# #############################################################################

# Space out a nice right-justified column. 
def col(x, width=12, unit=''):
  if isinstance(x, float):
    return format(x, str(width - len(unit) - 1) + '.1f') + unit + ' '
  else:
    return str(x).rjust(width - len(unit) - 1) + unit + ' '

# #############################################################################
# ########################################################### For Importability
# #############################################################################

if __name__=='__main__':
  main()


