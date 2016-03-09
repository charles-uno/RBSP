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

  # Nuke any previous event list to avoid double counting. 
  if os.path.exists('oddevents.txt'):
    os.remove('oddevents.txt')

  # There's one pickle directory for each date we're supposed to look at. 
  for date in sorted( os.listdir('/media/My Passport/rbsp/pkls/') )[:400]:

    # Check each date for both probes. Thirty minute chunks. 
    [ checkdate(probe, date, mpc=30) for probe in ('a', 'b') ]

  return

# #############################################################################
# ######################################################## Searching for Events
# #############################################################################

# =============================================================================
# ====================================================== Check a Day for Events
# =============================================================================

# The day is broken into chunks (mpc is minutes per chunk). Each chunk is run
# through a few different filters to seek out odd-harmonic poloidal Pc4 waves. 
def checkdate(probe, date, mpc=30):

  print date, probe

  # Load the day's data into a day object. 
  today = day(probe=probe, date=date)

  # Break the day into chunks and look at each chunk. 
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

  # This filter is a bit trickier than Lei's, because we're going a level deeper. We need not only the magnetic field, but also the electric field and how they relate to one another. 

  # Let's look for the maximum component as judged by the product of the poloidal electric and magnetic fields. This is something like the cross-spectral density, though there's no averaging over windows. 
  bfft, efft = ev.fft('bx'), ev.fft('ey')
  csd = np.abs(bfft*efft)
  imax = np.argmax(csd)

  # Threshold at 0.2 nT mV/m. This is pretty small. 
  if csd[imax] < 0.2:
    return False

  # Filter frequencies phenomenologically, rather than per the IAGA designations. Pgs are sometimes as slow as 5mHz, and top out at 17mHz (anything faster than that is probably a third harmonic). 
  if not ( 5 < ev.frq()[imax] < 17 ):
    return False

  # The coherence needs to be high so that we can compute a meaningful phase offset. 
  if not ev.coh('p')[imax] > 0.9:
    return False

  # The phase must fall between 60 and 120 degrees, with the correct sign. 
  if not ev.isodd('p')[imax]:
    return False

  # Let's also note the comparison between phase computed from CSD with phase computed from FFT components. 
  print '\tCoherence = ', ev.coh('p')[imax]
  print '\tMLAT      = ', ev.avg('mlat')
  print '\tPhase     = ', np.angle( ev.fft('ey')*np.conj( ev.fft('bx') ) , deg=True)[imax], 'degrees'
  print '\tCSD       = ', np.abs( ev.fft('ey')*ev.fft('bx') )[imax], ' nT mV/m'

  # If an event is found, return a line describing it. 
  return ( col(ev.probe) + col(ev.date) + col(ev.time) +
           col(ev.frq()[imax], unit='mHz') +
           col(np.max(csd), unit='nTmV/m') +
#           col( ev.coh('p')[imax], unit='COH' ) +
           col(ev.lag('p')[imax], unit='deg') + 
#           col(ev.avg('lshell'), unit='L') + 
#           col(ev.avg('mlt'), unit='MLT') + 
#           col(ev.avg('mlat'), unit='*MLAT') + 
#           col(ev.lpp, unit='LPP') 
           col(ev.coh('p')[imax], unit='Coh')
         )

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
  [ PW[i, 1].setLine( ev.coh(m), 'g' ) for i, m in enumerate(modes) ]
  bfft = [ np.abs( ev.fft('B' + m) ) for m in modes ]
  efft = [ np.abs( ev.fft('E' + m) ) for m in modes ]
  cfft = [ np.abs( ev.fft('B' + m)*ev.fft('E' + m) ) for m in modes ]
  [ PW[i, 1].setLine(c/np.max(c), 'm') for i, c in enumerate(cfft) ]
  [ PW[i, 1].setLine(b/np.max(b), 'r:') for i, b in enumerate(bfft) ]
  [ PW[i, 1].setLine(e/np.max(e), 'b:') for i, e in enumerate(efft) ]
  # Parity clogs up the plot too much. 
#  [ PW[i, 1].setLine( 0.9*ev.isodd(m), 'orange' ) for i, m in enumerate(modes) ]
  # Put title and labels on the plot. 
  rowlabels = ( notex('Poloidal'), notex('Toroidal'), notex('Parallel') )
  collabels = ( notex('B (Red) ; E (Blue)'), 
                notex('\\overset{\\sim}{E}\\,\\overset{\\sim}{B} (Magenta) ;' +
                      ' Coherence (Green)') )
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


