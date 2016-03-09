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
  for date in sorted( os.listdir('/media/My Passport/rbsp/pkls/') )[:200]:

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

  # RBSP is near (ish) to the equator, so an odd-mode poloidal FLR should have
  # an electric field antinode. Threshold out anything without a Fourier
  # component of at least 0.25mV/m. 
  efft = ev.fft('ey')
  if np.max(efft) < 0.25:
    return False

  # Filter by frequency. Insist that the peak frequency is in the Pc4 band. 
  ipeak = np.argmax(efft)
  if not ev.ispc4()[ipeak]:
    return False

  # Filter by cross-spectral density. Make sure the peak cross-spectral density
  # for the poloidal mode is near the peak for the poloidal electric field (in
  # Fourier space) -- though they need not match exactly! This is to filter out
  # the surprisingly-common events where Bx is dominated by an oscillation
  # clearly different from the one dominant in Ey... perhaps due to a
  # superposition of first and second harmonics? 
  if np.abs( ev.frq()[ np.argmax( ev.csd('p') ) ] - ev.frq()[ipeak] ) > 5:
    return False

  # Filter by coherence. This is how we ensure the phase lag is meaningful. 
  if not ev.iscoh('p')[ipeak]:
   return False

  # Filter by harmonic. We only want odd modes. 
  if not ev.harm('p')[ipeak] % 2:
    return False

  # If an event is found, return a line describing it. 
  return ( col(ev.probe) + col(ev.date) + col(ev.time) +
           col(ev.frq()[ipeak], unit='mHz') +
           col(np.max(efft), unit='mV/m') +
#           col( ev.coh('p')[ipeak], unit='COH' ) +
           col(ev.lag('p')[ipeak], unit='*') + 
           col(ev.avg('lshell'), unit='L') + 
#           col(ev.avg('mlt'), unit='MLT') + 
#           col(ev.avg('mlat'), unit='*MLAT') + 
           col(ev.lpp, unit='LPP') )

# #############################################################################
# ############################################################# Plotting Events
# #############################################################################

def plot(ev, save=False):
  global plotdir
  # Create plot window to hold waveforms and spectra. 
  PW = plotWindow(nrows=3, ncols=2)
  # Index the poloidal, toroidal, and field-aligned modes. 
  modes = ('p', 't', 'z')
  # Plot waveforms. 
  PW[:, 0].setParams( **ev.coords('waveform', cramped=True) )
  [ PW[i, 0].setLine(ev.get('B' + m), 'r') for i, m in enumerate(modes) ]
  [ PW[i, 0].setLine(ev.get('E' + m), 'b') for i, m in enumerate(modes) ]
  # Plot spectral properties. 
  PW[:, 1].setParams( **ev.coords('coherence', cramped=True) )
  [ PW[i, 1].setLine( ev.csd(m), 'g' ) for i, m in enumerate(modes) ]
  [ PW[i, 1].setLine( ev.coh(m), 'k' ) for i, m in enumerate(modes) ]
  [ PW[i, 1].setLine( 0.1 + 0.8*ev.isodd(m), 'm' ) for i, m in enumerate(modes) ]
  # Put title and labels on the plot. 
  rowlabels = ( notex('Poloidal'), notex('Toroidal'), notex('Parallel') )
  collabels = ( notex('B (Red) ; E (Blue)'), 
                notex('Coh. (Black) ; CSD (Green) ; Par. (Violet)') )
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


