#!/usr/bin/env python

# Charles McEachern

# Spring 2016

# #############################################################################
# #################################################################### Synopsis
# #############################################################################

# Lei's got some filters in his data which are designed to screen out
# fundamental mode Pc4 pulsations. Can we come up with our own list of events? 
# Specifically, a list of fundamental mode Pc4 pulsations? 

# #############################################################################
# ######################################################### Load Python Modules
# #############################################################################

from day import *

from plotmod import *


# #############################################################################
# ######################################################################## Main
# #############################################################################


def main():

  for date in sorted( os.listdir('/media/My Passport/rbsp/pkls/') )[:1]:

    print '\n' + date

    for probe in ('a', 'b')[:1]:

      checkdate(probe, date, mpc=30)

  return





def col(x, width=12, unit=''):
  if isinstance(x, float):
    return format(x, str(width - len(unit) - 2) + '.1f') + unit + ' '
  else:
    return str(x).rjust(width - len(unit) - 1) + unit + ' '





# Check a given date for one probe. Break the day into chunks of mpc minutes. 
def checkdate(probe, date, mpc=30):

  today = day(probe=probe, date=date)

  print '\n' + col('probe') + col('time') + col('frequency') + col('magnitude') + col('coherence')

  # Scroll through each 20 minute chunk. 
  for t in range(0, 86400, 60*mpc):
#  for t in ( timeint(time='20:40:00'), timeint(time='20:50:00') ):

    print col(probe) + col( timestr(t)[1] ),

    ev = today.getslice(t, duration=60*mpc)

    # Make sure the data is OK. 
    if not ev.isok():
      print col('X')
      continue

    # Find the strongest component of the poloidal electric field. 
    i = np.argmax( ev.fft('ey') )
    print col(ev.frq()[i], unit='mHz'),

    # If it's not in the Pc4 range, bail. 
    if not ev.ispc4()[i]:
      print col('X')
      continue

    # How strong is the electric field at this point? Note that Lei based his
    # cutoff instead on the magnetic field. That's awkward for us since the
    # magnetic field has a node at the magnetic equator, which we're near. 
    print col( np.abs( ev.fft('ey')[i] ), unit='mV/m'),

    # If it's too small, bail. 
    if np.abs( ev.fft('ey')[i] ) < 0.25:
      print col('X')
      continue

    # Is this spectral component coherent? That's important if we want to judge
    # the phase accurately. 
    print col( ev.coh('p')[i] ),

    if ev.coh('p')[i] < 0.75:
      print col('X')
      continue

    # How does the phase offset of this Fourier mode line up with the mlat? 
    # This gives the parity of the harmonic mode. 
    print col( ev.harm('p')[i] ),

    if not ev.isodd('p')[i]:
      print col('X')
      continue




    print col('ok')


# #############################################################################
# ############################################################### Plot an Event
# #############################################################################

def plot(ev):
  # Create plot window to hold waveforms and spectra. 
  PW = plotWindow(nrows=3, ncols=2)
  # Put title and labels on the plot. 
  title = notex(ev.name)
  rowlabels = ( notex('Poloidal'), notex('Toroidal'), notex('Parallel') )
  collabels = ( notex('B (Red) ; E (Blue)'), 
                notex('Coh. (Black) ; CSD (Green) ; Par. (Violet)') )
  PW.setParams(title=title, collabels=collabels, rowlabels=rowlabels)
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
  # Show the plot. 
  return PW.render()

# #############################################################################
# ########################################################### For Importability
# #############################################################################

if __name__=='__main__':
  main()


