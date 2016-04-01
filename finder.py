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

  # Append Dst information to the list of events. 
  return adddst()

  # What dates do we have data for? 
  dates = sorted( os.listdir('/media/My Passport/rbsp/pkls/') )

#  # If we're saving our data, nuke the previous list to avoid double-counting. 
#  if '-i' in argv and os.path.exists('pos.txt'):
#    print 'Removing old position listing'
#    os.remove('pos.txt')

#  # Tally the probes' positions. 
#  for date in dates:
#    print date
#    [ trackpos(probe, date, mpc=5) for probe in ('a', 'b') ]

  # If we're saving our data, nuke the previous list to avoid double-counting. 
  if '-i' in argv and os.path.exists('events_sensitive.txt'):
    print 'Removing old event listing'
    os.remove('events_sensitive.txt')

  # Search for events. Do the days in random order, for easier debugging. We
  # can just look at the first event we find. 
  for date in np.random.permutation(dates):
    print date
    # Check both probes. Thirty minute chunks. 
    [ checkdate(probe, date, mpc=30) for probe in ('a', 'b') ]
  return


# #############################################################################
# ################################################## Add Dst Data to Event List
# #############################################################################

def adddst():
  dst = read('dst.txt')
  # Load Dst values into an array. The first column is epoch time. 
  dstarr = np.zeros( (2*len(dst) - 1, 2), dtype=np.int )
  for i, d in enumerate(dst):
    date, time, val = d.split()
    dstarr[2*i, 0] = timeint(date=date, time=time)
    dstarr[2*i, 1] = val
  # Average to get half hours. 
  dstarr[1::2] = 0.5*( dstarr[:-2:2] + dstarr[2::2] )
  # Read in the events. 
  events = read('events.txt')
  for e in events:
    date, time = e.split()[1:3]
    t = timeint(date=date, time=time)
    # Find the line of Dst that matches this. 
    i = np.argmin( np.abs( t - dstarr[:, 0] ) )
    # For a consistent number of columns, non-double events are labeled ONLY. 
    only = col('ONLY') if 'BIG' not in e and 'SMALL' not in e else ''
    # Print out a new file which includes Dst along with the event info. 
    append(e + only + col( dstarr[i, 1] ), 'events_with_dst.txt')
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
    evdicts = [ ev.standing(m, pc4=True, thresh=0.001) for m in ('p', 't') ]
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
  mag = col( d['s'], digs=5 )
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
    print append(text, 'events_sensitive.txt')
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
  mag = tex('imag') + tex('L3S') + '\\!=\\!' + fmt(d['s'], digs=3) + tex('mW/m^2')
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
  stand = [ ev.standing(m, pc4=True, thresh=0.001) for m in ('p', 't') ]
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


