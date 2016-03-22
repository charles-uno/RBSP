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
from scipy.optimize import curve_fit

from warnings import catch_warnings, simplefilter

# #############################################################################
# ######################################################################## Main
# #############################################################################

# A timestamped directory, in case we want to save any plots.  
plotdir = '/home/user1/mceachern/Desktop/plots/' + now() + '/'

def main():

#  # Plot the location of the usable data. 
#  return posplot(save='-i' in argv)

#  # Plot the location of events... ?
#  return rateplot(save='-i' in argv)

  # What dates do we have data for? 
  dates = sorted( os.listdir('/media/My Passport/rbsp/pkls/') )

  # If we're saving our data, nuke the previous list to avoid double-counting. 
  if '-i' in argv and os.path.exists('events.txt'):
    print 'Removing old event listing'
    os.remove('events.txt')

  # Search for events. Do the days in random order, for easier debugging. We
  # can just look at the first event we find. 
  for date in np.random.permutation(dates):

    date = '2012-10-23'

    print date

    # Check both probes. Thirty minute chunks. 
    [ checkdate(probe, date, mpc=30) for probe in ('a', 'b') ]

#  # Tally the probes' positions. 
#  for date in dates:
#    print date
#    # Figure where the probes spent their time. Five minute chunks. 
#    [ trackpos(probe, date, mpc=5) for probe in ('a', 'b') ]

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
# ====================================================== Check a Day for Events
# =============================================================================

# The day is broken into chunks (mpc is minutes per chunk). Each chunk is run
# through a few different filters to seek out odd-harmonic poloidal Pc4 waves. 
def checkdate(probe, date, mpc=30):

  # Load the day's data into a day object. 
  today = day(probe=probe, date=date)

  # Iterate over each chunk of the day. 
  for t in range(0, 86400, 60*mpc):

    print probe, timestr(t)[1]

    # Check for poloidal and toroidal events independently. 
    ev = today.getslice(t, duration=60*mpc)
    evlines = [ checkevent(ev, mode) for mode in ('p', 't') ]

    # If there's anything to save, do so, then plot it. 
    keepevent(evlines)

    '''
    # If this chunk of time has no events, bail. 
    if len(evlines)==0:
      continue
    # Otherwise, 
    # If there's an event, print it out and plot it. 
    if len(evlines)==1:
      # If we're saving the events... 
      if '-i' in argv:
        print append( ''.join( evdescr[0] ), 'events.txt')
        plot(ev, save=True)
      # If we're just debugging... 
      else:
        print ''.join( evdescr[0] )
        plot(ev)
        exit()
#    evline = checkevent(ev)
#    if not evline:
#      continue
    # Show or save the event info. 
    if '-i' in argv:
      print append(evline, 'events.txt')
      plot(ev, save=True)
    else:
      print evline
      plot(ev)
      exit()
    '''
  return




# =============================================================================
# ============================================================== Store an Event
# =============================================================================


def keepevent(evlines):
  # Filter out non-events.
  events = [ line for line in evlines if line ]
  # If there are no events, bail. 
  if len(events)==0:
    return 0
  elif len(events)==1:

    


    if len(evlines)==1:
      # If we're saving the events... 
      if '-i' in argv:
        print append( ''.join( evdescr[0] ), 'events.txt')
        plot(ev, save=True)
      # If we're just debugging... 
      else:
        print ''.join( evdescr[0] )
        plot(ev)
        exit()





  for line in evlines:
    print evline

  exit()

  return





# =============================================================================
# ================================================== Fit Waveform to a Gaussian
# =============================================================================

# Define Gaussian function. 
def gauss(x, amp, avg, std):
  return amp*np.exp( -(x - avg)**2/(2.*std**2) )

# Fit a Gaussian to data. The guess is amplitude, mean, spread. 
def gaussfit(x, y, guess=None):
  # Try to fit a peak around the highest value. 
  if guess is None:
    guess = (np.max(y), x[ np.argmax(y) ], 1)
  # Suppress warnings, but don't return bad data. 
  try:
    with catch_warnings():
      simplefilter('ignore')
      return curve_fit(gauss, x, y, p0=guess)[0]
  except:
    return None, None, None

# =============================================================================
# ====================================================== Check a Chunk of a Day
# =============================================================================

# Check if the given event -- which is a slice of a day's worth of data -- has
# any odd-mode poloidal Pc4 events. 
def checkevent(ev, mode):
  global plotdir

  # If the event's data is no good, obviously there can be no event. 
  if not ev.isok():
    print '\t' + mode + ': bad data. '
    return False

  # This filter is a bit trickier than Lei's, because we're going a level
  # deeper. We need not only the magnetic field, but also the electric field
  # and how they relate to one another. 

  # Get the Fourier-domain Poynting flux. Really, get the absolute value of the
  # imaginary component... we want to know the strength of the standing wave. 
  freq = ev.frq()
  sfft = np.abs( np.imag( ev.sfft(mode) ) )

  # Fit a Gaussian to the spectrum. Make sure it works. 
  speak, fpeak, dfpeak = gaussfit(freq, sfft)
  if None in (speak, fpeak, dfpeak):
    print '\t' + mode + ': fit failed. '
    return False

  # Threshold on frequency: 7 to 25 mHz.  
  if not 7 < fpeak < 25:
    print '\t' + mode + ': not Pc4. '
    return False

  # Threshold on magnitude: at least 0.01 mW/m^2 in the standing wave. 
  if not speak > 1e-2:
    print '\t' + mode + ': too weak. '
    return False

  # Threshold on fit quality: if the peaks are off by 5mHz or more, bail. 
  imax = np.argmax(sfft)
  if not np.abs(freq[imax] - fpeak) < 5:
    print '\t' + mode + ': bad fit. '
    return False

  # Threshold on harmonic: anything too close to the equator is ambiguous. 
  harm = ev.harm(mode)
  ipeak = np.argmin( np.abs(freq - fpeak) )
  if harm[ipeak]==0:
    print '\t' + mode + ': ambiguous harmonic. '
    return False

  # Threshold on coherence: anything below 0.9 is too messy. 
  if not ev.coh(mode)[ipeak] > 0.9:
    print '\t' + mode + ': incoherent. '
    return False

  # Assemble a list of strings about this event. It'll get joined later. 
  pdt = col(ev.probe) + col(ev.date) + col(ev.time)
  pos = col( ev.avg('lshell') ) + col( ev.avg('mlt') ) + col( ev.avg('mlat') )
  pol = col( {'p':'POL', 't':'TOR'}[mode] )
  par = col( {1:'ODD', 2:'EVEN'}[ harm[ipeak] ] )
  fdf = col(fpeak) + col(2.355*dfpeak)
  return pdt + pos + col(speak) + par + pol + fdf + col(ev.lpp)


  '''
  modes = ('p', 't')
  spti = [ np.abs( np.imag( ev.sfft(m) ) ) for m in modes ]
  # Fit a Gaussian to each spectrum. Make sure it works. 
  args = [ gaussfit(freq, s) for s in spti ]
  if None in args:
    return False
  # Which Gaussian is larger? Consider only that one from here on. 
  im = 0 if args[0][0] > args[1][0] else 1
  speak, fpeak, dfpeak = args[im]
  mode = modes[im]
  # Threshold on frequency: 7 to 25 mHz.  
  if not 7 < fpeak < 25:
    return False
  # Threshold on magnitude: at least 0.01 mW/m^2 in the standing wave. 
  if not speak > 1e-2:
    return False
  # Threshold on fit quality: if the peaks are off by 5mHz or more, bail. 
  imax = np.argmax( spti[im] )
  if not np.abs(freq[imax] - fpeak) < 5:
    return False
  # Threshold on harmonic: anything too close to the equator is ambiguous. 
  harm = ev.harm(mode)
  ipeak = np.argmin( np.abs(freq - fpeak) )
  if harm[ipeak]==0:
    return False
  # Threshold on coherence: anything below 0.9 is too messy. 
  if not ev.coh(mode)[ipeak] > 0.9:
    return False
  '''

# #############################################################################
# ############################################################### Plot an Event
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
  # Grab the Fourier-domain Poynting flux, scaled by L^3. 
  modes = ('p', 't')
  sift = [ np.abs( np.imag( ev.sfft(m) ) ) for m in modes ]
  srft = [ np.abs( np.real( ev.sfft(m) ) ) for m in modes ]
  # Do a Gaussian fit of the imaginary component. Identify the larger mode. 
  freq = ev.frq()
  args = [ gaussfit(freq, s) for s in sift ]
  im = 0 if args[0][0] > args[1][0] else 1
  # Scale the spectra to the largest value, in the data or in a fit. 
  smax = max( np.max(sift), np.max(srft), args[im][0] )
  # Plot the real and imaginary spectral components. 
  PW[:, 1].setParams( **ev.coords('spectra', cramped=True) )
  [ PW[i, 1].setLine(s/smax, 'm') for i, s in enumerate(sift) ]
  [ PW[i, 1].setLine(s/smax, 'g') for i, s in enumerate(srft) ]
  # Plot the Gaussian fit of the imaginary component. 
  f = np.linspace(0, 50, 1000)
  [ PW[i, 1].setLine(f, gauss( f, *args[i] )/smax, 'k--') for i in range(2) ]
  # Plot labels. 
  rowlabels = ( notex('Poloidal'), notex('Toroidal'), notex('Parallel') )
  collabels = ( notex('B (Red) ; E (Blue)'), 
                tex('imag') + tex('L3S') + notex(' (Magenta) ; ') + 
                tex('real') + tex('L3S') + notex(' (Green)') )
  PW.setParams(collabels=collabels, sidelabel=ev.label(), rowlabels=rowlabels)
  # Dig up some information to put together the title. 
  ig = np.argmin( np.abs(args[im][1] - freq) )
  modename = notex( {0:'POL', 1:'TOR'}[im] )
  freqname = 'f\\!=\\!' + format(args[im][1], '.1f') + tex('mHz')
  harmname = notex( {1:'ODD', 2:'EVEN'}[ ev.harm( modes[im] )[ig] ] )
  fwhmname = '\\delta\\!f\\!=\\!' + format(args[im][2], '.1f') + tex('mHz')
  sfft = srft[im] + 1j*sift[im]
  sizename = tex('L3S') + '\\!=\\!(' + cmt(sfft[ig], digs=2) + ')' + tex('mW/m^2')
  title = '\\qquad{}'.join( (modename, harmname, freqname, fwhmname, sizename) )
  PW.setParams(title=title)

  # Show the plot, or save it as an image. 
  if save is True:
    if not os.path.exists(plotdir):
      os.makedirs(plotdir)
    return PW.render(plotdir + ev.name + '.png')
  else:
    return PW.render()

# #############################################################################
# ######################################################## Plot RBSP's Position
# #############################################################################

def posplot(save=False, unit='days'):
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
  # Parse the positions into an array of floats. 
  pos = np.array( [ [ float(x) for x in p.split()[3:6] ] for p in poslines ] )
  # Figure out appropriate bin sizes for the position data. 
  dl, dm = 0.5, 1
  lmin, lmax = np.floor( np.min( pos[:, 0] ) ), np.ceil( np.max( pos[:, 0] ) )
  mmin, mmax = np.floor( np.min( pos[:, 1] ) ), np.ceil( np.max( pos[:, 1] ) )
#  lmin, lmax = 4.5, 6.5
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
    if not os.path.exists(plotdir):
      os.makedirs(plotdir)
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


