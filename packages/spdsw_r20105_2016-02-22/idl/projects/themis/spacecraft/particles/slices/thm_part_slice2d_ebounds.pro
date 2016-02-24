;+
;Procedure:
;  thm_part_slice2d_ebounds
;
;
;Purpose:
;  Returns an array of gapless energy boundaries.  The number of 
;  elements returned will always be N+1 for N energy levels.
;
;
;Input:
;  dist: 3D particle data structure
;
;
;Output:
;  return value: Array of energy bin boundaries (# energy bins + 1)
;
;
;Notes:
;  Energy levels are ordered differently between ESA/SST
;
;
;$LastChangedBy: aaflores $
;$LastChangedDate: 2013-11-05 17:44:02 -0800 (Tue, 05 Nov 2013) $
;$LastChangedRevision: 13496 $
;$URL: svn+ssh://thmsvn@ambrosia.ssl.berkeley.edu/repos/spdsoft/trunk/projects/themis/spacecraft/particles/slices/thm_part_slice2d_ebounds.pro $
;
;-
function thm_part_slice2d_ebounds, dist

    compile_opt idl2, hidden

  n = dimen1(dist.energy)

  energies = fltarr(size(dist.energy,/dim)+[1,0])
  
  ; use midpoints
  energies[1:n-1,*] = (dist.energy[0:n-2,*] + dist.energy[1:n-1,*]) / 2.
  
  ; top/bottom energies
  energies[0,*] = dist.energy[0,*] + (dist.energy[0,*] - energies[1,*])
  energies[n,*] = dist.energy[n-1,*] + (dist.energy[n-1,*] - energies[n-1,*])
  
  return, energies
  
end
