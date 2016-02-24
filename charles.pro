; Written by Aaron Breneman
; Changes by Charles McEachern
; Spring 2016

; Here, we rework the crib sheet from Aaron into a form that's easier for the Python script to handle. 




pro charles, date, probe

	rbspx = 'rbsp' + probe
	timespan,date

        if ~keyword_set(ql) then ql = 0
        if ~keyword_set(bp) then bp = '12'
        
;Stuff to make plots pretty	
	rbsp_efw_init
	!p.charsize = 1.2
	tplot_options,'xmargin',[20.,15.]
	tplot_options,'ymargin',[3,6]
	tplot_options,'xticklen',0.08
	tplot_options,'yticklen',0.02
	tplot_options,'xthick',2
	tplot_options,'ythick',2


;load eclipse times
	rbsp_load_eclipse_predict,probe,date
	get_data, rbspx + '_umbra',data=eu
	get_data, rbspx + '_penumbra',data=ep


;; ;Load spice
;; 	if ~keyword_set(no_spice_load) then rbsp_load_spice_kernels
;; 	rbsp_load_spice_state,probe=probe,coord='gse',/no_spice_load
;; 	store_data,'rbsp'+probe+'_state_pos_gse',newname='rbsp'+probe+'_pos_gse'
;; 	store_data,'rbsp'+probe+'_state_vel_gse',newname='rbsp'+probe+'_vel_gse'


;Load spinfit MGSE Efield and Bfield
        if ~keyword_set(nospinfit) then begin
           rbsp_efw_spinfit_vxb_subtract_crib,probe,/noplot,ql=ql,boom_pair=bp
           evar = rbspx + '_efw_esvy_mgse_vxb_removed_spinfit'
        endif else begin
           rbsp_efw_vxb_subtract_crib,probe,/noplot,ql=ql
           evar = rbspx + '_efw_esvy_mgse_vxb_removed'
	endelse
        
	
	;tplot,'rbsp'+probe+'_' + ['vxb_mgse', 'sfit12_mgse']
;	if ~keyword_set(noplot) then tplot,'rbsp'+probe+'_' + ['vxb_mgse','rbsp'+probe+'_efw_esvy_mgse_vxb_removed_spinfit']
	if ~keyword_set(noplot) then tplot, rbspx + '_' + ['vxb_mgse',evar]



;Grab the spinfit Ew and Bw data
	split_vec, rbspx + '_mag_mgse'
;	get_data,'rbsp'+probe+'_efw_esvy_mgse_vxb_removed_spinfit',data=sfit
	get_data,evar,data=edata


        
        if is_struct(edata) then tinterpol_mxn, rbspx + '_mag_mgse',edata.x,newname=rbspx + '_mag_mgse'
        

                                ;smooth the background magnetic field
                                ;over 30 min for the E*B=0 calculation
        rbsp_detrend,rbspx + '_mag_mgse',60.*30.


	get_data,rbspx + '_mag_mgse',data=magmgse
        get_data,rbspx + '_mag_mgse_smoothed',data=magmgse_smoothed

	if ~is_struct(magmgse) then begin
		print,'NO MAG DATA FOR rbsp_efw_EdotB_to_zero_crib.pro TO USE...RETURNING'
		return		
	endif
		
	
	bmag = sqrt(magmgse.y[*,0]^2 + magmgse.y[*,1]^2 + magmgse.y[*,2]^2)
	bmag_smoothed = sqrt(magmgse_smoothed.y[*,0]^2 + magmgse_smoothed.y[*,1]^2 + magmgse_smoothed.y[*,2]^2)


;Replace axial measurement with E*B=0 version
	edata.y[*,0] = -1*(edata.y[*,1]*magmgse_smoothed.y[*,1] + edata.y[*,2]*magmgse_smoothed.y[*,2])/magmgse_smoothed.y[*,0]
	;; if ~keyword_set(suffix) then store_data,'rbsp'+probe+'_efw_esvy_mgse_vxb_removed_spinfit',data=edata
        ;; if keyword_set(suffix) then store_data,'rbsp'+probe+'_efw_esvy_mgse_vxb_removed_spinfit'+'_'+suffix,data=edata
	if ~keyword_set(suffix) then store_data,evar,data=edata
        if keyword_set(suffix) then store_data,evar+'_'+suffix,data=edata


;Find bad E*B=0 data (where the angle b/t spinplane MGSE and Bo is less than 15 deg) 
;Good data has By/Bx < 3.732   and  Bz/Bx < 3.732

	By2Bx = abs(magmgse_smoothed.y[*,1]/magmgse_smoothed.y[*,0])
	Bz2Bx = abs(magmgse_smoothed.y[*,2]/magmgse_smoothed.y[*,0])
	store_data,'B2Bx_ratio',data={x:edata.x,y:[[By2Bx],[Bz2Bx]]}
	ylim,'B2Bx_ratio',0,10
	options,'B2Bx_ratio','ytitle','By/Bx (black)!CBz/Bx (red)'
	badyx = where(By2Bx gt 3.732)
	badzx = where(Bz2Bx gt 3.732)


;calculate angles b/t despun spinplane antennas and Bo. 
	n = n_elements(edata.x)
	ang_ey = fltarr(n)
	ang_ez = fltarr(n)

	for i=0L,n-1 do ang_ey[i] = acos(total([0,1,0]*magmgse_smoothed.y[i,*])/(bmag_smoothed[i]))/!dtor
	for i=0L,n-1 do ang_ez[i] = acos(total([0,0,1]*magmgse_smoothed.y[i,*])/(bmag_smoothed[i]))/!dtor
	store_data,'angles',data={x:edata.x,y:[[ang_ey],[ang_ez]]}



;Calculate ratio b/t spinaxis and spinplane components
	e_sp = sqrt(edata.y[*,1]^2 + edata.y[*,2]^2)
	rat = abs(edata.y[*,0])/e_sp
	store_data,'rat',data={x:edata.x,y:rat}
	store_data,'e_sp',data={x:edata.x,y:e_sp}
	store_data,'e_sa',data={x:edata.x,y:abs(edata.y[*,0])}


;Check for Spinfit saturation
	;; get_data,'rbsp'+probe+'_efw_esvy_mgse_vxb_removed_spinfit',data=tmpp
	get_data,evar,data=tmpp
	badsatx = where(abs(tmpp.y[*,0]) ge 195.)
	badsaty = where(abs(tmpp.y[*,1]) ge 195.)
	badsatz = where(abs(tmpp.y[*,2]) ge 195.)



;Remove bad Efield data
;....saturated data from the rest of the tplot variables
;....saturated data from Ex
;....Ex data when the E*B=0 calculation is unreliable

	;; if ~keyword_set(suffix) then get_data,'rbsp'+probe+'_efw_esvy_mgse_vxb_removed_spinfit',data=tmpp
	;; if keyword_set(suffix) then  get_data,'rbsp'+probe+'_efw_esvy_mgse_vxb_removed_spinfit_'+suffix,data=tmpp
	if ~keyword_set(suffix) then get_data,evar,data=tmpp
	if keyword_set(suffix) then  get_data,evar+'_'+suffix,data=tmpp
	if badyx[0] ne -1 then tmpp.y[badyx,0] = !values.f_nan
	if badzx[0] ne -1 then tmpp.y[badzx,0] = !values.f_nan
	if badsatx[0] ne -1 then tmpp.y[badsatx,0] = !values.f_nan
	if badsaty[0] ne -1 then tmpp.y[badsaty,1] = !values.f_nan
	if badsatz[0] ne -1 then tmpp.y[badsatz,2] = !values.f_nan
	;; if ~keyword_set(suffix) then store_data,'rbsp'+probe+'_efw_esvy_mgse_vxb_removed_spinfit',data=tmpp
	;; if keyword_set(suffix) then  store_data,'rbsp'+probe+'_efw_esvy_mgse_vxb_removed_spinfit_'+suffix,data=tmpp
	if ~keyword_set(suffix) then store_data,evar,data=tmpp
	if keyword_set(suffix) then  store_data,evar+'_'+suffix,data=tmpp

	get_data,'rat',data=tmpp
	if badyx[0] ne -1 then tmpp.y[badyx] = !values.f_nan
	if badzx[0] ne -1 then tmpp.y[badzx] = !values.f_nan
	if badsatx[0] ne -1 then tmpp.y[badsatx] = !values.f_nan
	if badsaty[0] ne -1 then tmpp.y[badsaty] = !values.f_nan
	if badsatz[0] ne -1 then tmpp.y[badsatz] = !values.f_nan
	store_data,'rat',data=tmpp

	get_data,'e_sa',data=tmpp
	if badyx[0] ne -1 then tmpp.y[badyx] = !values.f_nan
	if badzx[0] ne -1 then tmpp.y[badzx] = !values.f_nan
	if badsatx[0] ne -1 then tmpp.y[badsatx] = !values.f_nan
	if badsaty[0] ne -1 then tmpp.y[badsaty] = !values.f_nan
	if badsatz[0] ne -1 then tmpp.y[badsatz] = !values.f_nan
	store_data,'e_sa',data=tmpp

	get_data,'e_sp',data=tmpp
	if badyx[0] ne -1 then tmpp.y[badyx] = !values.f_nan
	if badzx[0] ne -1 then tmpp.y[badzx] = !values.f_nan
	if badsatx[0] ne -1 then tmpp.y[badsatx] = !values.f_nan
	if badsaty[0] ne -1 then tmpp.y[badsaty] = !values.f_nan
	if badsatz[0] ne -1 then tmpp.y[badsatz] = !values.f_nan
	store_data,'e_sp',data=tmpp




;Remove corotation field
	;; dif_data,'rbsp'+probe+'_efw_esvy_mgse_vxb_removed_spinfit','rbsp'+probe+'_E_coro_mgse',newname='rbsp'+probe+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit'
        ;; if keyword_set(suffix) then dif_data,'rbsp'+probe+'_efw_esvy_mgse_vxb_removed_spinfit_'+suffix,'rbsp'+probe+'_E_coro_mgse',newname='rbsp'+probe+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_'+suffix

        if ~keyword_set(nospinfit) then begin
           dif_data,evar, rbspx + '_E_coro_mgse',newname=rbspx + '_efw_esvy_mgse_vxb_removed_coro_removed_spinfit'
           if keyword_set(suffix) then dif_data,evar+'_'+suffix, rbspx + '_E_coro_mgse',newname=rbspx + '_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_'+suffix
        endif else begin
           dif_data,evar, rbspx + '_E_coro_mgse',newname=rbspx + '_efw_esvy_mgse_vxb_removed_coro_removed'
           if keyword_set(suffix) then dif_data,evar+'_'+suffix, rbspx + '_E_coro_mgse',newname=rbspx + '_efw_esvy_mgse_vxb_removed_coro_removed_'+suffix
        endelse

;Plot results

	options,'rat','ytitle','|Espinaxis|/!C|Espinplane|'
	options,'e_sp','ytitle','|Espinplane|'
	options,'e_sa','ytitle','|Espinaxis|'
	options,'angles','ytitle','angle b/t Ey & Bo!CEz & Bo (red)'
;	ylim,'rbsp'+probe+'_efw_esvy_mgse_vxb_removed_spinfit',-10,10
	ylim,evar,-10,10
	ylim, rbspx + '_mag_mgse',-200,200
	ylim,['e_sa','e_sp','rat'],0,10
;	options,'rbsp'+probe+'_efw_esvy_mgse_vxb_removed_spinfit','labels',['xMGSE','yMGSE','zMGSE']
	options,evar,'labels',['xMGSE','yMGSE','zMGSE']


        tplot_options,'title','RBSP-'+probe + ' ' + date
        if ~keyword_set(noplot) then begin
           if ~keyword_set(nospinfit) then begin
              tplot,[rbspx + '_mag_mgse',$
                     rbspx + '_mag_mgse_smoothed',$
                     rbspx + '_efw_esvy_mgse_vxb_removed_spinfit',$
                     rbspx + '_efw_esvy_mgse_vxb_removed_coro_removed_spinfit',$
                     rbspx + '_E_coro_mgse',$
                     'angles',$
                     'rat',$
                     'e_sa',$
                     'e_sp']
           endif else begin
              tplot,[rbspx + '_mag_mgse',$
                     rbspx + '_mag_mgse_smoothed',$
                     rbspx + '_efw_esvy_mgse_vxb_removed',$
                     rbspx + '_efw_esvy_mgse_vxb_removed_coro_removed',$
                     rbspx + '_E_coro_mgse',$
                     'angles',$
                     'rat',$
                     'e_sa',$
                     'e_sp']

           endelse
           if keyword_set(eu) then timebar,eu.x
           if keyword_set(eu) then timebar,eu.x + eu.y
        endif



;--------------------------------------------------
;--------------------------------------------------
;--------------------------------------------------
;;;****Added stuff for Charles


get_data,rbspx + '_Lvec',data=wgse
wgse = wgse.y
rbsp_mgse2gse, rbspx + '_efw_esvy_mgse_vxb_removed_spinfit',wgse,newname='efield_spinfit_gse'
rbsp_mgse2gse, rbspx + '_mag_mgse',wgse,newname='bfield_gse'

;;interpolate to common set of times
tinterpol_mxn,'bfield_gse','efield_spinfit_gse',newname='bfield_gse'
tinterpol_mxn, rbspx + '_state_pos_gse','efield_spinfit_gse',newname=rbspx + '_state_pos_gse'


;;These are the three tplot variables you requested.
tplot,[rbspx + '_state_pos_gse','efield_spinfit_gse','bfield_gse']

;--------------------------------------------------
;--------------------------------------------------
;--------------------------------------------------

end



;--------------------------------------------------
;--------------------------------------------------
;--------------------------------------------------
;;;****Added stuff for Charles


get_data,'rbsp'+probe+'_Lvec',data=wgse
wgse = wgse.y
rbsp_mgse2gse,'rbsp'+probe+'_efw_esvy_mgse_vxb_removed_spinfit',wgse,newname='efield_spinfit_gse'
rbsp_mgse2gse,'rbsp'+probe+'_mag_mgse',wgse,newname='bfield_gse'

;;interpolate to common set of times
tinterpol_mxn,'bfield_gse','efield_spinfit_gse',newname='bfield_gse'
tinterpol_mxn,'rbsp'+probe+'_state_pos_gse','efield_spinfit_gse',newname='rbsp'+probe+'_state_pos_gse'


;;These are the three tplot variables you requested.
tplot,['rbsp'+probe+'_state_pos_gse','efield_spinfit_gse','bfield_gse']

;--------------------------------------------------
;--------------------------------------------------
;--------------------------------------------------

end

