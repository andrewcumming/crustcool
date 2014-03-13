pro pl

	readcol, 'data/fit_joint_wchandra_BB_2012july28_pu.out', time, T, T1,T2, uflux, uflux1, uflux2, A,Aerr, $
				pflux, pfluxerr, format=('X,D,X,X,X,X,X,X,X,D,D,D,D,D,D,X,X,X,D,D,X,X,D,D')

	time -= 55756.53318
	
	fnorm = 4.0*!dpi*(1.6*3.086d21)^2
	uflux *= fnorm
	uflux1 *= fnorm
	uflux2 *= fnorm

	pflux *=fnorm
	pfluxerr *= fnorm

	plot, time, uflux, psym=1, /xlog,/ylog, xtitle=textoidl('Time (d)'),$
			ytitle=textoidl('Unabsorbed 1-10 keV luminosity (erg s^{-1})'),$
				yrange = [1d33,1d36], ystyle=1
	oploterror, time, uflux, uflux2-uflux, /hibar, psym=1
	oploterror, time, uflux, uflux-uflux1, /lobar, psym=1

	oploterror, time, pflux, pfluxerr, psym=3

end



pro twist2, ps=ps
	
	!p.multi=[0,1,2,0,0]
	!p.charsize=1.4

	if keyword_set(ps) then begin
		open_ps, 'twist2.ps'
	endif
	
	readcol, 'data/fit_joint_wchandra_BB_2012july28_pu.out', time, T, T1,T2, uflux, uflux1, uflux2, A,Aerr, $
					format=('X,D,X,X,X,X,X,X,X,D,D,D,D,D,D,X,X,X,D,D')
	
	time -= 55756.53318
	
	Anorm=4.0*!dpi*(1.6*3.086d21/1d6)^2
	A*=Anorm
	Aerr*=Anorm

	fnorm = 4.0*!dpi*(1.6*3.086d21)^2
	uflux *= fnorm
	uflux1 *= fnorm
	uflux2 *= fnorm


	mtime = 1.0*dindgen(200)
	mA=0.12*(1.0-mtime/30.0)
	mA[where(mA<0.0)] = 0.0
	mL=1.3d35*50.0*mA^2

	plot, time, uflux, psym=1, /xlog,/ylog, xtitle=textoidl('Time (d)'),$
			ytitle=textoidl('Unabsorbed 1-10 keV luminosity (erg s^{-1})'),$
				yrange = [1d33,1d36], ystyle=1
	oploterror, time, uflux, uflux2-uflux, /hibar, psym=1
	oploterror, time, uflux, uflux-uflux1, /lobar, psym=1
	
	oplot, mtime, mL, linestyle=1

	ploterror, time, A, Aerr, psym=1,/xlog, $
			ytitle=textoidl('Area (10^{12} cm^2)'), $
			xtitle=textoidl('Time (d)'), yrange=[0.0,0.25]
	oplot, mtime, mA, linestyle=1

	if keyword_set(ps) then begin
		close_ps
	endif

end





pro twist,ps=ps
	
	!p.multi=[0,1,1,0,0]
	!p.charsize=1.4
	
	if keyword_set(ps) then begin
		;open_ps, 'twist.ps'
		set_plot,'ps'
		device,filename='twist.ps',/color
		loadct,39
		!p.thick=3
		!p.charthick=3
		!x.thick=3
		!y.thick=3
	endif
	
	readcol, 'data/j1647_fit_data.dat', time2, TT2, TTerr2, BBnorm2, BBnormerr2, BBrad2, BBraderr2,$
					format=('D,X,X,D,D,D,D,D,D')

	time2 -= time2[0] + 5.0
	ind = where(time2 gt 0.0)
	time2=time2[ind]
	BBnorm2=BBnorm2[ind]*1d39*(5.0/10.0)^2
	BBnormerr2=BBnormerr2[ind]*1d39*(5.0/10.0)^2
	BBrad2=sqrt(BBrad2[ind])*(5.0/10.0)  ; adjust to 5kpc
	K2=(0.1*BBrad2)^(-0.5)
	Te2=TT2[ind]/K2
	
	A2 = 4.0*!dpi*BBrad2^2/100.0
	
	readcol, 'data/fit_joint_wchandra_BB_2012july28_pu.out', time, T, T1,T2, uflux, uflux1, uflux2, A,Aerr, $
					format=('X,D,X,X,X,X,X,X,X,D,D,D,D,D,D,X,X,X,D,D')
	
	time -= 55756.53318
	
	Anorm=(1.6*3.086d21/1d6)^2
	A*=Anorm
	Aerr*=Anorm
	A*=4.0*!dpi
	Aerr*=4.0*!dpi

	fnorm = 4.0*!dpi*(1.6*3.086d21)^2
	uflux *= fnorm
	uflux1 *= fnorm
	uflux2 *= fnorm

	plot, A, uflux, psym=1, /xlog, /ylog, xtitle=textoidl('Area (10^{12} cm^2)'),$
			ytitle=textoidl('Unabsorbed 1-10 keV luminosity (erg s^{-1})'),$
				yrange = [1d33,1d35], ystyle=1, xrange=[1d-3,1.0]
	oploterror, A, uflux, Aerr, uflux2-uflux, /hibar, psym=1
	oploterror, A, uflux, Aerr, uflux-uflux1, /lobar, psym=1
	
	A=-3.0+2.0*dindgen(100)
	F=alog10(1.3)+35.0+2.0*A   
	oplot, 10^A, 50*10^F, linestyle=1

	oplot, A2, BBnorm2, psym=5,col=250

	if keyword_set(ps) then begin
		close_ps
	endif	

end



pro fc,ps=ps
	
	!p.multi=[0,2,2,0,0]
	!p.charsize=1.4
	
	
	if keyword_set(ps) then begin
		;open_ps, 'twist.ps'
		set_plot,'ps'
		device,filename='fc.ps'
		!p.thick=3
		!p.charthick=3
		!x.thick=3
		!y.thick=3
		!p.charsize=1.01
	endif
	
	readcol, 'data/NSmax.dat', time3, TT3, TTerr3, AA3, format=('D,X,D,D,X,X,D')
	time3 -= time3[0] + 5.0
	TT3 = 10^TT3
	TT3 /= 1.6e-9/1.38d-16
	print, AA3
	
	; Data for 1647
	; from Hongjun:
	; Time   Nh    err   kT   err   BBnorm   err   BB_rad   err   Gamma   err   PLnorm   err, 
	; where err means the fit uncertainty of the quantity in the previous column. I used XSPEC models 
	; to fit the data and thus the BBnorm has a specific definition (Luminosity in units of 10^39 erg/s 
	; at a distance of 10 kpc).	Note that the data start from 5 days before the 2006 outburst 
	; (the first row corresponds to the state of the star 5 days before the outburst).
	readcol, 'data/j1647_fit_data.dat', time2, TT2, TTerr2, BBnorm2, BBnormerr2, BBrad2, BBraderr2,$
					format=('D,X,X,D,D,D,D,D,D')
	time2 -= time2[0] + 5.0
	ind = where(time2 gt 0.0)
	time2=time2[ind]
	BBnorm2=BBnorm2[ind]*1d39*(5.0/10.0)^2
	BBnormerr2=BBnormerr2[ind]*1d39*(5.0/10.0)^2
	BBrad2=sqrt(BBrad2[ind])*(5.0/10.0)  ; adjust to 5kpc
	K2=(0.1*BBrad2)^(-0.5)
	Te2=TT2[ind]/K2
	
	readcol, 'data/fit_joint_wchandra_BB_2012july28_pu.out', time, T, T1,T2, uflux, uflux1, uflux2, A,Aerr, $
					format=('X,D,X,X,X,X,X,X,X,D,D,D,D,D,D,X,X,X,D,D')
	
	time -= 55756.53318
	
	Anorm=(1.6*3.086d21/1d6)^2
	A*=Anorm
	Aerr*=Anorm

	fnorm = 4.0*!dpi*(1.6*3.086d21)^2
	uflux *= fnorm
	uflux1 *= fnorm
	uflux2 *= fnorm
			
	K = A^(-0.25)
	Kerr = K*0.25 * Aerr/A	
	Terr = 0.5*(T2-T1)
	Te=T
	Teerr = Te
	Te=T/K
	Teerr = Te * sqrt((Terr/T)^2 + (Kerr/K)^2)
			
	plot, Te, K, psym=1, ytitle=textoidl('K^{-1/4}'), xtitle=textoidl('T (keV)'), yrange=[0,7],$
			xrange=[0,0.3], xstyle=1, ystyle=1
	oploterror, Te, K, Terr,Kerr, psym=1

	oplot, Te2, K2, psym=5, col=250
	;oplot, TT3*sqrt(AA3), 1.0/sqrt(AA3), psym=5, col=120

	;ploterror, 100.0-time, A, Aerr, psym=1, $
	;		ytitle=textoidl('(R/d)^2 (10 km/1.6 kpc)^2'), $
	;		xtitle=textoidl('Time (d)'), xrange=[0,100], xstyle=1
	
	R = sqrt(A)*10.0
	Rerr = R*0.5*(Aerr/A)
	ploterror, time, R, Rerr, psym=1,/xlog, $
			ytitle=textoidl('R (km)'), $
			xtitle=textoidl('Time (d)'), yrange=[0,3], xrange=[1,3000],xstyle=1
		oplot, time2, BBrad2, psym=5, col=250
		;oplot, time3, sqrt(AA3), psym=5, col=120

	ploterror, time, T, Terr, psym=1, /xlog, $
				ytitle=textoidl('T (keV)'), $
				xtitle=textoidl('Time (d)'), xrange=[1,3000], yrange=[0.4,0.9],xstyle=1
			oploterror, time2, TT2, TTerr2, psym=5, col=250, errcol=250
			;oplot, time3, TT3, psym=5, col=120

	plot, time, uflux, psym=1, /xlog,/ylog, xtitle=textoidl('Time (d)'),$
			ytitle=textoidl('L (1-10 keV) (erg s^{-1})'),$
			yrange = [1d33,1d35], ystyle=1, xrange=[1,3000],xstyle=1
	oploterror, time, uflux, uflux2-uflux, /hibar, psym=1
	oploterror, time, uflux, uflux-uflux1, /lobar, psym=1
	oplot, time2, BBnorm2, col=250, psym=1
			
			
;	A*=4.0*!dpi
;	Aerr*=4.0*!dpi
;	plot, A, uflux, psym=1, /xlog, /ylog, xtitle=textoidl('Area (10^{12} cm^2)'),$
;			ytitle=textoidl('Unabsorbed 1-10 keV luminosity (erg s^{-1})'),$
;;				yrange = [1d33,1d35], ystyle=1
;	oploterror, A, uflux, Aerr, uflux2-uflux, /hibar, psym=1
;	oploterror, A, uflux, Aerr, uflux-uflux1, /lobar, psym=1
;;	
;	A=-3.0+2.0*dindgen(100)
;	F=alog10(1.3)+35.0+2.0*A   
;;	oplot, 10^A, 50*10^F, linestyle=1
	
if keyword_set(ps) then begin
	close_ps
endif	


end



pro fce,ps=ps
	
;	!p.multi=[0,2,2,0,0]
	!p.multi=[0,1,1,0,0]
	!p.charsize=1.4
	
	if keyword_set(ps) then begin
		;open_ps, 'twist.ps'
		set_plot,'ps'
		device,filename='fce.ps', /color
		!p.thick=3
		!p.charthick=3
		!x.thick=3
		!y.thick=3
		!p.charsize=1.01
		
	endif
	
	readcol, 'data/j1647_fit_data.dat', time2, TT2, TTerr2, BBnorm2, BBnormerr2, BBrad2, BBraderr2,$
					format=('D,X,X,D,D,D,D,D,D')

	time2 -= time2[0] + 5.0
	ind = where(time2 gt 0.0)
	time2=time2[ind]
	BBnorm2=BBnorm2[ind]*1d39
	BBnormerr2=BBnormerr2[ind]*1d39
	BBrad2=sqrt(BBrad2[ind])*(5.0/10.0)  ; adjust to 5kpc
	K2=(0.1*BBrad2)^(-0.5)
	K2err = K2*0.25 * BBraderr2[ind]/BBrad2	
	Te2=TT2[ind]/K2
	Te2err = Te2 * sqrt((TTerr2/TT2)^2 + (K2err/K2)^2)
	
	readcol, 'data/fit_joint_wchandra_BB_2012july28_pu.out', time, T, T1,T2, uflux, uflux1, uflux2, A,Aerr, $
					format=('X,D,X,X,X,X,X,X,X,D,D,D,D,D,D,X,X,X,D,D')
	
	time -= 55756.53318
	
	Anorm=(1.6*3.086d21/1d6)^2
	A*=Anorm
	Aerr*=Anorm

	fnorm = 4.0*!dpi*(1.6*3.086d21)^2
	uflux *= fnorm
	uflux1 *= fnorm
	uflux2 *= fnorm
			
	K = A^(-0.25)
	Kerr = K*0.25 * Aerr/A	
	Terr = 0.5*(T2-T1)
	Te=T/K
	Teerr = Te * sqrt((Terr/T)^2 + (Kerr/K)^2)
			
	plot, Te, K/min(K), psym=1, ytitle=textoidl('K^{-1/4}'), $
				xtitle=textoidl('T (keV)/ K^{-1/4}'), $
				xrange=[0,0.3], xstyle=1, yrange=[0.8,3.5], ystyle=1
	oploterror, Te, K/min(K), Teerr,Kerr/min(K), psym=1
	oploterror, Te2, K2/min(K2), Te2err, K2err/min(K2), psym=5, col=250, errcol=250


	;ploterror, 100.0-time, A, Aerr, psym=1, $
	;		ytitle=textoidl('(R/d)^2 (10 km/1.6 kpc)^2'), $
	;		xtitle=textoidl('Time (d)'), xrange=[0,100], xstyle=1
	
	R = sqrt(A)*10.0
	Rerr = R*0.5*(Aerr/A)
;	ploterror, time, R, Rerr, psym=1,/xlog, $
;			ytitle=textoidl('R/d (km/1.6 kpc)'), $
;;			xtitle=textoidl('Time (d)'), yrange=[0,2.5]
;		oplot, time2, BBrad2, psym=5


;	ploterror, time, Te, Teerr, psym=1, /xlog, $
;				ytitle=textoidl('T (keV)/K^{-1/4}'), $
;				xtitle=textoidl('Time (d)')
;
;			oplot, time2, T2, psym=5

;			plot, time, uflux, psym=1, /xlog,/ylog, xtitle=textoidl('Time (d)'),$
;					ytitle=textoidl('L (1-10 keV) (erg s^{-1})'),$
;						yrange = [1d33,1d35], ystyle=1
;			oploterror, time, uflux, uflux2-uflux, /hibar, psym=1
;			oploterror, time, uflux, uflux-uflux1, /lobar, psym=1


;	A*=4.0*!dpi
;	Aerr*=4.0*!dpi
;	plot, A, uflux, psym=1, /xlog, /ylog, xtitle=textoidl('Area (10^{12} cm^2)'),$
;			ytitle=textoidl('Unabsorbed 1-10 keV luminosity (erg s^{-1})'),$
;;				yrange = [1d33,1d35], ystyle=1
;	oploterror, A, uflux, Aerr, uflux2-uflux, /hibar, psym=1
;	oploterror, A, uflux, Aerr, uflux-uflux1, /lobar, psym=1
;;	
;	A=-3.0+2.0*dindgen(100)
;	F=alog10(1.3)+35.0+2.0*A   
;;	oplot, 10^A, 50*10^F, linestyle=1
	if keyword_set(ps) then begin
		close_ps
	endif	
	

end




pro mygrid

	erase
	readcol, 'out/grid', y, T, F,rho, format=('D,D,D,D')

	y=10^y
	T=10^T
;	rho=10^rho

	plot, rho,T,/xlog,/ylog, psym=1, xrange=[0.1,1d12], yrange=[1d5,3d9],ystyle=1,xstyle=1
	
end





pro tc,source=source,ps=ps,noplot=noplot,noextras=noextras
	; plots the lightcurve as T_eff vs time
	; source ='1659', '1731', 'XTEJ'

		if keyword_set(ps) then open_ps, 'tc.ps'

!p.multi=[0,1,1,0,0]

	if (not keyword_set(source)) then source='1659'

	if (source eq '1659' or source eq '1731') then begin
		xr=[1.0,10000.0]
		yr=[40,140]
		xr=[30.0,10000.0]
		yr=[40,110]
	endif
	if (source eq '0748') then begin
		xr=[10.0,100000.0]
		yr=[80,140]
	endif
	if (source eq '0556') then begin
		xr=[1.0,5000.0]
		yr=[150,370]
	endif


	if (source eq 'xtej') then begin
		xr=[1.0,6000.0]
		yr=[100,185]
	endif
	if (source eq 'terz3' or source eq 'terz2' or source eq 'terz') then begin
		xr=[150.0,5000.0]
		yr=[67,104]
	endif

	if (source eq 'XTEJ2') then begin
		xr=[0.1,100.0]
		yr=[100,180]
	endif
	

	; read in the observations and plot
	readcol, 'data/'+source, t0, n, format=('D,I'), numline=1
	readcol, 'data/'+source, t, F, Fe, temp, tempe, skipline=1, format=('D,D,D,D,D')	
	t=t-t0[0]
;	ploterror, t, temp, tempe, psym=1, /xlog, xtitle=textoidl('Time (d)'),$
;			ytitle=textoidl('T_{eff} (keV)'), charsize=1.4, xrange=xr,xstyle=1,$
;			yrange=yr, ystyle=1,/ylog

		ploterror, t, temp, tempe, psym=1, /xlog, xtitle=textoidl('Time (d)'),$
				ytitle=textoidl('T_{eff} (eV)'), charsize=1.4, xrange=xr,xstyle=1,$
				yrange=yr, ystyle=1,/nohat


		if not keyword_set(noplot) then begin
		; plot lightcurve from simulation
		readcol, 'gon_out/prof', tt,Teff, format=('F,X,X,X,F')
		tt/=(24*3600.0)
		; the output is already redshifted
		; but needs to be converted to eV
		FF = 1.38d-16*Teff/(1.6d-12)
		oplot, tt,FF,linestyle=0,col=250
		endif
	
		if (source eq 'terz3' or source eq 'terz2' or source eq 'terz') then begin
			readcol, 'gon_out/prof_terz', tt,Teff, format=('F,X,X,X,F')
			tt/=(24*3600.0)
			; the output is already redshifted
			; but needs to be converted to eV
			FF = 1.38d-16*Teff/(1.6d-12)
			oplot, tt,FF,linestyle=0


			readcol, 'gon_out/prof_terz_2', tt,Teff, format=('F,X,X,X,F')
			tt/=(24*3600.0)
			; the output is already redshifted
			; but needs to be converted to eV
			FF = 1.38d-16*Teff/(1.6d-12)
			oplot, tt,FF,linestyle=0

		endif



	;	readcol, 'gon_out/prof_1659_new', tt,Teff, format=('F,X,X,X,F')
	;	tt/=(24*3600.0)
	;	FF = 1.38d-16*Teff/(1.6d-12)
	;	oplot, tt,FF,linestyle=0

	if (source eq '1731' and 0) then begin

		readcol, 'gon_out/prof_1731_chandra1', tt,Teff, format=('F,X,X,X,F')
		tt/=(24*3600.0)
		FF = 1.38d-16*Teff/(1.6d-12)
		oplot, tt,FF,linestyle=0
		
		openw, lun,'lc_1731_chandra1.dat', /get_lun
		for i=0,n_elements(tt)-1 do begin
			if tt[i] gt 0.0 then printf, lun, tt[i], FF[i], format='(g,g)'
		endfor
		free_lun, lun
		
		readcol, 'gon_out/prof_1731_chandra2', tt,Teff, format=('F,X,X,X,F')
		tt/=(24*3600.0)
		FF = 1.38d-16*Teff/(1.6d-12)
		oplot, tt,FF,linestyle=1

		openw, lun,'lc_1731_chandra2.dat', /get_lun
		for i=0,n_elements(tt)-1 do begin
			if tt[i] gt 0.0 then printf, lun, tt[i], FF[i], format='(g,g)'
		endfor
		free_lun, lun

		readcol, 'gon_out/prof_1731_chandra3', tt,Teff, format=('F,X,X,X,F')
		tt/=(24*3600.0)
		FF = 1.38d-16*Teff/(1.6d-12)
		oplot, tt,FF,linestyle=4

		openw, lun,'lc_1731_chandra3.dat', /get_lun
		for i=0,n_elements(tt)-1 do begin
			if tt[i] gt 0.0 then printf, lun, tt[i], FF[i], format='(g,g)'
		endfor
		free_lun, lun


		oplot, [5230.0,5230.0],[10.0,1000.0],linestyle=2

	endif



	if (source eq '1731' ) then begin

		readcol, 'gon_out/prof_mcmc1', tt,Teff, format=('F,X,X,X,F')
		tt/=(24*3600.0)
		FF = 1.38d-16*Teff/(1.6d-12)
		oplot, tt,FF,linestyle=0,col=80
		
		openw, lun,'lc_1731_mcmc1.dat', /get_lun
		for i=0,n_elements(tt)-1 do begin
			if tt[i] gt 0.0 then printf, lun, tt[i], FF[i], format='(g,g)'
		endfor
		free_lun, lun
		
		readcol, 'gon_out/prof_mcmc2', tt,Teff, format=('F,X,X,X,F')
		tt/=(24*3600.0)
		FF = 1.38d-16*Teff/(1.6d-12)
		oplot, tt,FF,linestyle=1,col=80

		openw, lun,'lc_1731_mcmc2.dat', /get_lun
		for i=0,n_elements(tt)-1 do begin
			if tt[i] gt 0.0 then printf, lun, tt[i], FF[i], format='(g,g)'
		endfor
		free_lun, lun

		readcol, 'gon_out/prof_mcmc3', tt,Teff, format=('F,X,X,X,F')
		tt/=(24*3600.0)
		FF = 1.38d-16*Teff/(1.6d-12)
		oplot, tt,FF,linestyle=0,col=250

		openw, lun,'lc_1731_mcmc3.dat', /get_lun
		for i=0,n_elements(tt)-1 do begin
			if tt[i] gt 0.0 then printf, lun, tt[i], FF[i], format='(g,g)'
		endfor
		free_lun, lun
		readcol, 'gon_out/prof_mcmc4', tt,Teff, format=('F,X,X,X,F')
		tt/=(24*3600.0)
		FF = 1.38d-16*Teff/(1.6d-12)
		oplot, tt,FF,linestyle=1,col=250

		openw, lun,'lc_1731_mcmc4.dat', /get_lun
		for i=0,n_elements(tt)-1 do begin
			if tt[i] gt 0.0 then printf, lun, tt[i], FF[i], format='(g,g)'
		endfor
		free_lun, lun


		oplot, [5230.0,5230.0],[10.0,1000.0],linestyle=2

	endif




		if (source eq '1659' or source eq '1731' or source eq 'XTEJ' or source eq '0748' or source eq '0556') and not keyword_set(noextras) then  begin
		; plot lightcurve from simulation
		readcol, 'gon_out/prof_'+source+'_A', tt,Teff, format=('F,X,X,X,F')
		tt/=(24*3600.0)
		; the output is already redshifted
		; but needs to be converted to eV
		FF = 1.38d-16*Teff/(1.6d-12)
		;oplot, tt,FF,linestyle=0
		; plot lightcurve from simulation
		readcol, 'gon_out/prof_'+source+'_B', tt,Teff, format=('F,X,X,X,F')
		tt/=(24*3600.0)
		; the output is already redshifted
		; but needs to be converted to eV
		FF = 1.38d-16*Teff/(1.6d-12)
		;oplot, tt,FF,linestyle=1

;		readcol, 'gon_out/prof_'+source+'_C', tt,Teff, format=('F,X,X,X,F')
;		tt/=(24*3600.0)
		; the output is already redshifted
		; but needs to be converted to eV
;		FF = 1.38d-16*Teff/(1.6d-12)
;		oplot, tt,FF,linestyle=2

		endif

		if (source eq '0556') then xyouts, 2,170,textoidl('MAXI J0556-332')
		if (source eq '1731') then xyouts, 700,120,textoidl('KS 1731-260')
		if (source eq '1659') then xyouts, 800,132,textoidl('MXB 1659-29')
;		if (source eq '1659') then xyouts, 200,152,textoidl('MXB 1659-29')
		if (source eq 'xtej') then begin
			xyouts, 500,170,textoidl('XTE J1701-462')
			oplot,[2900.0,2900.0],[10.0,1000.0],linestyle=1
		endif
		if (source eq 'XTEJ2') then xyouts, 10,170,textoidl('XTE J1709-267')
		if (source eq '0748') then xyouts, 400,130,textoidl('EXO 0748-676')
		if (source eq 'terz' or source eq 'terz2' or source eq 'terz3') then xyouts, 200,125,textoidl('IGR J17480-2446')

		if (source eq 'terz' or source eq 'terz2' or source eq 'terz3') then begin
			oplot, [1.0,1e5],[73.6-1.6,73.6-1.6],linestyle=1
			oplot, [1.0,1e5],[73.6+1.6,73.6+1.6],linestyle=1		
		endif
		

;		t=[2.0,3.5]
;		F=[alog10(140.0)-(t-2.0)/12.0]
;		oplot, 10^t, 10^F, linestyle=2, col=250, thick=2
		t=[1.0,2.5]
		F=alog10(140.0)-0.1*(t-1.0)
		;oplot, 10.0^t, 10.0^F, linestyle=2, col=250, thick=2
		t=[1.0,2.5]
		F=alog10(129.0)-0.09*(t-1.0)
		;oplot, 10.0^t, 10.0^F, linestyle=2, col=120, thick=2


	if keyword_set(ps) then close_ps
end


pro lcplot, namestring, ls, tscal=tscal, linecol=linecol,Lscale=Lscale,Lmin=Lmin,mu=mu
	if (namestring ne '') then namestring='_'+namestring
	if not keyword_set(Lmin) then Lmin = 1d33
;	if keyword_set(Lscale) then begin
;		readcol, 'gon_out/prof'+namestring, tm, F2,Fm, F3,format=('F,F,F,X,F')	
;		tm/=(24.0*3600.0)
;		Fm*=4.0*!dpi*1.12d6^2
;		Fmin=min(Fm)
;		print, namestring, Fmin
;		namestring = namestring + '_mu1'
;	endif
	readcol, 'gon_out/prof'+namestring, tm, F2,Fm, F3,format=('F,F,F,X,F')	
	print, 'Read '+'gon_out/prof'+namestring
	tm/=(24.0*3600.0)
	Fm*=4.0*!dpi*1.12d6^2
;	Fmin=min(Fm)

	if keyword_set(mu) then begin
		sin2theta = 4.0*mu*mu/(1.0+3.0*mu*mu)
		Fm*=sin2theta
		tm/=sin2theta
	endif
	
	if keyword_set(tscal) then begin
 	;	tm/=tscal
	endif



	if keyword_set(Lscale) then Fm=Fm*Lscale+Lmin*(1.0-Lscale)
;	print, min(Fm), Fmin
	;print, tm, Fm

	if keyword_set(linecol) then begin
		oplot, tm, Fm,linestyle=ls, col=linecol
	endif else begin
		oplot, tm, Fm,linestyle=ls
	endelse
	print, 'Plotted ', namestring
	
	openw, lun,'lc'+namestring+'.dat', /get_lun
	for i=0,n_elements(tm)-1 do printf, lun, tm[i], Fm[i], format='(g,g)'
	free_lun, lun
	
end


pro lcsum2, namestring, ls, mumin=mumin, Lscale=Lscale,source=source
	
	if not keyword_set(mumin) then mumin=0
	
	name = 'gon_out/prof_'+namestring
	if (namestring eq '') then name = 'gon_out/prof'
	print, 'Reading file ',name
	readcol, name, tm,Fm,format=('D,X,D'),/silent	
	tm/=(24.0*3600.0)  ; convert time to days
	Fm*=4.0*!dpi*1.12d6^2
	
	if keyword_set(source) then begin
		if (strcmp(source,'1627',4)) then begin
			yr=[1d33,1d35]
			xr=[10.0,10000.0]
		endif
		plot, tm, Fm, /nodata,/xlog,/ylog,yrange=yr,xrange=xr,xtitle=textoidl('Time (days)'),$
						ytitle=textoidl('Luminosity (erg s^{-1})')
		if (strcmp(source,'1627_1998',9)) then begin
			readcol, 'data/1998cooling.dat', t, F, dF, format=('D,D,D')
			F*=1d35
			dF*=1d35
			oploterror, t, F, dF, psym=6, /nohat, symsize=0.7,col=250
			dd=11.0
		;endif
		;if (strcmp(source,'1627_2008',9)) then begin
			readcol, 'data/2008cooling.dat', t, F, dF, format=('D,D,D')
			F*=1d35
			dF*=1d35
			oploterror, t[0:n_elements(t)-2], F[0:n_elements(t)-2], dF[0:n_elements(t)-2], psym=6, /nohat, symsize=0.7
			dd=11.0
			oploterror, t[n_elements(t)-1], F[n_elements(t)-1], dF[n_elements(t)-1], psym=6, /nohat, symsize=0.7,/lobar			
		endif
	endif	
	
	LL = [Fm[0]]
	for i=1L,n_elements(tm)-1 do begin
		II = 0.0
		nmu = 100
		for j=-nmu,nmu do begin
			mu=(1.0/nmu)*j
			; dipole l=1
			mup2 = 4.0*mu*mu/(1.0+3.0*mu*mu)
			; l=2
			;mup2 = (1.0 - 6.0*mu^2 + 9.0*mu^4)/(1.0 - 2.0*mu^2 + 5.0*mu^4)
			if (mu gt mumin) then begin
				II +=  (1.0/nmu) * mup2 * interpol(Fm,tm,tm[i]*mup2)
			endif else begin
				II += (1.0/nmu) * mup2 * Fm[-1]
			endelse
			
		endfor
		LL=[LL,II]
	endfor
	oplot, tm, LL, linestyle=ls
			
end




	
pro lcsum, namestring, ls, Rvec=Rvec, muup=muup, Lscale=Lscale
	
	nt=200
	
	nzones=10
	nstart=1
	
	tvec = dindgen(nt)*3.5/nt
	tvec = 10^tvec

	Fvec = dblarr(nt)
	Rvec = dblarr(nt)
	FF = dblarr(nzones,nt)
	Teff = dblarr(nzones,nt)
		
	nend = nzones
	if keyword_set(muup) then nend=muup*nzones
			
	for i=nstart,nend do begin
		name = 'gon_out/prof_'+namestring+'_mu'+strtrim(string(1.0*i/nzones,format='(f3.1)'))
		print, 'Reading file ',name
		readcol, name, tm,Fm,format=('D,X,D'),/silent	
		tm/=(24.0*3600.0)  ; convert time to days		
		Fvec1 = interpol(Fm,tm,tvec)  ; interpolate the flux onto the regular time grid

		; limb darkening law
		nn=0.0
		limb = (1.0+nn)*(1.0*i/nzones)^nn
		limb=1.0
		; effective temperature at each time for this zone
		Teff[i-1,*] = 1.38d-16*(Fvec1/5.67d-5)^0.25/1.6d-9   ; Teff in keV	
		
		FF[i-1,*] = Fvec1 * (1.0/nzones) * limb   ; (L/4piR^2) from patch i as a function of time

	endfor

	AA = 4.0*!dpi*1.12d6^2
	Fvec = total(FF,1)*AA  ; the luminosity is 4piR^2 * sum over all patches
	if keyword_set(Lscale) then Fvec *= Lscale;
	oplot, tvec, Fvec, linestyle=ls

	if (0) then begin

	fc=1.8

	nspec=1000
	E=(dindgen(nspec)+1.0)*5.0/(nspec+1)
	Epeak=dblarr(nt)
	for j=0,nt-1 do begin
		; loop through each time
		spec = dblarr(nspec)  ; reset spectrum
		for i=nstart,nzones do begin
			; loop over zones
			TT = Teff[i-1,j]   ; effective Temperature at this time and zone
			spec += E^3/(exp(E/(fc*TT))-1.0)   ; the corresponding spectrum with color corr
		endfor
		maxT=max(spec,ind)
		Epeak[j]=E[ind]   ; this is the energy at the peak of the spectrum at this time
	endfor
		
	; estimate the colour temperature from the peak energy
	tempvec=Epeak/2.28
	plot, tvec, tempvec, xrange=[1.0,1000.0],xstyle=1, /xlog,$
		xtitle=textoidl('Days after outburst (d)'), $
		ytitle=textoidl('Estimate of kT from spectral peak'), yrange=[0,1.0]

	xyouts, 300.0, 0.8, textoidl('f_c='+strtrim(string(fc,format='(g4.2)'),2)),charsize=1.2

	Rvec=dblarr(nt)
	for i=0,nt-1 do begin	
;		denom = max(FF[*,i])    ; one zone dominates and sets the temperature
;		meanT = total((FF[*,i]*nzones)^0.25,1)/nzones    ; this says the temperature averages
;		Rvec[i]=Fvec[i]/meanT^4
		Rvec[i]=sqrt(Fvec[i]/tempvec[i]^4)
	endfor
	Rvec/=max(Rvec)
	
	plot, tvec, Rvec, /xlog, xrange=[1.0,1000.0],xstyle=1,$
		xtitle=textoidl('Days after outburst (d)'), $
		ytitle=textoidl('Normalized blackbody radius'), yrange=[min(Rvec)-0.1,1.1]
endif
end







function Loneday, namestring, tt=tt
	
	readcol, 'gon_out/prof_'+namestring, tm, F2,Fm, F3,format=('F,F,F,X,F')	
		
	tm/=(24.0*3600.0)
	Fm*=4.0*!dpi*1.12d6^2
	
	if not keyword_set(tt) then tt=3.0
	L = interpol(Fm,tm,tt)
	
	return, L
	
end



pro tcb, ps=ps
	if keyword_set(ps) then begin
   	set_plot, 'ps'
   	device,filename='rhoE.ps',/color,/encapsul
  	endif

	readcol, 'LB.dat', B,rho,E,Eobs,Tc,name, format=('D,X,D,D,D,D,A'),delimiter=','

	!p.multi=[0,1,1,0,0]
	B*=1d14
	Tc*=1d8
	plot, B, Tc, /xlog, /ylog,psym=1

	if keyword_set(ps) then begin
     	device,/close
     	set_plot,'x'
  	endif
end



pro rhoT, ps=ps
	if keyword_set(ps) then open_ps,'rhoT.ps'
	
	readcol, 'LB.dat', rho,E,Eobs,name,T, format=('X,X,D,D,D,X,A,D'),delimiter=','
	
	rho*=1d10
	T*=1d9
;	E*=1d42
	print, Eobs
	
	plot, rho, T, /xlog, /ylog,/nodata, xrange=[5d9,1d12],xstyle=1,$
			xtitle=textoidl('Heating depth \rho_s (g cm^{-3})'),$
			ytitle=textoidl('Temperature (K)'), yrange=[5d8,8d9],ystyle=1
	for i=0,n_elements(name)-1 do begin
		
		oplot, [rho[i]], [T[i]], psym=5,symsize=1.1
		nametouse=name[i]
		fac=0.95
		if (strcmp(nametouse,'SGR 1627-41 (1998)',18) or $
				strcmp(nametouse,'1E 2259',7)) then fac=1.05
		xyouts, rho[i]*1.1, T[i]*fac, nametouse, charsize=1.0

	endfor
	
	if keyword_set(ps) then close_ps
end


pro rhoE, ps=ps
	if keyword_set(ps) then open_ps,'rhoE.ps'
	
	readcol, 'LB.dat', rho,E,Eobs,name, format=('X,X,D,D,D,X,A'),delimiter=','
	
	rho*=1d10
;	E*=1d42
	print, Eobs
	
	plot, rho, E, /xlog, /ylog,/nodata, xrange=[5d9,1d12],xstyle=1,$
			xtitle=textoidl('Heating depth \rho_s (g cm^{-3})'), yrange=[0.1,3000.0],ystyle=1,$
			ytitle=textoidl('Energy (10^{42} erg s^{-1})')
	for i=0,n_elements(name)-1 do begin
		
		oplot, [rho[i]], [E[i]], psym=5,symsize=1.1
		nametouse=name[i]
		;nametouse=strtrim(string(i),2)
		if (strcmp(nametouse,'SGR 0501',8)) then begin
			xyouts, rho[i]*0.25, E[i]*0.9, nametouse, charsize=1.0
		endif else begin
			if (strcmp(nametouse,'XTE',3)) then begin
				xyouts, rho[i]*0.27, E[i]*0.9, nametouse, charsize=1.0
			endif else begin
				if (i eq 0) then begin
					xyouts, rho[i]*0.16, E[i]*0.9, nametouse, charsize=1.0
				endif else begin
					xyouts, rho[i]*1.1, E[i]*0.9, nametouse, charsize=1.0
				endelse
			
			endelse
		endelse
		if (Eobs[i] gt 0) then begin
			oplot, [rho[i]],[Eobs[i]],psym=1
		;	oplot, [rho[i],rho[i]],[Eobs[i],E[i]],linestyle=1
		;	xyouts, rho[i]*1.05, Eobs[i]*0.9, nametouse, charsize=0.5
		endif
	endfor
	
	
	rho=9.0+3.0*0.01*dindgen(100)
	E=0.22*(10.0^(rho-10.0))^(5.0/3.0)
;	oplot, 10^rho,E*9.0,linestyle=0
	oplot, 10^rho,E*1.0,linestyle=0

	E=140.0*(10.0^(rho-10.0))^(1.0/3.0)
	oplot, 10^rho,1.7*E/E,linestyle=2
	oplot, 10^rho,E,linestyle=1

	E=846.0*(10.0^(rho-10.0))^(4.0/3.0)
	;oplot, 10^rho,E*0.003,linestyle=4
	
	if keyword_set(ps) then close_ps
end


pro LB, ps=ps
	
	if keyword_set(ps) then open_ps,'LB.ps'

	B=[1e13,1.8e13,3e13,5.6e13,1e14,1.8e14,3e14,5.6e14,1e15,1.8e15]
	time=10.0
	L=0.1*[Loneday('B1e+13_E30_Tc5e7_mu1',tt=time),Loneday('B1.8e+13_E30_Tc5e7_mu1',tt=time),$
		Loneday('B3e+13_E30_Tc5e7_mu1',tt=time),Loneday('B5.6e+13_E30_Tc5e7_mu1',tt=time),$
		Loneday('B1e+14_E30_Tc5e7_mu1',tt=time),Loneday('B1.8e+14_E30_Tc5e7_mu1',tt=time),$
		Loneday('B3e+14_E30_Tc5e7_mu1',tt=time),Loneday('B5.6e+14_E30_Tc5e7_mu1',tt=time),$
		Loneday('B1e+15_E30_Tc5e7_mu1',tt=time),$
		Loneday('B1.8e+15_E30_Tc5e7_mu1',tt=time)]

	plot, B, L, /xlog,/ylog, charsize=1.4, xtitle=textoidl('B (G)'),$
			ytitle=textoidl('Luminosity after 10 days (erg s^{-1})'), yrange=[1d34,1.5d36],ystyle=1,$
			xrange=[1d13,5d15],xstyle=1
		oplot, B,10*L
	time=3.0
;	L=[Loneday('B1e13E30.0',tt=time),Loneday('B3e13E30.0',tt=time),Loneday('B1e14E30.0',tt=time),$
;		Loneday('B1.8e14E30.0',tt=time),Loneday('B3e14E30.0',tt=time),Loneday('B1e15E30.0',tt=time),$
;		Loneday('B1.8e15E30.0',tt=time)]
L=0.1*[Loneday('B1e+13_E30_Tc5e7_mu1',tt=time),Loneday('B1.8e+13_E30_Tc5e7_mu1',tt=time),$
	Loneday('B3e+13_E30_Tc5e7_mu1',tt=time),Loneday('B5.6e+13_E30_Tc5e7_mu1',tt=time),$
	Loneday('B1e+14_E30_Tc5e7_mu1',tt=time),Loneday('B1.8e+14_E30_Tc5e7_mu1',tt=time),$
	Loneday('B3e+14_E30_Tc5e7_mu1',tt=time),Loneday('B5.6e+14_E30_Tc5e7_mu1',tt=time),$
	Loneday('B1e+15_E30_Tc5e7_mu1',tt=time),$
	Loneday('B1.8e+15_E30_Tc5e7_mu1',tt=time)]

;	oplot, B,L,linestyle=2


	readcol, 'out/LB.dat', B, L, name, format=('D,D,X,X,X,X,A'),delimiter=','
	B=2*B
	
	cs=0.8
	for i=0,n_elements(name)-1 do begin
		
		oplot, [B[i]*1d14], [L[i]*1d35], psym=1, symsize=1.5
		if (strcmp(name[i],'SGR 1627',8)) then begin
			xyouts, B[i]*1d14*0.3, L[i]*1d35*0.95, name[i], charsize=cs
		endif else begin
			if (strcmp(name[i],'CXOU',4)) then begin
				xyouts, B[i]*1d14*0.22, L[i]*1d35*1.0, name[i], charsize=cs
			endif else begin
				if (strcmp(name[i],'SGR 0501',8)) then begin
					xyouts, B[i]*1d14*0.35, L[i]*1d35*0.95, name[i], charsize=cs
				endif else begin
					xyouts, B[i]*1d14*1.1, L[i]*1d35*0.95, name[i], charsize=cs
				endelse
			endelse
		endelse

	endfor

	xyouts, 1.5e13, 1e35, textoidl('f_{surface}=0.1'), charsize=1.2
	xyouts, 1.5e13, 1e36, textoidl('f_{surface}=1'), charsize=1.2


	;B=13.0+0.1*dindgen(10)*2.5
	;L=35.9-0.88*(B-14.0)
	;oplot, 10^B,10^L,linestyle=1

	if keyword_set(ps) then close_ps
end

pro open_ps, name
	!p.font=0
	set_plot, 'ps'
	device,filename=name,/color,/times,/encapsul;,xsize=18.0,ysize=18.0
	!p.multi=[0,1,1,0,0]
	!p.charsize=1.2
;	!p.position=square()
	!p.thick=3
	!x.thick=3
	!y.thick=3
;	if keyword_set(xmarg) then !x.margin=xmarg else !x.margin=[7,2]
end

pro close_ps
	device,/close
    set_plot,'x'
	!p.font=-1
	!p.thick=1
	!x.thick=1
	!y.thick=1
		
	
end



pro lc_short,ps=ps

	if keyword_set(ps) then begin
			!p.font=0
			set_plot, 'ps'
			device,filename='lc_short.ps',/color,/times,/encapsul,xsize=20.0,ysize=20.0
			!p.multi=[0,3,3,0,0]
			!p.charsize=1.01
			!x.margin=[7,3]
;			!p.position=square()
			!p.thick=2
			!x.thick=2
			!y.thick=2
	endif
	
	cs=0.6
	!p.symsize=0.3
	
	lc,source='fluxes1547',/noplot,cs=cs
	lc,source='fluxes1822',/noplot,cs=cs
	lc,source='0501',/noplot,cs=cs
	lc,source='j1833',/noplot,cs=cs
	lc,source='2259',/noplot,cs=cs
	lc,source='1627_2008',/noplot,cs=cs

	if keyword_set(ps) then close_ps
end

pro lc_long,ps=ps
	

	if keyword_set(ps) then begin
			!p.font=0
			set_plot, 'ps'
			device,filename='lc_long.ps',/color,/times,/encapsul,xsize=20.0,ysize=20.0
			!p.multi=[0,3,3,0,0]
			!p.charsize=1.01
			!x.margin=[7,3]
;			!p.position=square()
			!p.thick=2
			!x.thick=2
			!y.thick=2
	endif
	
	cs=0.6
	!p.symsize=0.3
	
	lc,source='SGR1900',/noplot,cs=cs
	lc,source='1810',/noplot,cs=cs
	lc,source='1627_1998',/noplot,cs=cs
	lc,source='1647',/noplot,cs=cs
	lc,source='1048',/noplot,cs=cs

	if keyword_set(ps) then close_ps
end



pro lc_all,ps=ps

	if keyword_set(ps) then open_ps,'lc_all.ps'
	loadct,39

	lc,source='fluxes1822',/noplot,/nolabel,lcol=1, /all
	lc,source='2259',/noplot,/nolabel,/overplot,lcol=200, /all
	lc,source='fluxes1547',/noplot,/nolabel,/overplot,lcol=25, /all
	lc,source='1810',/noplot,/nolabel,/overplot,lcol=50, /all
	lc,source='j1833',/noplot,/nolabel,/overplot,lcol=75, /all
	lc,source='1647',/noplot,/nolabel,/overplot,lcol=100,/all
	lc,source='1627_2008',/noplot,/nolabel,/overplot,lcol=125, /all
	lc,source='1627_1998',/noplot,/nolabel,/overplot,lcol=150, /all
	lc,source='SGR1900',/noplot,/nolabel,/overplot,lcol=175, /all
	lc,source='1048',/noplot,/nolabel,/overplot,lcol=225, /all
	lc,source='0501',/noplot,/nolabel,/overplot,lcol=250, /all
	
	if keyword_set(ps) then close_ps
	
end

pro lc_slope, ps=ps
	
	if keyword_set(ps) then open_ps, 'lc_slope.ps'
	
	yr=[1d32,3d35]
	xr=[0.3,3000]
	plot, [1.0,1e5], [1d34,1d35], /xlog, /ylog, xtitle=textoidl('Time (d)'), $
		ytitle=textoidl('Luminosity (erg s^{-1})'), linestyle=0, $
		xrange=xr,xstyle=1, yrange=yr, ystyle=1, /nodata

;	lcplot, 'slope1_env', 0, Lscale=0.1,Lmin=PYL(1d8,1d14)
;	lcplot, 'slope2_env', 0, Lscale=0.1,Lmin=PYL(1d8,1d14)
;	lcplot, 'slope3_env', 0, Lscale=0.1,Lmin=PYL(1d8,1d14)
;	lcplot, 'slope4_env', 0, Lscale=0.1,Lmin=PYL(1d8,1d14)
;	lcplot, 'slope5_env', 0, Lscale=0.1,Lmin=PYL(1d8,1d14)

		lcplot, 'slope6_env', 0, Lscale=0.1,Lmin=PYL(5d7,1d14)
		lcplot, 'slope7_env', 0, Lscale=0.1,Lmin=PYL(5d7,1d14)
		lcplot, 'slope8_env', 0, Lscale=0.1,Lmin=PYL(5d7,1d14)
		lcplot, 'slope9_env', 0, Lscale=0.1,Lmin=PYL(5d7,1d14)
		lcplot, 'slope10_env', 0, Lscale=0.1,Lmin=PYL(5d7,1d14)
	;



;	lcplot, '', 1, Lscale=0.1
	
	readcol, 'data/fluxes1822', t0, format=('D')
	readcol, 'data/fluxes1822', t, F, dF,flag, format=('D,X,D,D,X,I'),skipline=1
	t-=t0[0]
	dd=1.6
	F*=4.0*!dpi*(3.086d21*dd)^2
	dF*=4.0*!dpi*(3.086d21*dd)^2
;	oploterror, t, F, dF, psym=1,/nohat, col=250, errcol=250
	ind = where(flag eq 1)			
;	oploterror, t[ind], F[ind], dF[ind], psym=6,/nohat,symsize=1, col=250, errcol=250
	
	if (1) then begin
		plot_anal, 4.0/3.0, 8d32
		plot_anal, 1.0, 5d33
		plot_anal, 0.833, 9d33
		plot_anal, 0.5, 5d34
		plot_anal, 0.0, 3d35


	endif
	
	xyouts, 100,1d35, textoidl('E_{dep}\propto \rho^\alpha'), charsize=1.4
	
	xyouts,2.0,8.5d32, textoidl('4/3'), charsize=1.01
	xyouts,2.0,4d33, textoidl('1'), charsize=1.01
	xyouts,2.0,1d34, textoidl('5/6'), charsize=1.01
	xyouts,2.0,2.5d34, textoidl('1/2'), charsize=1.01
	xyouts,1.7,8d34, textoidl('\alpha=0'), charsize=1.01
;	t = dindgen(10)*0.1*2.0 + 1.0
;	L = 35.0 - 1.5*(t-1.0)
;	t += 1.0
;	ind = where(L gt 33.0 and L lt 34.0)
;;	oplot, 10^t[ind], 10^L[ind], linestyle=2
;	xyouts, 800.0,1d34,textoidl("L\propto t^{-3/2}"),charsize=1.1
	
	
	if keyword_set(ps) then close_ps
	
end


pro plot_anal, beta, L0

	t=0.0+1.0*0.1*dindgen(10)
	L=alog10(L0) + t*0.21*(6.0*beta-5.0)
	oplot, 10^t, 10^L, linestyle=1

end

pro lc_grid_anal, rho, TT
	t = dindgen(10)*0.1*2.0 + 1.0
	L = 35.0 - 1.5*(t-1.0)

	t = t + alog10(rho) - 10.0
	L += 2.2*alog10(TT/2d9)
	oplot, 10^t, 10^L, linestyle=0, col=250 , thick=3
end


pro lc_grid,ps=ps

	if keyword_set(ps) then open_ps, 'lc_grid.ps'

	; set up the plot
	yr=[1d32,1d36]
	xr=[0.08,9e3]
	plot, [0.1,2e4], [1d34,1d36], /xlog, /ylog, xtitle=textoidl('Time (d)'), $
		ytitle=textoidl('Luminosity (erg s^{-1})'), linestyle=0, $
		xrange=xr,xstyle=1, yrange=yr, ystyle=1, /nodata

	;lcplot, 'slope1', 0

	; Plot models
	;lcplot, '1822', 0, Lscale=0.09
	;lcplot, '', 0, linecol=250, Lscale=0.1
	if (1) then begin
;		lcplot, 'B1e14_T2e9_rho1e+10',1, tscal=10.0, Lscale=0.1
;		lcplot, 'B1e14_T2e9_rho3e+09',1, tscal=3.0, Lscale=0.1
;		lcplot, 'B1e14_T2e9_rho3e+10',1, tscal=30.0, Lscale=0.1
;		lcplot, 'B1e14_T2e9_rho1e+11',1, tscal=100.0, Lscale=0.1
;		lcplot, 'B1e14_T2e9_rho3e+11',1, tscal=100.0, Lscale=0.1

;		lcplot, 'new_B1e14_T2e9_rho1e+10',0, tscal=10.0, Lscale=0.1
;		lcplot, 'new_B1e14_T2e9_rho3e+09',0, tscal=3.0, Lscale=0.1
;		lcplot, 'new_B1e14_T2e9_rho3e+10',0, tscal=30.0, Lscale=0.1
;		lcplot, 'new_B1e14_T2e9_rho1e+11',0, tscal=100.0, Lscale=0.1
;		lcplot, 'new_B1e14_T2e9_rho3e+11',0, tscal=100.0, Lscale=0.1

;	lcplot, 'B1e14_T2e9_Tc5e7_rho1e+10',1, tscal=3.0, Lscale=0.1
;	lcplot, 'B1e14_T2e9_Tc5e7_rho3e+10',1, tscal=30.0, Lscale=0.1
;	lcplot, 'B1e14_T2e9_Tc5e7_rho3e+09',1, tscal=30.0, Lscale=0.1
;	lcplot, 'B1e14_T2e9_Tc5e7_rho1e+11',1, tscal=100.0, Lscale=0.1
;	lcplot, 'B1e14_T2e9_Tc5e7_rho3e+11',1, tscal=100.0, Lscale=0.1


	lcplot, 'new_B1e14_T2e9_Tc5e7_rho1e+08_mu1_env',1, tscal=0.01, Lscale=0.1, Lmin=PYL(5d7,1d14)

	if (1) then begin
;	lcplot, 'new_B1e14_T2e9_Tc5e7_rho3e+08_mu1_env',1, tscal=0.03, Lscale=0.1, Lmin=PYL(5d7,1d14)
	lcplot, 'new_B1e14_T2e9_Tc5e7_rho1e+09_mu1_env',1, tscal=0.1, Lscale=0.1, Lmin=PYL(5d7,1d14)
;	lcplot, 'new_B1e14_T2e9_Tc5e7_rho3e+09_mu1_env',1, tscal=0.3, Lscale=0.1, Lmin=PYL(5d7,1d14)
	lcplot, 'new_B1e14_T2e9_Tc5e7_rho1e+10_mu1_env',1, tscal=1.0, Lscale=0.1, Lmin=PYL(5d7,1d14)
;	lcplot, 'new_B1e14_T2e9_Tc5e7_rho3e+10_mu1_env',1, tscal=3.0, Lscale=0.1, Lmin=PYL(5d7,1d14)
	lcplot, 'new_B1e14_T2e9_Tc5e7_rho1e+11_mu1_env',1, tscal=10.0, Lscale=0.1, Lmin=PYL(5d7,1d14)
;	lcplot, 'new_B1e14_T2e9_Tc5e7_rho3e+11_mu1_env',1, tscal=30.0, Lscale=0.1, Lmin=PYL(5d7,1d14)

	lcplot, 'new_B1e14_T1e9_Tc5e7_rho1e+08_mu1_env',0, tscal=0.01, Lscale=0.43, Lmin=PYL(5d7,1d14)
;	lcplot, 'new_B1e14_T1e9_Tc5e7_rho3e+08_mu1_env',0, tscal=0.03, Lscale=0.1, Lmin=PYL(5d7,1d14)
	lcplot, 'new_B1e14_T1e9_Tc5e7_rho1e+09_mu1_env',0, tscal=0.1, Lscale=0.43, Lmin=PYL(5d7,1d14)
;	lcplot, 'new_B1e14_T1e9_Tc5e7_rho3e+09_mu1_env',0, tscal=0.3, Lscale=0.1, Lmin=PYL(5d7,1d14)
	lcplot, 'new_B1e14_T1e9_Tc5e7_rho1e+10_mu1_env',0, tscal=1.0, Lscale=0.43, Lmin=PYL(5d7,1d14)
;	lcplot, 'new_B1e14_T1e9_Tc5e7_rho3e+10_mu1_env',0, tscal=3.0, Lscale=0.1, Lmin=PYL(5d7,1d14)
	lcplot, 'new_B1e14_T1e9_Tc5e7_rho1e+11_mu1_env',0, tscal=10.0, Lscale=0.43, Lmin=PYL(5d7,1d14)
;	lcplot, 'new_B1e14_T1e9_Tc5e7_rho3e+11_mu1_env',0, tscal=30.0, Lscale=0.1, Lmin=PYL(5d7,1d14)
	endif
	
	t = dindgen(10)*0.1*2.0 + 1.0
	L = 35.0 - 1.5*(t-1.0)
	t += 1.3
	ind = where(L gt 33.0 and L lt 34.0)
	oplot, 10^t[ind], 10^L[ind], linestyle=2
	xyouts, 1200.0,1d34,textoidl("L\propto t^{-3/2}"),charsize=1.1
	
;	lc_grid_anal, 1d8, 2d9
;	lc_grid_anal, 1d9, 2d9
;	lc_grid_anal, 1d10, 2d9
;	lc_grid_anal, 1d11, 2d9
;	lc_grid_anal, 1d8, 1d9
;	lc_grid_anal, 1d9, 1d9
;	lc_grid_anal, 1d10, 1d9
;	lc_grid_anal, 1d11, 1d9

;	lcplot, 'new_B1e14_T2e9_Tc1e8_rho3e+09_mu1_env',2, tscal=3.0, Lscale=0.1, Lmin=PYL(1d8,1d14)
;	lcplot, 'new_B1e14_T2e9_Tc1e8_rho1e+10_mu1_env',2, tscal=10.0, Lscale=0.1, Lmin=PYL(1d8,1d14)
;	lcplot, 'new_B1e14_T2e9_Tc1e8_rho3e+10_mu1_env',2, tscal=30.0, Lscale=0.1, Lmin=PYL(1d8,1d14)
;	lcplot, 'new_B1e14_T2e9_Tc1e8_rho1e+11_mu1_env',2, tscal=100.0, Lscale=0.1, Lmin=PYL(1e8,1d14)
;	lcplot, 'new_B1e14_T2e9_Tc1e8_rho3e+11_mu1_env',2, tscal=300.0, Lscale=0.1, Lmin=PYL(1e8,1d14)


	endif
			
	xyouts,100.0,3d35,textoidl('log \rho_c=8, 9, 10, 11'),charsize=1.1

	; Plot data
	readcol, 'data/1647', t, F, dF, format=('D,X,X,D,D')
	F*=1d35
	dF*=1d35
	dd=5.0
	Lq=7e32
;	oploterror, t, F, dF, psym=1, /nohat, symsize=0.7, col=120, errcol=120
	
	readcol, 'data/fluxes1822', t0, format=('D')
	readcol, 'data/fluxes1822', t, F, dF,flag, format=('D,X,D,D,X,I'),skipline=1
	t-=t0[0]
	dd=1.6
	F*=4.0*!dpi*(3.086d21*dd)^2
	dF*=4.0*!dpi*(3.086d21*dd)^2
	;oploterror, t, F, dF, psym=1,/nohat, col=250, errcol=250
	ind = where(flag eq 1)			
;	oploterror, t[ind], F[ind], dF[ind], psym=6,/nohat,symsize=1, col=250, errcol=250
	
	; 1627
	dd=11.0
	readcol, 'data/1998cooling.dat', t, F, dF, format=('D,D,D')
	F*=1d35
	dF*=1d35
;	oploterror, t, F, dF, psym=4, /nohat, col=80,errcol=80
	readcol, 'data/2008cooling.dat', t, F, dF, format=('D,D,D')
	F*=1d35
	dF*=1d35
;	oploterror, t, F, dF, psym=6, /nohat, col=100,errcol=100

	
	if keyword_set(ps) then close_ps

end




pro lc, source=source,ps=ps, nodata=nodata, nolabel=nolabel, noplot=noplot, overplot=overplot, lcol=lcol,$
			cs=cs, all=all, Lscale=Lscale, Lmin=Lmin,zoom=zoom
	
;	!p.multi=[0,1,3,0,0]
;	!p.charsize=2.5
	!p.multi=[0,1,1,0,0]
		!p.charsize=1.4
		
	if not keyword_set(source) then source = ''
		
	if keyword_set(ps) then open_ps, 'lc'+source+'.ps'
;		if keyword_set(ps) then begin
;		if keyword_set(source) then name = 'lc'+source+'.ps' else name='lc.ps'
   ;		open_ps, name
;		set_plot,'ps'
;		device, filename=name, /encapsul,/color
;		!p.thick=3
;		!x.thick=3
;		!y.thick=3
;		!p.charthick=3
 ; 	endif
	
	yr=[1d32,1d37]
	xr=[0.1,10000.0]
	if (strcmp(source,'0556',4)) then begin
		yr=[1d34,5d35]
		xr=[1.0,5000.0]
	endif
	if (strcmp(source,'2259',4)) then yr=[1d34,1d35]
	if (strcmp(source,'fluxes1822',10)) then begin
		yr=[3d30,1.5d35]
		xr=[0.5,3000]
		if keyword_set(zoom) then begin
			yr=[3d30,2d33]
			xr=[100.0,3000]
		endif
	endif
	if (strcmp(source,'fluxes1547',10)) then yr=[1d33,1d36]
	if (strcmp(source,'1659',4)) then yr=[1d31,1d36]
	if (strcmp(source,'0501',4)) then begin
		yr=[1d34,7d35]
		xr=[1.0,6000.0]
	endif
	if (strcmp(source,'j1833',5)) then yr=[7d33,2d35]
	if (strcmp(source,'1810',4)) then begin
		yr=[5d32,5d35]
		xr=[10.0,6000.0]
	endif
	if (strcmp(source,'1048',4)) then begin
		yr=[1d32,1d35]
		xr=[100.0,10000.0]
	endif
	if (strcmp(source,'1647',4)) then begin
		yr=[3d32,4d35]
		xr=[0.3,6000.0]
	endif
	if (strcmp(source,'1627_1998',9)) then begin
		yr=[3d32,1d35]
		xr=[10.0,10000.0]
	endif
	if (strcmp(source,'1627_2008',9)) then begin
		yr=[8d32,1d35]
		xr=[1.0,10000.0]
	endif
	if (strcmp(source,'1900',4)) then yr=[1d34,1d38]
	if (strcmp(source,'SGR1900',7)) then begin
		 yr=[1d34,1d38]
		xr=[0.1,10000.0]
	endif

	if (strcmp(source,'1731',4)) then yr=[1d32,3d35]

	if keyword_set(all) then yr=[2d32,3d36]
	
	Lq=0.0
	
	; F= flux at grid point 1
	; F2= flux at grid point 2
	; F3= flux computed using Teff-Tb relation
	readcol, 'gon_out/prof', tm, F2, Fm,F3,Fin,Lnu,format=('F,F,D,X,F,X,X,X,D,D')	
	tm/=(24.0*3600.0)

	Fm*=4.0*!dpi*1.12d6^2
	F2*=4.0*!dpi*1.12d6^2
	Fin*=4.0*!dpi*1.12d6^2
	
	;if keyword_set(Lscale) then Fm=Fm*Lscale+Lmin
	if keyword_set(Lscale) then Fm=Fm*Lscale+Lmin*(1.0-Lscale)
	
;		plot, t, F, /xlog, /ylog, xtitle=textoidl('Time (d)'), $
;				ytitle=textoidl('Luminosity (erg s^{-2})'), linestyle=0, charsize=1.5, $
;				xrange=[0.01,1000.0],xstyle=1, yrange=[1d35,1d38]
	if not keyword_set(overplot) then begin
		plot, tm, Fm, /xlog,/ylog, xtitle=textoidl('Time (d)'), $
			;xtitle=textoidl('Days since BAT trigger (d)'), $
			ytitle=textoidl('Luminosity (erg s^{-1})'), linestyle=0, $
			xrange=xr,xstyle=1, yrange=yr, ystyle=1, nodata=noplot
	endif
;	oplot, t,F,linestyle=1

	;if not keyword_set(noplot) then oplot, tm, abs(Fin), col=80
	;if not keyword_set(noplot) then oplot, tm, abs(Lnu), col=180

	ind = where(Fin lt 0.0)
	Lin=0.0
	Lout=0.0
	for i=1,n_elements(ind)-1 do begin
		j = ind[i]
		Lout += F2[j]*(tm[j]-tm[j-1])
		Lin += Fin[j]*(tm[j]-tm[j-1])		
	endfor
	print, "Lin=", Lin * 3600.0*24.0
	print, "Lout=", Lout * 3600.0*24.0
	
	Lout = 0.0
	for i=1,n_elements(tm)-1 do begin
		Lout += (F2[i]-Fin[i])*(tm[i]-tm[i-1])* 3600.0*24.0	
	endfor
	print, "Total energy leaving crust = ", Lout
	
	
; PL
;	L=2.5d35/(tm/10.0)^0.5
;	oplot, tm, L, linestyle=2



;	readcol, 'out/nu.dat', rhonu, TTnu, format=('D,D')
;	tnu = (rhonu/3d9)
;	Lnu = 1.7d35 * 0.1 * (TTnu/10.0)^2.1
;	oplot, tnu, Lnu, linestyle=1


;oplot, tm, 5d33 * (tm)^(7.0*21.0/100.0), col=120
;oplot, tm, 8d35 * (tm)^(-3.0/4.0), col=120
;oplot, tm, 2d34 * (tm)^(63.0/140.0), col=120

if not keyword_set(all) then begin



if (0) then begin
lcplot, 'B15E3.0Q1',1
lcplot, 'B14E3.0Q1',1
lcplot, 'B13E3.0Q1',1
lcplot, 'B16E3.0Q1',1
xyouts, 80.0, 1d35, textoidl('B=10^{13},10^{14},10^{15},10^{16}G')
endif

if (0) then begin

	lcplot, 'B15E10.0Q1Tc4',3
	lcplot, 'B15E10.0Q1Tc1',0
	lcplot, 'B15E10.0Q1Tc2',2
	lcplot, 'B15E10.0Q1Tc5',1
	
endif

if (0) then begin

	lcplot, 'B1e14E30.0',0
	lcplot, 'B3e14E30.0',0
	lcplot, 'B1e15E30.0',0
	lcplot, 'B3e15E30.0',0


endif


if (0) then begin
;lcplot, 'maxL_15', 0, Lscale=0.1, Lmin = PYL(5e7,1e15), linecol=250
lcplot, 'maxL_14', 0, Lscale=0.1, Lmin = PYL(5e7,1e14), linecol=80
;lcplot, 'maxL_3e14', 0, Lscale=0.1, Lmin = PYL(5e7,3e14), linecol=120
;lcplot, 'maxL_3e15', 0, Lscale=0.1, Lmin = PYL(5e7,3e15), linecol=250
;lcplot, 'maxL_0', 0, Lscale=0.1, Lmin = PYL(5e7,0)

xyouts, 0.3,1d35, textoidl('B=0')
xyouts, 0.3,10^34.8, textoidl('B=10^{14}G'), col=80
xyouts, 0.3,10^34.6, textoidl('B=10^{15}G'), col=250
xyouts, 0.3,10^34.4, textoidl('B=3\times 10^{15}G'), col=120
;xyouts, 0.3,10^34.2, textoidl('B=3\times 10^{15}G'), col=250

readcol, 'out/fcontour_0.dat', Fcont, tcont, format=('X,X,X,D,X,X,D')
tcont/=24.0*3600.0
Fcont=10^Fcont
Fcont*=4.0*!dpi*1.2d6^2/(1.24^2)
;oplot,  tcont*1.24/2, 0.5*(1.6/2.28)*Fcont*0.1, linestyle=2

readcol, 'out/fcontour_14.dat', Fcont, tcont, format=('X,X,X,D,X,X,D')
tcont/=24.0*3600.0
Fcont=10^Fcont
Fcont*=4.0*!dpi*1.2d6^2/(1.24^2)
oplot,  tcont*1.24/2, 0.5*(1.6/2.28)*Fcont*0.1, linestyle=2, col=80

endif


;lcplot, 'temp_1',1
;lcplot, 'temp_2',2
;lcplot, 'temp_3',3
if (0) then begin
;	lcplot, 'B1e15E3.0_1e9',2
;	lcplot, 'B1e15E1.0_1e9',2
;	lcplot, 'B1e15E0.3_1e9',2
;	lcplot, 'B1e15E10.0_1e9',2
;	lcplot, 'B1e15E100.0_1e9',2
;	lcplot, 'B1e15E30.0_1e9',2
;lcplot, 'B1e14E3.0_1e9',0
;lcplot, 'B1e14E1.0_1e9',0
;lcplot, 'B1e14E0.3_1e9',0
;lcplot, 'B1e14E10.0_1e9',0
;lcplot, 'B1e14E100.0_1e9',0
;lcplot, 'B1e14E30.0_1e9',0

if (0) then begin
areaf=0.1
lfloor=PYL(1e8,1e14)
lcplot, 'B1e14E3_1e9_mu1_env',0,Lscale=areaf,Lmin=lfloor
lcplot, 'B1e14E1_1e9_mu1_env',0,Lscale=areaf,Lmin=lfloor
lcplot, 'B1e14E0.3_1e9_mu1_env',0,Lscale=areaf,Lmin=lfloor
lcplot, 'B1e14E10_1e9_mu1_env',0,Lscale=areaf,Lmin=lfloor
lcplot, 'B1e14E100_1e9_mu1_env',0,Lscale=areaf,Lmin=lfloor
lcplot, 'B1e14E30_1e9_mu1_env',0,Lscale=areaf,Lmin=lfloor
endif

if (0) then begin
lcplot, 'B1e15E3_1e9',2,Lscale=areaf
lcplot, 'B1e15E1_1e9',2,Lscale=areaf
lcplot, 'B1e15E0.3_1e9',2,Lscale=areaf
lcplot, 'B1e15E10_1e9',2,Lscale=areaf
lcplot, 'B1e15E100_1e9',2,Lscale=areaf
lcplot, 'B1e15E30_1e9',2,Lscale=areaf
endif

;lcplot, '',2,Lscale=0.1

;lcplot, 'B1e14E3_1e9_new2',2
;lcplot, 'B1e14E1_1e9_new2',2
;lcplot, 'B1e14E0.3_1e9_new2',2
;lcplot, 'B1e14E10_1e9_new2',2
;lcplot, 'B1e14E100_1e9_new2',2
;lcplot, 'B1e14E30_1e9_new2',2

xyouts, 40.0,10^(36.4), textoidl('E_{25}=0.3, 1, 3, 10, 30, 100'),charsize=1.1
;xyouts, 40.0,10^(36.7), textoidl('B=10^{14} G'),charsize=1.1
;xyouts, 40.0,10^(36.1), textoidl('10% of the surface'),charsize=1.1
;xyouts, 40.0,10^(35.8), textoidl('E_{total}=5\times 10^{41} ergs E_{25}'),charsize=1.1
;xyouts, 80.0,4d35, textoidl('B=10^{14},10^{15} G'),charsize=1.2
endif

if (0) then begin
t = dindgen(10)*0.1*2.0
L = 34.5 - 1.5*(t-1.0)
t+=0.2
ind = where(L gt 33.0 and L lt 35.0)
oplot, 10^t[ind], 10^L[ind], linestyle=2
xyouts, 300.0,3d34,textoidl("L\propto t^{-1.05}"),charsize=1.1
endif
;plot_anal, 0.0, 5d35



if (0) then begin
	lcplot, 'B1e15E3.0',2
	lcplot, 'B1e15E1.0',2
	lcplot, 'B1e15E0.3',2
	lcplot, 'B1e15E10.0',2
	lcplot, 'B1e15E100.0',2
	lcplot, 'B1e15E30.0',2
lcplot, 'B1e14E3.0',0
lcplot, 'B1e14E1.0',0
lcplot, 'B1e14E0.3',0
lcplot, 'B1e14E10.0',0
lcplot, 'B1e14E100.0',0
lcplot, 'B1e14E30.0',0
xyouts, 80.0,4d35, textoidl('E_{25}=0.3,1,3,10,30,100'),charsize=1.2
xyouts, 80.0,2d35, textoidl('B=10^{14},10^{15} G'),charsize=1.2
endif
if (0) then begin
lcplot, 'B15E3.0Q1',0
lcplot, 'B15E1.0Q1',0
lcplot, 'B15E0.3Q1',0
lcplot, 'B15E10.0Q1',0
lcplot, 'B15E100.0Q1',0
lcplot, 'B15E300.0Q1',0
xyouts, 80.0,6d35, textoidl('E_{25}=0.3,1,3,10,100,300')
endif

if (0) then begin
	lcplot, 'deep1',0, Lscale=0.1, Lmin=PYL(3e7,1d14)  ; no sph, Tc=3e7  gap 1  Q=1
	lcplot, 'deep2',2, Lscale=0.1, Lmin=PYL(3e7,1d14)  ; no sph, Tc=3e7  gap 5
	lcplot, 'deep3',1, Lscale=0.1, Lmin=PYL(3e7,1d14)  ; no sph, Tc=3e7  gap 4
	lcplot, 'deep4',3, Lscale=0.1, Lmin=PYL(3e7,1d14)  ; no sph, Tc=3e7  gap 0

	lcplot, 'deep5',0, Lscale=0.1, Lmin=PYL(3e7,1d14),linecol=250  ; sph, Tc=3e7  gap 1
	lcplot, 'deep6',2, Lscale=0.1, Lmin=PYL(3e7,1d14),linecol=250  ; sph, Tc=3e7  gap 5
	lcplot, 'deep7',1, Lscale=0.1, Lmin=PYL(3e7,1d14),linecol=250  ; sph, Tc=3e7  gap 4
	lcplot, 'deep8',3, Lscale=0.1, Lmin=PYL(3e7,1d14),linecol=250  ; sph, Tc=3e7  gap 0

	lcplot, 'deep9',0, Lscale=0.1, Lmin=PYL(3e7,1d14),linecol=80  ; no sph, Tc=3e7  gap 1  Q=0
	lcplot, 'deep10',0, Lscale=0.1, Lmin=PYL(3e7,1d14),linecol=180  ; no sph, Tc=3e7  gap 1  Q=10


	lcplot, 'deep11',0, Lscale=0.1, Lmin=PYL(1e8,1d14),linecol=80  ; no sph, Tc=3e7  gap 1  Q=0
	lcplot, 'deep12',0, Lscale=0.1, Lmin=PYL(1e8,1d14)  ; no sph, Tc=3e7  gap 1  Q=1
	lcplot, 'deep13',0, Lscale=0.1, Lmin=PYL(1e8,1d14),linecol=180  ; no sph, Tc=3e7  gap 1  Q=10
	
	
	

;	xyouts,500.0,2d32,textoidl('with'),charsize=1.2
;	xyouts,500.0,1.2d32,textoidl('SF phonons'),charsize=1.2
;	xyouts,2000.0,1d33,textoidl('SF phonons'),charsize=1.2
;	xyouts,2000.0,1.8d33,textoidl('without'),charsize=1.2

;xyouts, 1500.0, 8d33, textoidl('Q_{imp}=0,1,10')

;xyouts,2500.0,1.2d32,textoidl('without SF phonons'),charsize=0.9
;xyouts,400.0,1d32,textoidl('with SF phonons'),charsize=0.9

endif


if (0) then begin
	;lcplot, 'B14rho1',0,tscal=3.0
	;lcplot, 'B14rho2',0,tscal=10.0
	;lcplot, 'B14rho3',0,tscal=30.0
	;lcplot, 'B14rho4',0,tscal=100.0
;	lcplot, 'B14rho5',2,tscal=10.0
	;lcplot, 'B14rho6',2,tscal=3.0
	;lcplot, 'B14rho7',2,tscal=30.0
;	lcplot, 'B14rho8',2,tscal=100.0
	
;	lcplot, 'B1e15E30.0_1e9',4

	lcplot, 'B1e14E3.0_1e9',3
	lcplot, 'B1e14E3.0_1e9_constantT4.2e8',0
	lcplot, 'B1e14E3.0_1e9_E10.8constantpergram',2
	
;	xyouts,30.0,1.3d35,textoidl('log \rho_c=9.5,10,10.5,11'),charsize=1.0
endif


if (0) then begin
;	lcplot, 'B14rho1',0,tscal=3.0
;	lcplot, 'B14rho2',0,tscal=10.0
;	lcplot, 'B14rho3',0,tscal=30.0
;	lcplot, 'B14rho4',0,tscal=100.0
;;	lcplot, 'B14rho5',2,tscal=10.0
;	lcplot, 'B14rho6',2,tscal=3.0
;	lcplot, 'B14rho7',2,tscal=30.0
;	lcplot, 'B14rho8',2,tscal=100.0
	
	
if (0) then begin
	lcplot, 'B1e14_T1e9_rho1e+10',0,linecol=250, tscal=10.0, Lscale=0.1
	lcplot, 'B1e14_T1e9_rho3e+09',0,linecol=250, tscal=3.0, Lscale=0.1
	lcplot, 'B1e14_T1e9_rho3e+10',0,linecol=250, tscal=30.0, Lscale=0.1
	lcplot, 'B1e14_T1e9_rho1e+11',0,linecol=250, tscal=100.0, Lscale=0.1
endif

if (0) then begin
		lcplot, 'B1e14_T2e9_rho1e+10',0, tscal=10.0, Lscale=0.1
		lcplot, 'B1e14_T2e9_rho3e+09',0, tscal=3.0, Lscale=0.1
		lcplot, 'B1e14_T2e9_rho3e+10',0, tscal=30.0, Lscale=0.1
		lcplot, 'B1e14_T2e9_rho1e+11',0, tscal=100.0, Lscale=0.1
		lcplot, 'B1e14_T2e9_rho3e+11',0, tscal=100.0, Lscale=0.1
			lcplot, 'B1e14_T2e9_Tc5e7_rho1e+10',1, tscal=10.0, Lscale=0.1
			lcplot, 'B1e14_T2e9_Tc5e7_rho3e+09',1, tscal=3.0, Lscale=0.1
			lcplot, 'B1e14_T2e9_Tc5e7_rho3e+10',1, tscal=30.0, Lscale=0.1
			lcplot, 'B1e14_T2e9_Tc5e7_rho1e+11',1, tscal=100.0, Lscale=0.1
			lcplot, 'B1e14_T2e9_Tc5e7_rho3e+11',1, tscal=100.0, Lscale=0.1
		
	endif
	
;	lcplot, 'B1e14E3000000000',0,linecol=250,Lscale=0.1
;	lcplot, 'B1e14E10000000000',0,linecol=250,Lscale=0.1
;	lcplot, 'B1e14E30000000000',0,linecol=250,Lscale=0.1
;	lcplot, 'B1e14E100000000000',0,linecol=250,Lscale=0.1
	
	
	xyouts,30.0,1d35,textoidl('log \rho_c=9.5,10,10.5,11'),charsize=1.2
	
;lcplot, 'B15E3.0Q1_rho1',0 ;3e9-1e10    1.7d42
;lcplot, 'B15E0.3Q1_rho1',0 ;3e9-1e10    1.7d41
;lcplot, 'B15E10.0Q1_rho2',1 ;1e10-3e10   8d42
;lcplot, 'B15E1.0Q1_rho2',1 ;1e10-3e10   8d41
;lcplot, 'B15E100.0Q1_rho4',2  ;1e11-3e11  1.1d44
;lcplot, 'B15E10.0Q1_rho4',2  ;1e11-3e11  1.2d43
;lcplot, 'B15E1000.0Q1_rho5',3  ;1e12-3e12  1.0d45
;lcplot, 'B15E100.0Q1_rho5',3 ;1e12-3e12  1.0d44
endif

;lcplot, 'B15E1.0Q1_rho1',1 ;3e9-1e10    5.8e41
;lcplot, 'B15E1.0Q1_rho2',1  ;1e10-3e10  8.1e41ergs
;lcplot, 'B15E1.0Q1_rho3',1  ;3e10-1e11  1.1e42ergs
;lcplot, 'B15E1.0Q1_rho4',1  ;1e11-3e11  1.2e42ergs
;lcplot, 'B15E10.0Q1_rho4',2  ;1e11-3e11  1.2e42ergs
;lcplot, 'B15E1.0Q1_rho5',1  ;1e12-3e12  1.0e42ergs


;	lcplot, '1627', 2
;	lcplot, '1627_2', 1

if (strcmp(source,'1627_1998',9)) then begin
;	lcplot, '1627_1998_step', 0

;	lcplot, '1627_4', 1
;	lcplot, '1627_8', 2
;	lcplot, '1627_7', 0
;	lcplot, '1627_N100', 0

lcplot, '1627_1998_new', 0
lcplot, '1627_1998_new_lowTc_all', 1
lcplot, '1627_1998_new_lowTc', 2

endif
if (strcmp(source,'1627_2008',9)) then begin
	;	lcplot, '1627_9', 0
;		lcplot, '1627_2008_step', 0
;	lcplot, '1627_10', 2

;lcplot, '1627_new_2', 0
;lcplot, '1627_new_1', 2
;lcplot, '1627_new_3', 1

lcsum2, '1627_chandra1',0,mumin=0.78
lcsum2, '1627_chandra2',1,mumin=0.58


oplot, [2700.0,2700.0],[1d32,1d36],linestyle=2

endif
	
	
if (strcmp(source,'0501',4)) then begin
	lcplot, '0501_step',0	
endif
		
if (strcmp(source,'1647',4)) then begin
;	lcplot, '1647_step', 0
;	lcplot, '1647_1', 0
lcplot, '1647_new_1', 2
lcplot, '1647_new_2', 1
lcplot, '1647_new_3', 0
;	lcplot, '1647_2', 1
;	lcplot, '1647_3', 2
;	lcplot, '1647_4', 2

;lcplot, '1647_mu1.0', 1
;lcplot, '1647_mu0.9', 1
;lcplot, '1647_mu0.8', 1
;lcplot, '1647_mu0.7', 1
;l;cplot, '1647_mu0.6', 1
;lcplot, '1647_mu0.5', 1
;lcplot, '1647_mu0.4', 1
;lcplot, '1647_mu0.3', 1
;lcplot, '1647_mu0.2', 1
;lcplot, '1647_mu0.1', 1

endif

		
;	lcplot, '1627_3', 1
;	lcplot, '1627_6', 3
	
;	lcplot, '1900_noneut', 1   ; turned off neutrinos, Q=0,B=1e15,Tc=7e8,E25=17,
					; grid has yt=1e10, heats everywhere
	
	
if (0) then begin
	lcplot, 'T1e9_y1e10_mu1.0',1
	lcplot, 'T1e9_y1e10_mu1.0',0,Lscale=0.1,Lmin=PYL(1e8,1d14) 
	lcsum2, 'T1e9_y1e10', 2

	y=35.7
	oplot, [100.0,200.0],[10^y,10^y], linestyle=0
	xyouts, 250.0, 10^(y-0.01), textoidl('Radial magnetic field'), charsize=1.01
	xyouts, 250.0, 10^(y-0.15), textoidl('(f_{surf}=0.1)'), charsize=1.01
	y=35.35
	oplot, [100.0,200.0],[10^y,10^y], linestyle=1
	xyouts, 250.0, 10^(y-0.01), textoidl('Radial magnetic field'), charsize=1.01
	xyouts, 250.0, 10^(y-0.15), textoidl('(f_{surf}=1)'), charsize=1.01
	y=35.0
	oplot, [100.0,200.0],[10^y,10^y], linestyle=2
	xyouts, 250.0, 10^(y-0.01), textoidl('Dipole magnetic field'), charsize=1.01
	xyouts, 250.0, 10^(y-0.15), textoidl('(f_{surf}=1)'), charsize=1.01

endif
	
if (strcmp(source,'fluxes1822',10)) then begin		
;	lcsum2,'1822_new_B', 2  ; 2e9 to 3e10  E25=2.3  Q=1   B=6e13   Tc=2e7
;	lcplot, '1822_new_C',0,Lscale=0.15,Lmin=PYL(1.5d7,1d14) 
;	lcplot, '1822_new_G',1,Lscale=0.15,Lmin=PYL(4d7,1d14) 
;	lcplot, '1822_new_F',0,Lscale=0.15,Lmin=PYL(9d7,1d14) 
	
;	lcsum2,'1822_new_D', 2  ; 2e9 to 3e10  E25=2.3  Q=1   B=6e13   Tc=2e7
;	lcsum2,'1822_new_E', 2  ; 2e9 to 3e10  E25=2.3  Q=1   B=6e13   Tc=2e7
		
	lcplot, '1822_chandra1',0,Lscale=0.15,Lmin=PYL(1d7,1d14) 
	lcplot, '1822_chandra2',1,Lscale=0.15,Lmin=PYL(2d7,1d14) 
		
;	arrow,60,6.8d32,90,6.8d32,/data
;	xyouts,27,6.2d32,'ROSAT',charsize=0.9
	arrow,60,6.8d32/20.0,90,6.8d32/20.0,/data
	xyouts,27,6.2d32/20.0,'ROSAT',charsize=0.9


;	lcsum2,'temp',0
	;	lcplot, '1822_A', 1  ; 2e9 to 3e10  E25=2.3  Q=1   B=6e13   Tc=2e7
;	lcplot, '1822_B', 0	   ; 2e9 to 3e10  E25=2.3  Q=10   B=6e13   Tc=2e7
;	lcplot, '1822_C', 3   ; 2e9 to 1e10  E25=1.7  Q=1   B=1e15   Tc=1.5e7
;	lcplot, '1822_D', 1   ; 3.5e9 to 5.5e9  E25=11  Q=10   B=1e15   Tc=1.5e7
;lcplot, '1822_step', 0  
;lcplot, '1822_edep', 2  

if (1) then begin
;lcplot, '1822_new_withB_mu0.9',0
;lcplot, '1822_new_withB_mu0.8',0
;lcplot, '1822_new_withB_mu0.7',0
;lcplot, '1822_new_withB_mu0.6',0
;lcplot, '1822_new_withB_mu0.5',0
;lcplot, '1822_new_withB_mu0.4',0
;lcplot, '1822_new_withB_mu0.3',0
;lcplot, '1822_new_withB_mu0.2',0
;lcplot, '1822_new_withB_mu0.1',0

;lcplot, '1822_new_withB_mu1.0',1,mu=0.9
;lcplot, '1822_new_withB_mu1.0',1,mu=0.8
;lcplot, '1822_new_withB_mu1.0',1,mu=0.7
;lcplot, '1822_new_withB_mu1.0',1,mu=0.6
;lcplot, '1822_new_withB_mu1.0',1,mu=0.5
;lcplot, '1822_new_withB_mu1.0',1,mu=0.4
;lcplot, '1822_new_withB_mu1.0',1,mu=0.3
;lcplot, '1822_new_withB_mu1.0',1,mu=0.2
;lcplot, '1822_new_withB_mu1.0',1,mu=0.1

;lcplot, '1822_new_withB_Tc2e7_mu1.0',0
;lcsum, '1822_new_withB_Tc2e7', 2
;lcsum2, '1822_new_withB_Tc2e7', 0


;lcplot, '1822_new_withB_mu1.0',0
;lcsum, '1822_new_withB', 2
;lcsum2, '1822_new_withB', 0

;lcplot, 'deep1',0, Lscale=1.0, Lmin=PYL(3e7,1d14)
;lcplot, 'deep12',0, Lscale=1.0, Lmin=PYL(1e8,1d14)

;lcplot, '1822_new_step_Tc2e7',linecol=250


endif else begin
	
	if (0) then begin
lcplot, '1822_new_withB_Tc2e7_mu1.0',1, linecol=80
lcplot, '1822_new_withB_Tc2e7_mu0.9',1, linecol=80
lcplot, '1822_new_withB_Tc2e7_mu0.8',1, linecol=80
lcplot, '1822_new_withB_Tc2e7_mu0.7',1, linecol=80
lcplot, '1822_new_withB_Tc2e7_mu0.6',1, linecol=80
lcplot, '1822_new_withB_Tc2e7_mu0.5',1, linecol=80
lcplot, '1822_new_withB_Tc2e7_mu0.4',1, linecol=80
lcplot, '1822_new_withB_Tc2e7_mu0.3',1, linecol=80
lcplot, '1822_new_withB_Tc2e7_mu0.2',1, linecol=80
lcplot, '1822_new_withB_Tc2e7_mu0.1',1, linecol=80
lcsum, '1822_new_withB_Tc2e7', 2, muup=1.0, Lscale=1.0
lcplot, '1822_new_step_Tc2e7'

endif else begin
;	for i=1,10 do begin
;		name='rn2_mu'+strtrim(string(0.1*i,format='(F3.1)'),1)
;  for the case of steps of 0.05 in mu, need to do odd and even separately
;		if (i mod 2 eq 0) then begin
;			name='rn_mu'+strtrim(string(0.1*i,format='(F3.1)'),1)
;		endif else begin
;			name='rn_mu'+strtrim(string(0.05*i,format='(F4.2)'),1)
;		endelse
;		print, 'Opening file ',name
;		lcplot, name,1, linecol=80
;	endfor
;	lcsum, 'rn2',2
;	lcplot, 'rn2'
endelse
endelse
;xyouts, 4.0, 1d35, textoidl('\mu=1 (magnetic pole)'), col=80, charsize=1.2
;xyouts, 4.0, 3d33, textoidl('\mu=0.1 (near equator)'), col=80, charsize=1.2

;lcplot, '1822_mu0.9',1
;lcplot, '1822_mu0.8',1
;lcplot, '1822_mu0.7',1
;lcplot, '1822_mu0.6',1
;lcplot, '1822_mu0.5',1
;lcplot, '1822_mu0.4',1
;lcplot, '1822_mu0.3',1
;lcplot, '1822_mu0.2',1
;lcplot, '1822_mu0.1',1
	
endif

if (strcmp(source,'j1833',5)) then begin
	lcplot, 'j1833_step', 0  
endif
if (strcmp(source,'2259',4)) then begin
	lcplot, '2259_step', 0  
endif

if (strcmp(source,'0556',4)) then begin
	lcplot, '0556_A', 0  
endif
if (strcmp(source,'1810',4)) then begin
	lcplot, '1810_step', 0  
	lcplot, '1810_step_nosph', 2
	lcplot, '1810_step_highg', 1
endif
if (strcmp(source,'1048',4)) then begin
	lcplot, '1048_step', 0  
endif	
if (strcmp(source,'fluxes1547',10)) then begin
	lcplot, '1547_step', 0  
	lcplot, '1547_2_step', 0  
endif
;if (strcmp(source,'SGR1900',7)) then begin
;	lcplot, '1900_step', 0  
;endif

endif
if (0) then begin
;if not keyword_set(overplot) then begin
	lcplot, '1547_step', 0  
	lcplot, '1822_step', 0  
	lcplot, '0501_step', 0  
	lcplot, 'j1833_step', 0  
	lcplot, '2259_step', 0  
	lcplot, '1627_2008_step', 0  
	lcplot, '1627_1998_step', 0
	lcplot, '1900_step', 0
	lcplot, '1810_step', 0  
	lcplot, '1048_step', 0  
	lcplot, '1647_step', 0  
endif


	if keyword_set(lcol) then !p.color=lcol else !p.color=255

	if not keyword_set(noplot) then	oplot, tm,Fm,col=250
	
	
	if not keyword_set(nodata) then begin
		
		if (strcmp(source,'0556',4)) then begin
			
				readcol, 'data/0556', t, F, dF,Fp,dFp, format=('D,D,D,X,X,D,D')
				fac = 1e-13 * 4.0*!dpi*(3.086d21*42.7)^2 
				dF*=fac
				F2=fac*(F-Fp)
				F=fac*F
				oploterror, t, F, dF, psym=1, /nohat, symsize=1
				;oploterror, t, F2, dF, psym=1, /nohat, symsize=1,col=120,errcol=120
				dd=42.7
		endif
		

		if (strcmp(source,'1647',4)) then begin
		;	readcol, 'data/1647', t0, format=('D')
			readcol, 'data/1647', t, F, dF, format=('D,X,X,D,D')
			F*=1d35
			dF*=1d35
			dd=5.0
			Lq=7e32
			Lq*=(dd/5.0)^2
			dF*=(dd/5.0)^2
			F*=(dd/5.0)^2
			oploterror, t, F, dF, psym=1, /nohat, symsize=0.7

	
		endif

		if (strcmp(source,'1627',4) and 0) then begin
			readcol, 'data/1998cooling.dat', t, F, dF, format=('D,D,D')
			F*=1d35
			dF*=1d35
			oploterror, t, F, dF, psym=6, /nohat, symsize=0.5
			readcol, 'data/2008cooling.dat', t1, F1, dF1, format=('D,D,D')
			F1*=1d35
			dF1*=1d35
			oploterror, t1, F1, dF1, psym=1, /nohat
			dd=11.0
		endif

		if (strcmp(source,'1627_1998',9)) then begin
			readcol, 'data/1998cooling.dat', t, F, dF, format=('D,D,D')
			F*=1d35
			dF*=1d35
			oploterror, t, F, dF, psym=6, /nohat, symsize=0.7
			dd=11.0
		endif
		if (strcmp(source,'1627_2008',9)) then begin
			readcol, 'data/2008cooling.dat', t, F, dF, format=('D,D,D')
			F*=1d35
			dF*=1d35
			oploterror, t[0:n_elements(t)-2], F[0:n_elements(t)-2], dF[0:n_elements(t)-2], psym=6, /nohat, symsize=0.7
			dd=11.0
			oploterror, t[n_elements(t)-1], F[n_elements(t)-1], dF[n_elements(t)-1], psym=6, /nohat, symsize=0.7,/lobar			
		endif
		
		if (strcmp(source,'fluxes',6)) then begin
			readcol, 'data/'+source, t0, format=('D')
			readcol, 'data/'+source, t, F, dF,dF2,format=('D,X,D,D,D'),skipline=1		
			t-=t0[0]
			dd=8.0
			if (strcmp(source,'fluxes1822',10)) then dd=1.6
			if (strcmp(source,'fluxes1547',10)) then dd=3.9
			F*=4.0*!dpi*(3.086d21*dd)^2
			dF*=4.0*!dpi*(3.086d21*dd)^2
			oploterror, t, F, dF, psym=1,/nohat, /lobar;, col=250, errcol=250
			oploterror, t, F, dF2, psym=1,/nohat, /hibar;, col=250, errcol=250

				
			if (0)then begin ;(strcmp(source,'fluxes1822',10) and not strcmp(source,'fluxes1822bolo',14)) then begin
				readcol, 'data/'+source, t2, F2, dF2,flag, format=('D,X,D,D,X,I'),skipline=1		
				t2-=t0[0]
				ind = where(flag eq 1)			
				F2*=4.0*!dpi*(3.086d21*dd)^2
				dF2*=4.0*!dpi*(3.086d21*dd)^2
				oploterror, t2[ind], F2[ind], dF2[ind], psym=6,/nohat,symsize=1 ;, col=250, errcol=250
				t = [t,t2[ind]]
				F = [F,F2[ind]]
				ind = sort(t)
				t=t[ind]
				F=F[ind]
			endif	
			
			if (strcmp(source,'fluxes1547',10)and 1) then begin
				readcol, 'data/'+source+'_2', t0, format=('D')
				readcol, 'data/'+source+'_2', t, F, dF,format=('D,X,D,D'),skipline=1		
						t-=t0[0]
								F*=4.0*!dpi*(3.086d21*dd)^2
								dF*=4.0*!dpi*(3.086d21*dd)^2
								oploterror, t, F, dF, psym=5,/nohat
				
				
			endif
			
			
			
			
			
		endif
			
		if (strcmp(source,'0501',4) or strcmp(source,'j1833',5)) then begin
			readcol, 'data/'+source+'lum', t, F, format=('D,D'),skipline=1		
			F*=1d34
			dF=F*0.1
			dd=5
			if (strcmp(source,'j1833',5)) then dd=10.0
			oplot, t, F, psym=1			
		endif
		
		if (strcmp(source,'1810',4)||strcmp(source,'1048',4)) then begin
			readcol, 'data/'+source, t0, format=('D')
			readcol, 'data/'+source, t, F, format=('D,D'),skipline=1
			t-=t0[0]
			dd=3.5
			if (strcmp(source,'1048',4)) then dd=2.7
			F*=4.0*!dpi*(3.086d21*dd)^2
			dF=F*0.1
			oplot, t, F, psym=1				
		endif

		if (strcmp(source,'SGR1900',7)) then begin
			readcol, 'data/SGR1900', t0, format=('D'), numline=1
			readcol, 'data/SGR1900', t, F, format=('D,D'),skipline=1
			F*=1d-11
			t-=t0[0]
			dd=13.5 ;;  for 1900
			F*=4.0*!dpi*(3.086d21*dd)^2
			dF=F*0.1
		;	print, t, F, dF
		;	oplot, t, F, psym=6
			
				readcol, 'data/1900', t, F, format=('D,D')
				F*=1d35
				dF=F*0.1
				dd=13.5
				Lq=1.0d35   ; Hongjun said 0.8 to 1.3 e35
				
				F*=(dd/13.5)^2
				dF*=(dd/13.5)^2
				
				oplot, t, F, psym=6, symsize=0.5
			
						dd=5.5
						Lq=1.0d35   ; Hongjun said 0.8 to 1.3 e35

						F*=(dd/13.5)^2
						dF*=(dd/13.5)^2

						oplot, t, F, psym=6, symsize=0.5
			
						xyouts, 1.0,1d35,textoidl('d=5.5 kpc'),charsize=0.55
						xyouts, 100.0,8d35,textoidl('d=13.5 kpc'),charsize=0.55
			
		endif
		
		if (strcmp(source,'2259',4) ) then begin
			readcol, 'data/'+source, t0, format=('D'), numline=1
			readcol, 'data/'+source, t, F, dF, format=('D,D,D'),skipline=1
			F*=1d-11
			dF*=1d-11
			t-=t0[0]

			;	dd=11.0
			dd=16.0 ;;  for 1900
			if (strcmp(source,'2259',4)) then dd=3.0
			print, 'Assuming distance=',dd

			F*=4.0*!dpi*(3.086d21*dd)^2
			dF*=4.0*!dpi*(3.086d21*dd)^2

			print, t, F, dF

			if (strcmp(source,'2259',4)) then begin
				source='2259_xtepf'
				readcol, 'data/'+source, t0, format=('D'), numline=1
				readcol, 'data/'+source, t1, F1, dF1, format=('D,D,D'),skipline=1
				t1-=t0[0]
				dF1*=F[n_elements(F)-1]/F1[n_elements(F1)-1]
				F1*=F[n_elements(F)-1]/F1[n_elements(F1)-1]
				oploterror, t1, F1, dF1, psym=1,/nohat
			endif
			oploterror, t, F, dF, psym=5,/nohat
			F=F1  
			t=t1
			dF=dF1
		endif
	
	if keyword_set(source) then begin
	; now calculate chi-squared
		Lmodel = interpol(Fm,tm,t)
		chi = total(((Lmodel-F)/dF)^2)
		print, 'chisq=', chi, ' chisq/N=', chi/n_elements(t)

		print, 'distance=', dd

		ind = sort(t)
		L10 = interpol(F[ind],t[ind],10.0)
		print, 'L10=',L10
;		print,t[ind],F[ind]

endif
	
	endif



	if not keyword_set(nolabel) then begin
		
		if not keyword_set(cs) then cs=1.3
	; put a label
	if (strcmp(source,'2259',4)) then xyouts, 100.0, 7d34, '1E 2259+586', charsize=cs
;	if (strcmp(source,'fluxes1822',10)) then xyouts, 100.0, 7d34, 'Swift J1822-1606', charsize=cs
	if (strcmp(source,'fluxes1547',10)) then xyouts, 60.0, 4d35, '1E 1547.0-5408', charsize=cs
	;if (strcmp(source,'1627',4)) then xyouts, 200.0, 4d35, 'SGR 1627-41', charsize=1.5
	;if (strcmp(source,'1627_1998',9)) then xyouts, 200.0, 2d35, 'SGR 1627-41', charsize=cs
	if (strcmp(source,'1627_1998',9)) then xyouts, 200.0, 2d35, '1998 outburst', charsize=cs
	;if (strcmp(source,'1627_1998',9)) then xyouts, 200.0, 1.4d35, '1998 outburst', charsize=cs
	;if (strcmp(source,'1627_2008',9)) then xyouts, 200.0, 6d34, 'SGR 1627-41', charsize=cs
;	if (strcmp(source,'1627_2008',9)) then xyouts, 200.0, 4.5d34, '2008 outburst', charsize=cs
	if (strcmp(source,'1627_2008',9)) then xyouts, 200.0, 6d34, '2008 outburst', charsize=cs
	if (strcmp(source,'0501',4)) then xyouts, 100.0, 4d35, 'SGR 0501+4516', charsize=cs
	if (strcmp(source,'j1833',5)) then xyouts, 100.0, 1d35, 'SGR 1833-0832', charsize=cs
	if (strcmp(source,'1810',4)) then xyouts, 100.0, 2d35, 'XTE J1810-197', charsize=cs
	if (strcmp(source,'1048',4)) then xyouts, 1200.0, 6d34, '1E 1048.1-5937', charsize=cs
	if (strcmp(source,'1647',4)) then xyouts, 20.0, 1.5d35, 'CXOU J164710.2-455216', charsize=cs
	if (strcmp(source,'SGR1900',7)) then xyouts, 100.0, 1d37, 'SGR1900+14', charsize=cs
	
	endif
	

	if keyword_set(nolabel) then begin
	; put a label
	if (strcmp(source,'2259',4)) then xyouts, 0.2, 3d32, '1E 2259+586', charsize=1.0
	if (strcmp(source,'fluxes1822',10)) then xyouts, 0.2, 1.5*3d32, 'Swift J1822-1606', charsize=1.0
	if (strcmp(source,'fluxes1547',10)) then xyouts, 0.2, 1.5^2*3d32, '1E 1547.0-5408', charsize=1.0
	;if (strcmp(source,'1627',4)) then xyouts, 200.0, 435, 'SGR 1627-41', charsize=1.5
	if (strcmp(source,'1627_1998',9)) then xyouts, 0.2, 1.5^3*3d32, 'SGR 1627-41 (1998)', charsize=1.0
	if (strcmp(source,'1627_2008',9)) then xyouts, 0.2, 1.5^4*3d32, 'SGR 1627-41 (2008)', charsize=1.0
	if (strcmp(source,'0501',4)) then xyouts, 0.2, 1.5^5*3d32, 'SGR 0501+4516', charsize=1.0
	if (strcmp(source,'j1833',5)) then xyouts, 0.2, 1.5^6*3d32, 'SGR 1833-0832', charsize=1.0
	if (strcmp(source,'1810',4)) then xyouts, 0.2, 1.5^7*3d32, 'XTE J1810-197', charsize=1.0
	if (strcmp(source,'1048',4)) then xyouts, 0.2, 1.5^8*3d32, '1E 1048.1-5937', charsize=1.0
	if (strcmp(source,'1647',4)) then xyouts, 0.2, 1.5^9*3d32, 'CXOU J164710.2-455216', charsize=1.0
	if (strcmp(source,'SGR1900',7)) then xyouts, 0.2, 1.5^10*3d32, 'SGR1900+14', charsize=1.0
	endif



	if (Lq gt 0.0) then begin

		oplot, [0.1,1d5],[Lq,Lq],linestyle=1

	endif

	if not keyword_set(nodata) then begin
	
		ind=sort(t)
		t=t[ind]
		F=F[ind]
	
		I1=0.0
		I2=0.0
		I3=0.0
		lastt = t[n_elements(t)-1]
		for i=1, n_elements(t[where(t le lastt)])-1 do begin
			I1+=F[i]*(t[i]-t[i-1])
			I2+=F[i-1]*(t[i]-t[i-1])
			I3+=0.5*(F[i-1]+F[i])*(t[i]-t[i-1])
		endfor
		I1*=24.0*3600.0
		I2*=24.0*3600.0
		I3*=24.0*3600.0
		print, 'Total energy in lightcurve=',I1,I2,I3
		t=tm
		F=Fm
		I1=0.0
		I2=0.0
		I3=0.0
		lastt=10000.0
		for i=1, n_elements(t[where(t le lastt)])-1 do begin
			I1+=F[i]*(t[i]-t[i-1])
			I2+=F[i-1]*(t[i]-t[i-1])
			I3+=0.5*(F[i-1]+F[i])*(t[i]-t[i-1])
		endfor
		I1*=24.0*3600.0
		I2*=24.0*3600.0
		I3*=24.0*3600.0
		print, 'Total energy in lightcurve=',I1,I2,I3

endif
	
		if (0) then begin
		; show a scaling
		t=[alog10(10.0),alog10(50.0)]
		L=33.75 + (t-1.0)*0.15
		oplot, 10^t, 10^L, thick=3	
		xyouts, 20.0,4.3e33,textoidl('\propto t^{0.15}')

		t=[alog10(10.0),alog10(50.0)]
		L=34.8 - (t-1.5)*1.32*(15.0/29.0)
		oplot, 10^t, 10^L, linestyle=0, thick=3	
		xyouts, 20.0,1e35,textoidl('\propto t^{-0.7}')


		;t=[alog10(10.0),alog10(50.0)]
		;L=34.8 - (t-1.5)*0.4
		;oplot, 10^t, 10^L, linestyle=0, thick=3	
		;xyouts, 20.0,1e35,textoidl('\propto t^{-0.7}')
	endif

	
;	if (strcmp(source,'fluxes1822',10)) then begin		
;		lcsum, '1822', 2
;	endif

;	lcsum, '1647', 2


		readcol, 'data/RCS_pars_cflux_2.txt', t, F, format=('D,X,X,X,X,D,X,X,X,X,X,X,X,X')
		F=10^F
		F*=(3.086d21*5.0)^2*4.0*!dpi
		t-=t[0]
		;oplot, t, F, psym=4, col=120

		readcol, 'data/RCS_pars_cflux_0.5.txt', t, F, format=('D,X,X,X,X,D,X,X,X,X,X,X,X,X')
		F=10^F
		F*=(3.086d21*5.0)^2*4.0*!dpi
		t-=t[0]
		;oplot, t, F/1.7, psym=4, col=80

		readcol, 'data/Flux_bb_pow_0.5_10.txt', t, F, format=('D,D')
		;F=10^F
		F*=(3.086d21*5.0)^2*4.0*!dpi
		t-=t[0]
		;oplot, t, F/3.0, psym=5, col=250



		if keyword_set(ps) then close_ps	
end








pro surf
	; plots the Teff-Tb (or actually plots flux-Tb) relation
	; and compares to the T_top - Flux from the simulation
	!p.multi=[0,1,1,0,0]

	ytop=12.0   ; log10 of column depth of outer grid point
	;ytop=13.6   ; log10 of column depth of outer grid point

	readcol, 'out/grid', y, T, F, format=('D,D,D')
	F=10^(F-21.0)
	T=10^T
	ind=where(abs(y-ytop) lt 0.01)
	plot, F[ind], T[ind],linestyle=0,/xlog,/ylog,$
			ytitle=textoidl("T_b"), xtitle=textoidl('F_{21}'), charsize=1.4

		readcol, 'grid_sorty', y, T, F, format=('D,D,D')
		F=10^(F-21.0)
		T=10^T
		ind=where(abs(y-ytop) lt 0.01)
		oplot, F[ind], T[ind], linestyle=1

	readcol, 'gon_out/prof', Teff,T, format=('X,X,X,X,F,F')
	F*=1d-21
	print, 'max F21 in cooling curve=',max(F)
	oplot, F,T, thick=3


	g14=2.28
	T=dindgen(100)*0.03 + 6.5
	
	; This is the Gudmundsson Tb-Teff
	F=1d-21*1d24*g14*5.67e-5*(10^T/1.288e8)^(1.0/0.455)
	oplot, F,10^T, col=120
	
	; This is the Potekhin fit
	Faccr=alog10(g14)+2.42*alog10(18.1)+2.42*(T-9.0)+24.0
	T9=10^(T-9.0)
	C=T9-sqrt(7.0*T9*sqrt(g14))/1d3
	Firon=alog10(g14)+alog10((7.0*C)^2.25+(0.333*C)^1.25)+24.0
	eta=g14^2*(10^ytop)*4.0*!dpi*1d12/(1.4*2d33)
	a=alog10(1.2+(5.3d-6/eta)^0.38)+1.667*(T-9.0)
	a=10^a
	F=alog10(10^Faccr + a*10^Firon)-alog10(1+a)
	F=10^(F-21.0)
	F*=5.67e-5
	T=10^T
;	print, T, F
	oplot, F, T, col=250
	
			
end


pro surf2
	; plots the Teff-Tb (or actually plots flux-Tb) relation
	; and compares to the T_top - Flux from the simulation
	!p.multi=[0,1,1,0,0]

	ytop=12.0   ; log10 of column depth of outer grid point
	;ytop=13.6   ; log10 of column depth of outer grid point

	readcol, 'out/grid', y, T, F, format=('D,D,D')
	F=10^F
	Teff = (F/5.67e-5)^0.25
	T=10^T
	ind=where(abs(y-ytop) lt 0.01)

	plot, Teff[ind], T[ind],linestyle=1,/xlog,/ylog,$
			ytitle=textoidl("T_b"), xtitle=textoidl('T_{eff}'), charsize=1.4

		readcol, 'out/grid', y, T, F, format=('D,D,D')
		F=10^F
		Teff = (F/5.67e-5)^0.25
		T=10^T
		ind=where(abs(y-ytop) lt 0.01)
		
	;	oplot, Teff[ind], T[ind]

	readcol, 'gon_out/prof', T,Teff, format=('X,X,X,X,X,F,F')
	oplot, Teff,T, thick=3

	
	g14=2.28
	T=dindgen(100)*0.03 + 6.5
	
	; This is the Gudmundsson Tb-Teff
	F=1d24*g14*5.67e-5*(10^T/1.288e8)^(1.0/0.455)
	Teff = (F/5.67e-5)^0.25
	oplot, Teff,10^T, col=120
	
	; This is the Potekhin fit
	yhe=9.0
	Faccr=alog10(g14)+2.42*alog10(18.1)+2.42*(T-9.0)+24.0
	T9=10^(T-9.0)
	C=T9-sqrt(7.0*T9*sqrt(g14))/1d3
	Firon=alog10(g14)+alog10((7.0*C)^2.25+(0.333*C)^1.25)+24.0
	eta=g14^2*(10^yhe)*4.0*!dpi*1.12d6^2/(1.4*2d33)
	a=alog10(1.2+(5.3d-6/eta)^0.38)+1.667*(T-9.0)
	a=10^a
	F=alog10(10^Faccr + a*10^Firon)-alog10(1+a)
	F=10^F
	F*=5.67e-5
	T=10^T
	Teff = (F/5.67e-5)^0.25
	oplot, Teff, T, col=250
	a=0.0
	F=10^Faccr
	F*=5.67e-5
	Teff = (F/5.67e-5)^0.25
	oplot, Teff, T, col=250

	
;	readcol, 'TeTb-8-9-2.10.data', Teff, Tb, format=('D,D')
;	oplot, Teff, Tb, linestyle=3
end


pro tbteff

	readcol, 'gon_out/TbTeff', T, F, format=('X,X,X,X,D,D')
	radius=11.2
	ZZ=1.32
	F/=5.67d-5
	F=F^0.25
	F/=ZZ
	F*=1.38d-16/1.6d-12
	plot, T,F, /xlog,/ylog
	
	T=[7.0,8.0]
	F=alog10(80.0)+0.5*(T-7.0)
	oplot, 10^T, 10^F, linestyle=1
	
end




pro surf3, ps=ps

	if keyword_set(ps) then open_ps,'surf3.ps'
	
;	surf3plot, 1d13,/nocorr
;	surf3plot, 1d16,/overplot, ls=0
;	surf3plot, 1d15,/overplot, ls=0,/nocorr
;	surf3plot, 1d14,/overplot, ls=3
surf3plot, 1d14, ls=2,/nocorr
surf3plot, 1d15, ls=2,/nocorr,/overplot
	surf3plot, 1.0, ls=2,/overplot, col=80
;	surf3plot, 1.0,/overplot, /usearras,ls=2
;	surf3plot, 1.0,/overplot, /usegpe,ls=1

	TT=7.0+3.0*0.01*dindgen(100)
	LL=1.7d35 * (10^(TT-9.0))^2.1
;	oplot, LL,10^TT, col=250
	TT=7.0+3.0*0.01*dindgen(100)
	LL=0.93d35 * (10^(TT-9.0))^2.1
;	oplot, LL,10^TT, col=250

	oplot, [2d31,3.3d31],[2.7d9,2.7d9], linestyle=2
	xyouts, 3.7d31, 2.6d9, textoidl('Potekhin & Yakovlev 2001'), charsize=1.01
	oplot, [2d31,3.3d31],[4.1d9,4.1d9],linestyle=0
	xyouts, 3.7d31,4d9, textoidl('This work'), charsize=1.01
;;	oplot, [2d31,3d31],[1.65d9,1.65d9], linestyle=1
;	xyouts, 3.5d31, 1.55d9, textoidl('This work B=0 G, y=10^{14} g cm^{-2}'), charsize=1.2
;	oplot, [2d31,3d31],[1.00d9,1.00d9], linestyle=2
;	xyouts, 3.5d31, 0.95d9, textoidl('Arras   B=10^{15}G solid surface'), charsize=1.2

xyouts, 7d31, 1d8, textoidl('B=0'), charsize=1.2, col=80
xyouts, 1d33, 2.8d7, textoidl('B=10^{14} G'), charsize=1.2

xyouts, 3d35, 2d8, textoidl('log_{10} y = 9'), charsize=1.1
xyouts, 4d36, 1.7d9, textoidl('12'), charsize=1.1
xyouts, 5d35, 3d9, textoidl('14'), charsize=1.1, col=80
xyouts, 4d36, 2.7d9, textoidl('14'), charsize=1.1




;	readcol, 'gon_out/prof', Teff,T, format=('X,X,D,X,X,D')
;	F=Teff*4.0*!dpi*(11.2*1d5)^2
;	print, 'max F21 in cooling curve=',max(F)
;	print, T
;	oplot, F,T, thick=3, col=250

	if (1) then begin
	print, 'gon_out/TbTeff'
	readcol, 'gon_out/TbTeff', T, F, format=('X,X,X,X,D,D')
	radius=11.2
		ZZ=1.32
	F *= 4.0*!dpi*1d10*radius^2
	F/=ZZ^2
;	oplot, F, T, col=250

		print, 'out/grid_1e14_nopotek'
		readcol, 'out/grid_1e14_nopotek', y, T, F, format=('D,D,D')
		F=10^F
		T=10^T
		radius=11.2
		ZZ=1.32
		F *= 4.0*!dpi*1d10*radius^2
		F/=ZZ^2

		ytop=14.0
		ind=where(abs(y-ytop) lt 0.01)
		oplot, F[ind], T[ind]
		ytop=12.0
		ind=where(abs(y-ytop) lt 0.01)
		oplot, F[ind], T[ind]
		ytop=9.0
		ind=where(abs(y-ytop) lt 0.01)
		oplot, F[ind], T[ind]

		print, 'out/grid_1e15_nopotek'
		readcol, 'out/grid_1e15_nopotek', y, T, F, format=('D,D,D')
		F=10^F
		T=10^T
		radius=11.2
		ZZ=1.32
		F *= 4.0*!dpi*1d10*radius^2
		F/=ZZ^2

		ytop=14.0
		ind=where(abs(y-ytop) lt 0.01)
		oplot, F[ind], T[ind], col=120, linestyle=0
		ytop=12.0
		ind=where(abs(y-ytop) lt 0.01)
		oplot, F[ind], T[ind], col=120, linestyle=0
		ytop=9.0
		ind=where(abs(y-ytop) lt 0.01)
		oplot, F[ind], T[ind], col=120, linestyle=0


		print, 'out/grid_He4'
		readcol, 'out/grid_He4', y, T, F, format=('D,D,D')
		F=10^F
		T=10^T
		radius=11.2
		ZZ=1.32
		F *= 4.0*!dpi*1d10*radius^2
		F/=ZZ^2
		ytop=14.0
		ind=where(abs(y-ytop) lt 0.01)
		oplot, F[ind], T[ind], linestyle=0, col=80




		B=1d16
		chi = 1.0 + 0.0492 * (1d-12*B)^0.292 / (T[ind]*1d-9)^0.24
		;oplot, F[ind]*chi^4, T[ind], col=120
		B=1d15
		chi = 1.0 + 0.0492 * (1d-12*B)^0.292 / (T[ind]*1d-9)^0.24
		;oplot, F[ind]*chi^4, T[ind], col=120
		B=1d14
		chi = 1.0 + 0.0492 * (1d-12*B)^0.292 / (T[ind]*1d-9)^0.24
		;oplot, F[ind]*chi^4, T[ind], col=120


	endif
	
	if keyword_set(ps) then begin
     	device,/close
     	set_plot,'x'
  	endif
	
end


function PYL, Tc, B, radius=radius, ZZ=ZZ, grav=grav
	if not keyword_set(radius) then radius=12.0
	if not keyword_set(ZZ) then ZZ=1.2364
	if not keyword_set(grav) then grav=1.60355d14

	T9 = Tc*1d-9
	xi = T9 - 0.001*(1d-14*grav)^0.25*sqrt(7.0*T9);
	flux = 5.67d-5 * 1d24 * grav*1d-14 * ((7*xi)^2.25+(0.333*xi)^1.25);
	; now correct for B ... 
	; use the enhancement along the field direction;
;	chi = 1.0 + 0.0492 * (1d-12*B)^0.292 / T9^0.24
;	flux *= chi^4

	; or use eq. (31) or PY2001  which gives F(B)/F(0)
	bet = 0.074*sqrt(1d-12*B)*T9^(-0.45);
	a1=5059.0*T9^0.75/sqrt(1.0 + 20.4*sqrt(T9) + 138.0*T9^1.5 + 1102.0*T9*T9);
	a2=1484.0*T9^0.75/sqrt(1.0 + 90.0*T9^1.5+ 125.0*T9*T9);
	a3=5530.0*T9^0.75/sqrt(1.0 + 8.16*sqrt(T9) + 107.8*T9^1.5+ 560.0*T9*T9);
	fac = (1.0 + a1*bet^2 + a2*bet^3 + 0.007*a3*bet^4)/(1.0+a3*bet^2);
 	flux *= fac


	flux *= 4.0*!dpi*1d10*radius^2
	flux/=ZZ^2

	return, flux
end

pro surf3plot, B, overplot=overplot, ls=ls, usearras=usearras, nocorr=nocorr, usegpe=usegpe, col=col

	

;	B=1d15
	Tc=7.0 + dindgen(60)*0.025*2.5
	Tc=10^Tc
	grav = 2.28d14
	radius=11.2
	ZZ=1.32
;; R=12km, M=1.4M_sun
;	grav = 1.6d14
;	radius=12.0
;	ZZ=1.24
		
	; Potekhin & Yakovlev 2001 eq.(27)
	T9 = Tc*1d-9
	xi = T9 - 0.001*(1d-14*grav)^0.25*sqrt(7.0*T9);
	flux = 5.67d-5 * 1d24 * grav*1d-14 * ((7*xi)^2.25+(0.333*xi)^1.25);
	; now correct for B ... 
	if not keyword_set(nocorr) then begin
	; or use eq. (31) or PY2001  which gives F(B)/F(0)
	bet = 0.074*sqrt(1d-12*B)*T9^(-0.45);
	a1=5059.0*T9^0.75/sqrt(1.0 + 20.4*sqrt(T9) + 138.0*T9^1.5 + 1102.0*T9*T9);
	a2=1484.0*T9^0.75/sqrt(1.0 + 90.0*T9^1.5+ 125.0*T9*T9);
	a3=5530.0*T9^0.75/sqrt(1.0 + 8.16*sqrt(T9) + 107.8*T9^1.5+ 560.0*T9*T9);
	fac = (1.0 + a1*bet^2 + a2*bet^3 + 0.007*a3*bet^4)/(1.0+a3*bet^2);
 	flux *= fac
	endif else begin
	
	; use the enhancement along the field direction;
	;double chi;
	;chi = 1.0 + 0.0492*pow(EOS.B*1e-12,0.292)/pow(T9,0.24);
	;flux *= pow(chi,4.0);
	
	chi = 1.0 + 0.0492 * (1d-12*B)^0.292 / T9^0.24
	flux *= chi^4
	
	endelse

	flux *= 4.0*!dpi*1d10*radius^2
	flux/=ZZ^2

	print, interpol(flux,Tc,[2.5d8]), 'for B=', B
	print, interpol(flux,Tc,[1.5d7]), 'for B=', B



	if keyword_set(usearras) then begin
	
		Ts=4.0d6*(grav/1d14)^0.25*T9^0.5
		flux = 5.67d-5*Ts^4
		flux *= 4.0*!dpi*1d10*radius^2
		flux/=ZZ^2
		
	
	
	endif
	
	
	if keyword_set(usegpe) then begin
	
			Ts=1.0d6*(grav/1d14)^0.25*T9^0.55/(0.1288)^0.55
			flux = 5.67d-5*Ts^4
			flux *= 4.0*!dpi*1d10*radius^2
			flux/=ZZ^2
		
			print, T9, Ts, flux
		
		
	endif
	


	if (0) then begin

	if keyword_set(overplot) then begin
		oplot, flux,Tc/1d7, linestyle=ls
	endif else begin
		plot, flux, Tc/1d7, /xlog, ytitle=textoidl('T_c (10^7 K)'),$
	 		xtitle=textoidl('L_\infty (erg s^{-1})'), charsize=1.5, xrange=[1d31,3d34],$
			yrange=[1,15], ystyle=1
	endelse

endif else begin

	if keyword_set(overplot) then begin
		if keyword_set(col) then begin
			oplot, flux,Tc, linestyle=ls, col=col
		endif else begin
			oplot, flux,Tc, linestyle=ls, col=col
		endelse
	endif else begin
		plot, flux, Tc, /xlog,/ylog, ytitle=textoidl('T_c (K)'),yrange=[1d7,8e9],ystyle=1,$
	 		xtitle=textoidl('L_\infty (erg s^{-1})'), charsize=1.2, xrange=[1d31,1d37],xstyle=1,linestyle=ls
	endelse

endelse
	
end

pro init_plotone,name, ls=ls,lcol=lcol,ener=ener

	fname = 'gon_out/initial_condition'
	if (name ne '') then fname+='_'+name 

	readcol, fname, zone, y, T, rho, CV, K, tt,E,TC,$
						format=('I,D,D,D,D,D,X,X,D,D,X,X,D')
						
	;if keyword_set(ener) then begin
	;	T*=(1d-10*rho)^0.6/ener^0.666
	;	print,T[where(rho lt 1d11)]
	;endif					
						
	if keyword_set(lcol) then begin
		oplot,y/2.28d14,T,linestyle=ls,col=lcol
	endif else begin
		oplot,y/2.28d14,T,linestyle=ls
	endelse


end


pro initialT, ps=ps, ener=ener
	
		if keyword_set(ps) then open_ps,'initialT.ps'
	
		!p.multi=[0,1,1,0,0]
	
		plot, [6d8,3d18], [5d7,1e9], /xlog,/ylog, /nodata, xtitle=textoidl('Column depth (g cm^{-2})'),$
				ytitle=textoidl('T (K)'), ystyle=1, xstyle=1,charsize=1.3

			if (0) then begin
			init_plotone,'B1e14E100.0',ener=100.0
			init_plotone,'B1e14E30.0',ener=30.0
			init_plotone,'B1e14E10.0',ener=10.0
			init_plotone,'B1e14E3.0',ener=3.0
			init_plotone,'B1e14E1.0',ener=1.0
			init_plotone,'B1e14E0.3',ener=0.3

			init_plotone,'B1e15E0.3',ls=2;,lcol=250
			init_plotone,'B1e15E1.0',ls=2;,lcol=250
			init_plotone,'B1e15E3.0',ls=2;,lcol=250
			init_plotone,'B1e15E10.0',ls=2;,lcol=250
			init_plotone,'B1e15E30.0',ls=2;,lcol=250
			init_plotone,'B1e15E100.0',ls=2;,lcol=80
		endif
		if (0) then begin
			init_plotone,'B1e14E100.0_1e9',ener=100.0
			init_plotone,'B1e14E30.0_1e9',ener=30.0
			init_plotone,'B1e14E10.0_1e9',ener=10.0
			init_plotone,'B1e14E3.0_1e9',ener=3.0
			init_plotone,'B1e14E1.0_1e9',ener=1.0
			init_plotone,'B1e14E0.3_1e9',ener=0.3

			init_plotone,'B1e15E0.3_1e9',ls=2;,lcol=250
			init_plotone,'B1e15E1.0_1e9',ls=2;,lcol=250
			init_plotone,'B1e15E3.0_1e9',ls=2;,lcol=250
			init_plotone,'B1e15E10.0_1e9',ls=2;,lcol=250
			init_plotone,'B1e15E30.0_1e9',ls=2;,lcol=250
			init_plotone,'B1e15E100.0_1e9',ls=2;,lcol=80
		endif else begin
		;	init_plotone,'B1e14E3.0_1e9',ener=3.0,ls=3
		;	init_plotone, 'B1e14E3.0_1e9_constantT4.2e8',ls=0
		;	init_plotone, 'B1e14E3.0_1e9_E10.8constantpergram',ls=2
		endelse
		
		init_plotone,'',ls=0,lcol=250
		init_plotone,'1731_1',ls=0
		init_plotone,'1731_2',ls=0
							
							
		readcol, 'ign.CC', T, y,format=('D,X,D,X')
		oplot, y, T, linestyle=2
							
		; plot TB(rho)
		B=1d15
		rho = dindgen(10)*0.1 *3 + 8.5
		rho = 10^rho
		x = 10.0*(1d-9*rho*0.5)^0.333
		T = 1.34e9*1d-13*B/x
		;oplot, rho, T, linestyle=4
		
		
		;T = alog10(3.5d8) - 0.6*(rho-10.0)
		;if keyword_set(ener) then T+=2.0*alog10(ener)/3.0
		;T+=2.0*alog10(10.8*10^(rho-11.0))/3.0
		;oplot, 10^rho, 10^T, linestyle=3
							
							
		; read Gamma/T for the initial model, used to plot melting curve
		readcol, 'gon_out/grid_profile', ym, GammaT, format=('X,X,D,X,X,X,X,X,X,D')
		Tm = 1e8 * GammaT/175.0
		print, ym,Tm
		;oplot, ym, Tm, linestyle=1
		;xyouts, 1.5e11,8e8,textoidl('\Gamma=175')	,charsize=0.9, orientation=10		
		;xyouts, 6d10,6d9, textoidl('t=1 hour'),charsize=1.1
								
		if keyword_set(ps) then close_ps
end





pro expansion,ps=ps
	
		if keyword_set(ps) then begin
	   	set_plot, 'ps'
	   	device,filename='expansion.ps',/color
	  	endif
	 
	!p.charsize=3.0
	
	readcol, 'gon_out/initial_condition', zone, P, T, rho, CV, K, tt,E,TC,Kp,tanal,cve,cvion,chirho,chiT,H,$
						format=('I,D,D,D,D,D,X,X,D,D,X,X,D,X,D,D,D,D,D,D,D')

	!p.multi=[0,1,3,0,0]

	plot, P, 2.303*H*chiT/chirho,/xlog,/ylog, ytitle=textoidl('d\delta z/d\log_{10} P (cm)'),xtitle=textoidl('Pressure (erg cm^{-3})')


	n = n_elements(P)
	z = dblarr(n)
	z[n-1] = 0.0
	for i=2,n_elements(P) do begin
		z[n-i]=z[n-(i-1)] + alog(P[n-(i-1)]/P[n-i])*H[n-i]*chiT[n-i]/chirho[n-i]
	endfor
	plot, P, z, /xlog, /ylog, yrange=[0.1,1e4]

	II = dblarr(n)
	II[n-1] = 0.0
	for i=2,n_elements(P) do begin
		II[n-i]=II[n-(i-1)] + alog10(P[n-(i-1)]/P[n-i])*rho[n-i]*H[n-i]*2*z[n-i]/1e6
	endfor
	II*=8.0*!dpi*1d24/3.0
	II/=2d33*1d12
	plot, P, II, /xlog, /ylog, yrange=[1d-11,1d-7]
	
	
	if keyword_set(ps) then begin
     	device,/close
     	set_plot,'x'
  	endif
end
	


pro initial,ps=ps
	; plot some properties of the profile at t=0
	; shows T(y), CV, K
	
		if keyword_set(ps) then begin
	   	set_plot, 'ps'
	   	device,filename='initial.ps',/color
	  	endif
	 
	!p.charsize=2.0
	
	readcol, 'gon_out/initial_condition', zone, y, T, rho, CV, K, tt,E,TC,Kp,tanal,cve,cvion,$
						format=('I,D,D,D,D,D,X,X,D,D,X,X,D,X,D,D,D,D')
	readcol, 'gon_out/initial_condition', zone2, y2, T2, rho2, CV2, K2, tt2,$
						format=('I,D,D,D,D,D,X,X,D')

			;	y2*=2.28d14


			;rho=y/2.28d14
			;rho2=y/2.28d14

	!p.multi=[0,2,2,0,0]
	plot, rho, T,/xlog,/ylog, ytitle=textoidl('T (K)'),xtitle=textoidl('Column depth (g cm^{-2})'),xrange=[1d7,1d14]
	oplot,rho2,T2,linestyle=1
	oplot, rho,TC+1e5, linestyle=2

	plot, rho, CV,/xlog,/ylog, ytitle=textoidl('C_V'),xtitle=textoidl('Column depth (g cm^{-2})'), yrange=[1d3,5d7],xrange=[1d11,1d14]
;	oplot,rho2,CV2,linestyle=1
	oplot,rho,cve,linestyle=2
;	oplot,rho,cve*(rho/1e9)^0.3333*1e9/T,linestyle=2
	oplot,rho,cvion,linestyle=1

	plot, rho, K/(1e-9*T*(rho*1d-9)^0.333), /xlog,/ylog, ytitle=textoidl('K'),xtitle=textoidl('Column depth (g cm^{-2})'),$
				yrange=[1d13,1d20]
	oplot,rho2,K2,linestyle=1
	oplot,rho,Kp,linestyle=2

	plot, rho, tt, /xlog,/ylog,xtitle=textoidl('Column depth (g cm^{-2})'),$
				ytitle=textoidl('Thermal time (d)'), yrange=[1.0,1e4]
	oplot,rho2,tt2,linestyle=1	
	
	oplot, rho, tanal*4.0, col=250
	oplot, rho, 2.0*(rho/1e9)^1.0, col=80
	oplot, rho, (y/(1.6d14*rho*1000.0))^2, linestyle=1
	
	oplot, rho,  0.05*(y/(1.6d14*rho*1000.0))^2 * rho*1e-9*1d-7*CV/(1d-18*K), col=120
;	plot, rho, E, /xlog,/ylog, yrange=[1d42,1d45]
	
	if keyword_set(ps) then begin
     	device,/close
     	set_plot,'x'
  	endif
end
	
	
pro crust_model
	; plot some properties of the profile at t=0
	; shows T(y), CV, K 
	readcol, 'data/crust_model', P, rho, mue, Qheat, Z, A, Qimp, Yn, $
						format=('D,D,D,D,D,D,D,D')

	!p.multi=[0,2,2,0,0]
	plot, rho, Qimp+1e-3,/xlog,/ylog, ytitle=textoidl('Q_{imp}'),xtitle=textoidl('\rho (g cm^{-3})')

	plot, rho, A,/xlog, ytitle=textoidl('A, Z'),xtitle=textoidl('\rho (g cm^{-3})')
	oplot, rho, Z
	
	plot, rho, Qheat, /xlog, ytitle=textoidl('Q_{heat} (MeV)'),xtitle=textoidl('\rho (g cm^{-3})')

	plot, rho, Yn, /xlog, ytitle=textoidl('Yn'),xtitle=textoidl('\rho (g cm^{-3})')

	
end	
	
	
pro prof2, delay=delay, png=png, source=source
	; We need   window, retain=2   to avoid an X-windows error

	; Animates the crust cooling
	; delay - specifies a delay to slow down the animation
	; png  - if set output png files to make a movie (see makefile)
	; source - a string '1659' '1731' 'XTEJ' determines which data points to plot
	!p.multi=[0,1,2,0,0]

	; read Gamma/T for the initial model, used to plot melting curve
	readcol, 'gon_out/grid_profile', ym, Tmelt,GammaT, format=('X,D,X,X,X,X,X,X,D,D')
	Tm = 1e8 * GammaT/175.0
	; read initial profile
	readcol, 'gon_out/initial_condition', y0, T0, format=('X,D,D')
	; read lightcurve
	readcol, 'gon_out/prof', tt,Teff, format=('F,X,X,X,F')
	tt/=(24*3600.0)
	; the output is already redshifted
	FF = 1.38d-16*Teff/(1.6d-12)
	;FF = 1.38d-16*(FF/5.67d-5)^0.25/(1.6d-12*1.31)
	;tt*=1.31
	; read lightcurve for 1659
	readcol, 'gon_out/prof_1659_1_y9', ttx,FFx, format=('D,X,X,X,D')
	ttx/=(24*3600.0)
	FFx = 1.38d-16*FFx/1.6d-12  ;(FFx/5.67d-5)^0.25/(1.6d-12)
	;ttx/=1.31
	; read observations
	if not keyword_set(source) then source='tc'
	readcol, 'data/'+source, t00, n, format=('D,I'), numline=1
	readcol, 'data/'+source, t, temp, tempe, skipline=1, format=('D,X,X,D,D')
	tobs=t-t00[0]+(5.0-1.6)
	Fobs=temp
	Fobse=tempe

	; now look at the detailed profiles as a function of time
	openr, lun, 'gon_out/out', /get_lun

	; first find out how many grid points
	ngrid=0
	readf, lun, ngrid, format='(I0)'
	print, 'Number of grid points=',ngrid

	multiplot,/reset
	!p.multi=[0,1,3,0,0]

	; animate
	count=0
	oldtime = 0 
	starttime = -999.0
	while (not eof(lun)) do begin

		; get the time
		readf, lun, time
		if not keyword_set(png) then begin
			print, 'time=', time
		endif
		if (starttime eq -999.0) then starttime=time

		; read in next batch of data
		data=dblarr(15,ngrid)
		readf, lun, data
		y=data(0,*)
		T=data(1,*)
		F=data(2,*)
		beta=data(12,*)
		gamma=data(13,*)

		if (1) then begin ;time gt 0.0) then begin

		; plot upper panel
		erase
		if (1) then begin 
		plot, y, T, /xlog, /ylog,charsize=1.5, ytitle=textoidl('T (K)'),$
				xtitle=textoidl('P (g cm^{-2})'), yrange=[1e7,1e9],ystyle=1, $
				xrange=[1d23,8d32], xstyle=1,/nodata
		ind = where(gamma le 175.0,cnt)
		if (cnt gt 0) then begin
			oplot, y[ind], T[ind], thick=3, linestyle=0, col=80
		endif
		oplot, y[where(gamma ge 175.0)], T[where(gamma ge 175.0)], thick=3, linestyle=0, col=120
		oplot, y0, T0, linestyle=1, thick=2
		oplot, ym, Tm, linestyle=2
		;		oplot, y,T*gamma/175.0, linestyle=2, col=250
		;oplot, ym,Tmelt, linestyle=2   ; plot Tmelt which is at the grid centers
		xyouts, 1d30, 5d9, textoidl('t='+string(time/(3600.0*24.0))+' days')
		endif 
		if (1) then begin
			plot, y, abs(F), /xlog, /ylog, charsize=1.5, ytitle=textoidl('Flux (erg cm^{-2} s^{-1})'),$
					xtitle=textoidl('P (g cm^{-2})'),/nodata, yrange=[1d19,3d22], ystyle=1,$
					xrange=[1d23,8d32], xstyle=1
			ind = where(gamma le 175.0 and F gt 0.0,cnt)
			if (cnt gt 0) then begin
				plotsym, 0, /fill
				oplot, y[ind], F[ind], thick=3, psym=8, col=80
			endif
			ind = where(gamma le 175.0 and F lt 0.0,cnt)
			if (cnt gt 0) then begin
				plotsym, 0
				oplot, y[ind], -F[ind], thick=3, linestyle=0, col=80, psym=8
			endif
			ind = where(gamma gt 175.0 and F gt 0.0,cnt)
			if (cnt gt 0) then begin
				plotsym, 0,/fill
				oplot, y[ind], F[ind], thick=3, psym=8, col=120
			endif
			ind = where(gamma gt 175.0 and F lt 0.0,cnt)
			if (cnt gt 0) then begin
				plotsym, 0
				oplot, y[ind], -F[ind], thick=3, psym=8, col=120
			endif



			oplot, ym, Tm, linestyle=2
		endif

		if (1) then begin
		; plot lower panel
		tt2=tt[where(tt*24*3600.0 le time,ntt)]
		ff2=ff[where(tt*24*3600.0 le time)]
		ploterror, tobs, Fobs, Fobse,/xlog, xtitle=textoidl('Time (d)'), $
			ytitle=textoidl('T_{eff} (eV)'), charsize=1.5, $
			xrange=[1.0,7d3],xstyle=1,psym=2,yrange=[min(Fobs)-30.0,max(Fobs)+70.0], ystyle=1
		if (ntt gt 1) then begin
			oplot,tt2,FF2,linestyle=0
		endif
		oplot, ttx, FFx, linestyle=1
		endif


		if keyword_set(delay) then begin
			for i=1L,delay do begin
			endfor
		endif

		if keyword_set(png) then begin
			if (abs(alog10(abs(time))-alog10(abs(oldtime))) gt 0.01 or oldtime eq 0) then begin
				count++	
				image = TVRD(0,0,!D.X_Size,!D.Y_Size,True=1, Order=order)
				filename=string(format='("png/",I03,".png")',count)
				Write_PNG,filename,image
				print, 'Writing png for t=',time, alog10(time), alog10(oldtime), count
				oldtime=time
			endif
		endif

	endif

	endwhile
	
	free_lun,lun
	!p.multi=[0,1,1,0,0]
end

	
	
	pro profl, delay=delay, png=png, source=source
		; We need   window, retain=2   to avoid an X-windows error

		; Animates the crust cooling
		; delay - specifies a delay to slow down the animation
		; png  - if set output png files to make a movie (see makefile)
		; source - a string '1659' '1731' 'XTEJ' determines which data points to plot
		!p.charsize=2


openw, lun2, 'out/kippen.dat', /get_lun


;	readcol, 'out/fcontour_15.dat', Pneut, Tneut, format=('X,D,X,X,D')
	readcol, 'out/fcontour.dat', Pneut, Tneut, format=('X,D,X,X,D')
	Pneut=10.0^Pneut * 2.28d14
	Tneut=10.0^Tneut
	readcol, 'out/fcontour_14.dat', Pneut2, Tneut2, format=('X,D,X,X,D')
	Pneut2=10.0^Pneut2 * 2.28d14
	Tneut2=10.0^Tneut2


;	name = 'out/grid_1e15_potek'
;name = 'out/grid_3e15_potek'
;name = 'out/grid_1e14_potek'
	name = 'out/grid_He9'
;	name = 'out/grid_He4'
;	name = 'out/grid'
	print, 'Reading envelope models from '+name+'... (this takes a while)'
	readcol, name, ysteady, Tsteady, Fsteady,rhosteady, format=('D,D,D,D')

		; read Gamma/T for the initial model, used to plot melting curve
		readcol, 'gon_out/grid_profile', ym, rhom, GammaT, format=('X,D,D,X,X,X,X,X,X,D')
		Tm = 1e8 * GammaT/175.0
		print, ym,Tm
		; read initial profile
		readcol, 'gon_out/initial_condition', y0, T0, format=('X,D,D')
		; read lightcurve
		
		readcol, 'gon_out/prof', tt, F2, FF,F3,format=('F,F,F,X,F')	
		tt/=(24.0*3600.0)
		FF*=4.0*!dpi*1.12d6^2
	
		if (0) then begin
		; read lightcurve for 1659
		readcol, 'gon_out/prof_1659', ttx,FFx, format=('F,X,X,X,F')
		ttx/=(24*3600.0)
		FFx = 1.38d-16*(FFx/5.67d-5)^0.25/(1.6d-12)
		;ttx*=1.31
		; read observations
		if not keyword_set(source) then source='1659'
		readcol, 'data/'+source, t00, n, format=('D,I'), numline=1
		readcol, 'data/'+source, t, temp, tempe, skipline=1, format=('D,X,X,D,D')
		tobs=t-t00[0]
		Fobs=temp
		Fobse=tempe
	endif

	if (strcmp(source,'fluxes1822',10)) then begin
	;	if (strcmp(source,'fluxes',6)) then begin
			readcol, 'data/'+source, t0, format=('D')
			readcol, 'data/'+source, tobs, Fobs, dFobs,format=('D,X,D,D'),skipline=1		
			tobs-=t0[0]
		;	dd=8.0
		;	if (strcmp(source,'fluxes1822',10)) then 
		dd=1.6
		;	if (strcmp(source,'fluxes1547',10)) then dd=3.9
			Fobs*=4.0*!dpi*(3.086d21*dd)^2
			dFobs*=4.0*!dpi*(3.086d21*dd)^2
		;	oploterror, t, F, dF, psym=1,/nohat
				
		;	if (strcmp(source,'fluxes1822',10)) then begin
				readcol, 'data/'+source, tobs2, Fobs2, dFobs2,flag, format=('D,X,D,D,X,I'),skipline=1		
				tobs2-=t0[0]
				ind = where(flag eq 1)			
				Fobs2*=4.0*!dpi*(3.086d21*dd)^2
				dFobs2*=4.0*!dpi*(3.086d21*dd)^2
				Fobs2=Fobs2[ind]
				tobs2=tobs2[ind]
				dFobs2=dFobs2[ind]
			;	oploterror, t2[ind], F2[ind], dF2[ind], psym=6,/nohat,symsize=1
			;endif	
		;endif
endif
if (strcmp(source,'1627_1998',9)) then begin
	readcol, 'data/1998cooling.dat', tobs, Fobs, dFobs, format=('D,D,D')
	Fobs*=1d35
	dFobs*=1d35
	oploterror, tobs, Fobs, dFobs, psym=6, /nohat, symsize=0.7
	dd=11.0
endif
if (strcmp(source,'1627_2008',9)) then begin
	readcol, 'data/2008cooling.dat', tobs, Fobs, dFobs, format=('D,D,D')
	Fobs*=1d35
	dFobs*=1d35
	oploterror, tobs, Fobs, dFobs, psym=6, /nohat, symsize=0.7
	dd=11.0
endif

if (strcmp(source,'1647',4)) then begin
	readcol, 'data/1647', t0, format=('D')
	readcol, 'data/1647', tobs, Fobs, format=('D,D')
	Fobs*=1d35
	dFobs=Fobs*0.1
	dd=5.0
	Lq=4.4e32
;	oplot, t, F, psym=6, symsize=0.5
endif
	
	if (strcmp(source,'0501',4) or strcmp(source,'j1833',5)) then begin
		readcol, 'data/'+source+'lum', tobs, Fobs, format=('D,D'),skipline=1		
		Fobs*=1d34
		dFobs=Fobs*0.1
		dd=5
;		oplot, t, F, psym=2				
	endif


	if (strcmp(source,'1810',4)||strcmp(source,'1048',4)) then begin
		readcol, 'data/'+source, t0, format=('D')
		readcol, 'data/'+source, tobs, Fobs, format=('D,D'),skipline=1
		tobs-=t0[0]
		dd=3.5
		if (strcmp(source,'1048',4)) then dd=2.7
		Fobs*=4.0*!dpi*(3.086d21*dd)^2
		dFobs=Fobs*0.1
;		oplot, t, F, psym=2				
	endif



		; now look at the detailed profiles as a function of time
		openr, lun, 'gon_out/out', /get_lun

		; first find out how many grid points
		ngrid=0
		readf, lun, ngrid, grav, format='(I0,F)'
		print, 'Number of grid points=',ngrid
		print, 'Gravity = ', grav

		multiplot,/reset
		!p.multi=[0,1,1,0,0]

		; animate
		count=0
		oldtime = 0 
		while (not eof(lun)) do begin

			; get the time
			readf, lun, time

			; read in next batch of data
			data=dblarr(12,ngrid)
			readf, lun, data
			y=data(0,*)
			T=data(1,*)
			F=data(2,*)
			rho = data(5,*)
			beta=data(10,*)
			epsnu=data(11,*)

			if (time gt 0.0) then begin

			if (time le 1d5 or ((alog10(abs(time))-alog10(abs(oldtime))) gt 0.02)) then begin
					F1 = 0.1*fix((alog10(F[0])-alog10(grav/2.28d14))/0.1)
					F2 = F1+0.1
					Finterp=((alog10(F[0])-alog10(grav/2.28d14))-F1)/(F2-F1)
					print, 'time=', time, ' log F =', alog10(F[0]), ' rounded F=', F1,F2 


			; plot upper panel
			erase
			plot, y, T, /xlog, /ylog, ytitle=textoidl('T (K)'),$
					xtitle=textoidl('P (erg cm^{-3})'), yrange=[1e6,1e9],ystyle=1, $
					xrange=[1d19,7d32], xstyle=1
;				plot, rho, T, /xlog, /ylog, ytitle=textoidl('T (K)'),$
;						xtitle=textoidl('\rho (g cm^{-3})'), yrange=[1e7,1e10],ystyle=1, $
;						xrange=[5d5,1d14], xstyle=1
			oplot, y0, T0, linestyle=1
			oplot, ym, Tm, linestyle=2
			oplot, Pneut2, Tneut2, linestyle=0,col=200
			;oplot, Pneut, Tneut, linestyle=0,col=80

			Fmin=20.0
			deltaF=7.0
			oplot, y, 10^(7.0 + 3.0 * (alog10(F)-Fmin)/deltaF), col=250, linestyle=3
			TF0 = 7.1 + 3.0 * (alog10(F[0])-Fmin)/deltaF
			if (TF0 gt 7.0 and TF0 lt 10.0) then begin
				xyouts, 5d21, 10^(TF0-0.1), textoidl('log_{10} F='+string(alog10(F[0]),format='(D5.1)')), charsize=1.1, col=250
			endif

			oplot, y, 10^(7.0 + 3.0 * (alog10(epsnu*y/grav)-Fmin)/deltaF), col=80, linestyle=4

			ind = where(abs(Fsteady-F1) lt 0.05 and 2.28e14*10.0^ysteady le y[0],nind)
			ind2 = where(abs(Fsteady-F2) lt 0.05 and 2.28e14*10.0^ysteady le y[0],nind2)
			if (nind gt 0 and nind2 gt 0) then begin
				oplot, 2.28e14*10.0^ysteady[ind], 10.0^(Tsteady[ind]+(Tsteady[ind2]-Tsteady[ind])*Finterp), col=120
			endif
			ind = where(abs(Fsteady-F1) lt 0.05 and 2.28e14*10.0^ysteady gt y[0],nind)
			ind2 = where(abs(Fsteady-F2) lt 0.05 and 2.28e14*10.0^ysteady gt y[0],nind2)
			if  (nind gt 0 and nind2 gt 0) then begin
				oplot, 2.28e14*10.0^ysteady[ind], 10.0^(Tsteady[ind]+(Tsteady[ind2]-Tsteady[ind])*Finterp), col=120, linestyle=1
			endif

			if (0) then begin  ; plot the steady-state envelope for the bracketing fluxes
			ind = where(abs(Fsteady-F2) lt 0.05 and grav*10.0^ysteady le y[0],nind)
			if nind gt 0 then begin
				oplot, grav*10.0^ysteady[ind], 10.0^Tsteady[ind], col=180
			endif
			ind = where(abs(Fsteady-F2) lt 0.05 and grav*10.0^ysteady gt y[0],nind)
			if nind gt 0 then begin
				oplot, grav*10.0^ysteady[ind], 10.0^Tsteady[ind], col=180, linestyle=1
			endif
			ind = where(abs(Fsteady-F1) lt 0.05 and grav*10.0^ysteady le y[0],nind)
			if nind gt 0 then begin
				oplot, grav*10.0^ysteady[ind], 10.0^Tsteady[ind], col=180
			endif
			ind = where(abs(Fsteady-F1) lt 0.05 and grav*10.0^ysteady gt y[0],nind)
			if nind gt 0 then begin
				oplot, grav*10.0^ysteady[ind], 10.0^Tsteady[ind], col=180, linestyle=1
			endif
			endif
				
			if (time lt 3.15d7) then begin
				if (time lt 3600.0*24.0*1.0) then begin
					caption='t='+string(time/(3600.0),format='(D5.1)')+' hours'
				endif else begin
					if (time lt 3600.0*24.0*3.0) then begin
						caption='t='+string(time/(3600.0*24.0),format='(D6.2)')+' days'
					endif else begin
						caption='t='+string(time/(3600.0*24.0),format='(D5.1)')+' days'					
					endelse
				endelse
			endif else begin
				caption='t='+string(time/3.15d7,format='(D5.1)')+' yrs'
			endelse
			xyouts, 1d20,4.5d8, textoidl(caption)

			if (0) then begin
			; plot lower panel
			tt2=tt[where(tt*24*3600.0 le time,ntt)]
			ff2=ff[where(tt*24*3600.0 le time)]
			plot, tt2, ff2,/xlog, xtitle=textoidl('Time (d)'), $
				ytitle=textoidl('L_\infty (erg s^{-1})'), linestyle=0,$
				xrange=[0.1,5000.0],xstyle=1,yrange=[5d32,2d35], ystyle=1,/ylog
			;if (ntt gt 1) then begin
			oplot,tt,FF,linestyle=1
			;endif
			;	oploterror, t, F, dF, psym=1,/nohat
						
				oploterror, tobs, Fobs, dFobs, psym=1,/nohat
				
				if (strcmp(source,'fluxes1822',10)) then begin
					oploterror, tobs2, Fobs2, dFobs2, psym=6,/nohat,symsize=1
				endif
			
			;oplot, ttx, FFx, linestyle=1

		endif

			if keyword_set(delay) then begin
			;	for i=1L,delay do begin
			;	endfor
			wait,delay
			endif

			if keyword_set(png) then begin
				if ((time le 1d10) or ((alog10(abs(time))-alog10(oldtime)) gt 0.02)) then begin
					count++	
					image = TVRD(0,0,!D.X_Size,!D.Y_Size,True=1, Order=order)
					filename=string(format='("png/",I03,".png")',count)
					Write_PNG,filename,image
					print, 'Writing png for t=',time, alog10(abs(time)), alog10(oldtime), count
					oldtime=abs(time)
				
				
				if (0) then begin
					count=0
					ind = where(abs(Fsteady-F1) lt 0.01 and 2.28e14*10.0^ysteady le y[0],nind)
					ind2 = where(abs(Fsteady-F2) lt 0.01 and 2.28e14*10.0^ysteady le y[0],nind)
					if nind gt 0 then begin
						for jj=0,(nind-1)/10 do begin
							ii=jj*10
							printf,lun2,time,10.0^rhosteady[ind[ii]], 10.0^(Tsteady[ind[ii]]+(Tsteady[ind2[ii]]-Tsteady[ind[ii]])*Finterp)
							count++
						endfor
					endif
				
						
					for ii=0,n_elements(rho)-1 do begin
						printf,lun2,time, rho[ii],T[ii]
						count++
					endfor

					print, nind,count
				endif
				
				
				
				endif
			endif
			oldtime = abs(time)
			endif
			
			
		endif
		endwhile

		free_lun,lun2
		free_lun,lun

	end


pro initialT10,ps=ps
	if keyword_set(ps) then open_ps,'initialT10.ps'

;	T10,source='B1e15E100.0_1e9',ls=2
;	T10,source='B1e15E30.0_1e9',/overplot,ls=2
;	T10,source='B1e15E10.0_1e9',/overplot,ls=2
;	T10,source='B1e15E3.0_1e9',/overplot,ls=2
;	T10,source='B1e15E1.0_1e9',/overplot,ls=2
;	T10,source='B1e15E0.3_1e9',/overplot,ls=2
	
;T10,source='slope3',ls=0, nread=14, tt=0.01
;T10,source='slope4',ls=0, nread=14, tt=0.01, /overplot
;T10,source='slope5',ls=0, nread=14, tt=0.01, /overplot
;T10,source='',ls=1, nread=15, tt=0.01
;T10,source='',ls=0, nread=15, tt=1.0, /overplot
;T10,source='',ls=0, nread=15, tt=3.0, /overplot
;T10,source='',ls=0, nread=15, tt=10.0, /overplot
;T10,source='',ls=0, nread=14, tt=30.0, /overplot
;T10,source='',ls=0, nread=14, tt=100.0, /overplot

	
	if (1) then begin
	ttshow=0.000001
T10, source='B1e14E10_1e9_mu1_env',ls=1, tt=ttshow, nread=15
T10, source='B1e14E1_1e9_mu1_env',ls=1, tt=ttshow, nread=15,/overplot
T10, source='B1e14E100_1e9_mu1_env',ls=1, tt=ttshow, nread=15,/overplot
T10, source='B1e14E30_1e9_mu1_env',ls=1, tt=ttshow, nread=15,/overplot
T10, source='B1e14E3_1e9_mu1_env',ls=1, tt=ttshow, nread=15,/overplot
T10, source='B1e14E0.3_1e9_mu1_env',ls=1, tt=ttshow, nread=15,/overplot


	ttshow=1.0
T10, source='B1e14E10_1e9_mu1_env',ls=0, tt=ttshow, nread=15,/overplot
T10, source='B1e14E1_1e9_mu1_env',ls=0, tt=ttshow, nread=15,/overplot
T10, source='B1e14E100_1e9_mu1_env',ls=0, tt=ttshow, nread=15,/overplot
T10, source='B1e14E30_1e9_mu1_env',ls=0, tt=ttshow, nread=15,/overplot
T10, source='B1e14E3_1e9_mu1_env',ls=0, tt=ttshow, nread=15,/overplot
T10, source='B1e14E0.3_1e9_mu1_env',ls=0, tt=ttshow, nread=15,/overplot

endif

	if (0) then begin
	ttshow=0.0001
T10, source='B1e14E100_1e9_mu1',ls=1, tt=ttshow, nread=14
T10, source='B1e14E1_1e9_mu1',ls=1, tt=ttshow, nread=14,/overplot
T10, source='B1e14E10_1e9_mu1',ls=1, tt=ttshow, nread=14,/overplot
T10, source='B1e14E30_1e9_mu1',ls=1, tt=ttshow, nread=14,/overplot
T10, source='B1e14E3_1e9_mu1',ls=1, tt=ttshow, nread=14,/overplot
T10, source='B1e14E0.3_1e9_mu1',ls=1, tt=ttshow, nread=14,/overplot


	ttshow=1.0
T10, source='B1e14E100_1e9_mu1',ls=0, tt=ttshow, nread=14,/overplot
T10, source='B1e14E1_1e9_mu1',ls=0, tt=ttshow, nread=14,/overplot
T10, source='B1e14E10_1e9_mu1',ls=0, tt=ttshow, nread=14,/overplot
T10, source='B1e14E30_1e9_mu1',ls=0, tt=ttshow, nread=14,/overplot
T10, source='B1e14E3_1e9_mu1',ls=0, tt=ttshow, nread=14,/overplot
T10, source='B1e14E0.3_1e9_mu1',ls=0, tt=ttshow, nread=14,/overplot

endif

	E = [0.3,1.0,3.0,10.0,30.0,100.0]	
	for i=0, n_elements(E)-1 do begin
	rho = 9.0+2.0*0.01*dindgen(100)
	rho = 10^rho
	T = 4.2d8 * (rho * 1d-10)^(-0.5) * E[i]^0.6
	oplot, rho, T, col=250, thick=3
	endfor
		
	rho = 9.0+2.0*0.01*dindgen(100)
	rho = 10^rho
	T = 4.2d8 * (rho * 1d-10)^(-0.5) * (10.0*((rho*1d-10)^1.333))^0.63
;	oplot, rho, T, col=250
		
	
		
;	T10,source='B1e14E100.0_1e9',/overplot
;	T10,source='B1e14E30.0_1e9',/overplot
;	T10,source='B1e14E10.0_1e9',/overplot
;	T10,source='B1e14E3.0_1e9',/overplot
;	T10,source='B1e14E1.0_1e9',/overplot
;	T10,source='B1e14E0.3_1e9',/overplot

;	T10,/overplot,nread=14, ls=1
	
	if (0) then begin
		T10,source='B1e15E100.0',ls=2
	T10,source='B1e15E30.0',/overplot,ls=2
	T10,source='B1e15E10.0',/overplot,ls=2
	T10,source='B1e15E3.0',/overplot,ls=2
	T10,source='B1e15E1.0',/overplot,ls=2
	T10,source='B1e15E0.3',/overplot,ls=2
	
	T10,source='B1e14E100.0',/overplot
	T10,source='B1e14E30.0',/overplot
	T10,source='B1e14E10.0',/overplot
	T10,source='B1e14E3.0',/overplot
	T10,source='B1e14E1.0',/overplot
	T10,source='B1e14E0.3',/overplot
	endif

	readcol, 'out/fcontour_14.dat', Pneut2, Tneut2, rhoneut2, format=('X,D,X,X,D,D')
	Pneut2=10.0^Pneut2 * 2.28d14
	Tneut2=10.0^Tneut2
	rhoneut2=10^rhoneut2
	oplot, rhoneut2,Tneut2,linestyle=3, col = 160
	xyouts, 1e8,7e9,textoidl('Neutrino cooling'), col=160, charsize=1.01


	oplot, [3e9,3e9], [1e7,1e10], linestyle=2


	if keyword_set(ps) then close_ps
end


pro times

	T10, nread=15, tt=0.01, ls=2
	T10, nread=15, tt=0.02, ls=0,/overplot
	T10, nread=15, tt=0.1, ls=0,/overplot
	T10, nread=15, tt=1.0, ls=0,/overplot
	T10, nread=15, tt=3.0, ls=0,/overplot
	T10, nread=15, tt=10.0, ls=0,/overplot
	T10, nread=15, tt=30.0, ls=0,/overplot
	T10, nread=15, tt=100.0, ls=0,/overplot
	T10, nread=15, tt=300.0, ls=0,/overplot

end

	
pro T10, delay=delay, png=png, source=source,overplot=overplot, ls=ls, nread=nread, tt=tt

	if not keyword_set(tt) then tt=10.0

	; read Gamma/T for the initial model, used to plot melting curve
	readcol, 'gon_out/grid_profile', ym, Zm, Am, GammaT, format=('X,X,D,X,D,X,D,X,X,D')
	Tm = 1e8 * GammaT/175.0
	TD = 1.4e8 * sqrt(ym*1d-9) * Zm / (0.4*Am)
	readcol, 'out/cv.dat', rhoc, cv, format=('D,D')
	readcol, 'out/nu.dat', rhonu, Tnu, format=('D,D')

	; now look at the detailed profiles as a function of time
	fname='gon_out/out'
	if keyword_set(source) then begin
		fname += '_'+source
	endif
	openr, lun, fname, /get_lun

	; first find out how many grid points
	ngrid=0
	readf, lun, ngrid, format='(I0)'
	print, 'Number of grid points=',ngrid

	multiplot,/reset
	!p.multi=[0,1,1,0,0]

	flag = 1
	
	while (not eof(lun) and flag) do begin

		; get the time
		readf, lun, time

		; read in next batch of data
		if not keyword_set(nread) then nread=13
		data=dblarr(nread,ngrid)
		readf, lun, data
		y=data(0,*)
		T=data(1,*)
		F=data(2,*)
		rho = data(5,*)
		beta=data(12,*)


		if (time gt tt*1.0*24.0*3600.0 and flag) then begin
			print, 'time=', time/(24.0*3600.0)

			; plot upper panel
			;erase
			if keyword_set(overplot) then begin
				oplot, rho, T,linestyle=ls
			endif else begin
				plot, rho, T, /xlog, /ylog, ytitle=textoidl('T (K)'),$
					xtitle=textoidl('\rho (g cm^{-3})'), yrange=[5e7,1e10],ystyle=1, $
					xrange=[1d7,5d11], xstyle=1,linestyle=ls, charsize=1.3
				oplot, ym, Tm, linestyle=2, col=80
				;oplot, ym, TD/3.5, linestyle=4
				oplot, rhoc, cv*1e8, linestyle=3, col=250
			;	oplot, rhonu, Tnu*1e8, linestyle=4
				xyouts, 1.5e11,3e9,textoidl('\Gamma=175')	,charsize=1.01, orientation=20, col=80	
				xyouts, 1.2e11,1.1e9,textoidl('C_{V,e}=C_{V,ion}')	,charsize=1.01, orientation=20, col=250	
			;	xyouts, 1.5e11,1.2e9,textoidl('t_{therm}=t_\nu')	,charsize=1.01, orientation=0	
				xyouts, 1e10,6.5e7,textoidl('E_{25}=0.3,1,3,10,30,100'),charsize=1.2
				;xyouts, 6d10,6d9, textoidl('t=1 day'),charsize=1.1
			endelse
			flag=0
		endif
		
	endwhile

	free_lun,lun

end

	
	
	
	
; ================= Routines to look at output of MCMC runs ==========================
	
pro mc, linearQ=linearQ
	; linearQ => plot linear Q axis rather than log
	readcol, 'out/settle2.mcmc', Tc, Tt, Q, M,R, g, mdot, Qin, format=('D,D,D,D,D,X,D,X,D,D')
	Tc/=1d7
	Tt/=1d8
	if not keyword_set(linearQ) then Q=alog10(Q)
	if not keyword_set(linearQ) then Qin=alog10(Qin)
	
	readcol, 'out/settle2.mcmc.1659_neutrons', Tc0, Tt0, Q0, M0,R0, mdot0, Qin0, format=('D,D,D,D,D,X,X,X,D,D')
	Tc0/=1d7
	Tt0/=1d8
	if not keyword_set(linearQ) then Q0=alog10(Q0)
	if not keyword_set(linearQ) then Qin0=alog10(Qin0)
	
	multiplot,/reset
	!p.multi=[0,1,1,0,0]
	!p.charsize=1.3

	if keyword_set(linearQ) then begin


		plot, Q, Tt, psym=3, $
				ytitle=textoidl('T_t (10^8 K)'), xtitle=textoidl('Q_{imp}'), $
				yrange=[1.0,8.0], ystyle=1, xrange=[0.0,100.0], xstyle=1
;		plot,  Tc, Tt, psym=3, $
;				xtitle=textoidl('T_c (10^7 K)'), ytitle=textoidl('T_t (10^8 K)'), $
;				yrange=[1.0,8.0], ystyle=1, xrange=[4.0,10.0], xstyle=1
;		plot, Q, Qin, psym=3, $
;		 		ytitle=textoidl('Q_{inner}'), xtitle=textoidl('Q_{imp}'), $
;				yrange=[2.0,8.0], ystyle=1, xrange=[0.0,100.0], xstyle=1
	;		plot, Q, mdot, psym=3, $
	 ;			ytitle=textoidl('accretion rate'), xtitle=textoidl('log_{10} Q_{imp}'), $
	;			 xrange=[-3.0,2.0], xstyle=1
				
	endif else begin
		
		plot, M, R, psym=3, $
			ytitle=textoidl('R (km)'), xtitle=textoidl('M (solar masses)'),$
			yrange=[9.0,16.0], ystyle=1, xrange=[1.0,3.0], xstyle=1
	;	oplot, M0[1000:*],R0[1000:*], col=250, psym=3
		plot,  Tc, Tt, psym=3, $
			xtitle=textoidl('T_c (10^7 K)'), ytitle=textoidl('T_t (10^8 K)');, $
			;yrange=[1.0,6.0], ystyle=1, xrange=[3.0,8.0], xstyle=1
	;	oplot, Tc0[1000:*], Tt0[1000:*], col=250,psym=3
		plot, 10^Q, g, psym=3;, $
 		;	ytitle=textoidl('log_{10} Q_{inner}'), xtitle=textoidl('log_{10} Q_{imp}');, $
;			yrange=[3.0,8.0], ystyle=1, xrange=[-3.0,2.0], xstyle=1
		plot, Q, mdot, psym=3, $
 			ytitle=textoidl('accretion rate'), xtitle=textoidl('log_{10} Q_{imp}'), $
			 xrange=[-3.0,2.0], xstyle=1
	;	oplot, Q0, mdot0, col=250,psym=3

	endelse
	
	;plot, M,R, psym=3
	
end
	
		
pro mct
; plots the temperature profile .. this is for MCMC with T(y) as parameters
	nlines=file_lines('out/settle2.mcmc')
	nlines-=2

	openr,lun,'out/settle2.mcmc',/get_lun
	readf,lun,format='(d)',num
	print, num, nlines

	y=dblarr(num)
	readf,lun,y
	print, y

	T=dblarr(num)
	readf, lun, T
	multiplot,/reset
	!p.multi=[0,1,1,0,0]
	!p.charsize=1.3
	plot, y, T, /xlog, /ylog, xrange=[1d10,3d18], xstyle=1, xtitle=textoidl('Column depth (g cm^{-2})'),$
		ytitle=textoidl('Temperature (K)'), linestyle=0, yrange=[1d7,5d8]

	data = dblarr(num,nlines-3)
	readf,lun,data
	Tbar=total(data,2)/(nlines-3)
	T2=total(data^2,2)/(nlines-3)
	sigT=sqrt(T2-Tbar^2)
	oploterror, y, Tbar, sigT,psym=1,/nohat

	print, 'Number of lines read=',nlines

	free_lun, lun
end
		
	
	
pro heat

	readcol, 'heat.dat', E,T,format=('D,D')

	plot, E, T, /xlog,/ylog, xrange=[1d-2,1d3]

end




pro CV, ps=ps

	if keyword_set(ps) then begin
		set_plot, 'ps'
		name = 'cv.ps' 
		device,filename=name,/color
	endif  

	readcol, 'cv.dat', T, e1,e2,e3,e4,e5,e6,e7, format=('D,D,D,D,D,D,D,D')

	efac=60.0*1.67d-24/1.38d-16

	e1*=efac
	e2*=efac
	e3*=efac
	e4*=efac
	e5*=efac
	e6*=efac
	
	T*=1d8

	plot, T, e1, /xlog, xrange=[1d8,1d10], charsize=1.4, xtitle=('T (K)'),$
		ytitle=textoidl('C_V (k_B/60 m_p)'), yrange=[0.5,30.0],ystyle=1,/ylog
	oplot, T,e2,linestyle=1
	oplot, T,e3, linestyle=2
	oplot, T,e4, linestyle=3
	oplot, T,e5, linestyle=4
	;oplot, T,e6, linestyle=1
	;oplot, T,e7, linestyle=3

	e8 = 0.6+0.666*(alog10(T)-9.0)
	e9 = 0.6+1.0*(alog10(T)-9.0)
	e8=10^e8
	e9=10^e9
	oplot, T, e8, col=250

	if keyword_set(ps) then begin
		device,/close
		set_plot,'x'
	endif

end


pro neut, ps=ps

		if keyword_set(ps) then begin
	   	set_plot, 'ps'
		name = 'neut.ps' 
	   	device,filename=name,/color
	  	endif  


	readcol, 'neut.dat', T, e1,e2,e3,e4,e5,e6,e7, format=('D,D,D,D,D,D,D,D')
	
	T*=1d8
	
	plot, T, e1, /xlog,/ylog, xrange=[1d9,1d10], charsize=1.4, xtitle=('T (K)'),$
		ytitle=textoidl('Q_\nu (erg cm^{-3} s^{-1})'), yrange=[1d17,1d25],ystyle=1
	oplot, T,e2,linestyle=1
	oplot, T,e3, linestyle=2
	oplot, T,e4, linestyle=3
	oplot, T,e5, linestyle=4
;	oplot, T,e6
;	oplot, T,e7, linestyle=3
		
	;e8 = 9.04d14 * 100.0 * (T/1d9)^5
	;oplot, T, e8, linestyle=2

	if (0) then begin
	readcol, 'neut_nosynch.dat', T, e1,e2,e3,e4,e5,e6,e7, format=('D,D,D,D,D,D,D,D')
	T*=1d8
	oplot, T, e1, linestyle=1
	oplot, T, e2, linestyle=1
	oplot, T, e3, linestyle=1
	oplot, T, e4, linestyle=1
	oplot, T, e5, linestyle=1
	oplot, T, e6, linestyle=1
	oplot, T, e7, linestyle=1
	endif
	



	
	xx=1.1e9
	xx=4e9
	
	xyouts, 1.2e9, 2d24, textoidl('B=10^{13} G')
	
	ev=1d25/3600.0
	oplot, [1e8,1e10], [ev,ev], linestyle=1
	xyouts, xx,ev*1.3,textoidl('10^{25} erg cm^{-3} in 1 hour')

	ev=1d25/(24.0*3600.0)
	oplot, [1e8,1e10], [ev,ev], linestyle=1
	xyouts, xx,ev*1.3,textoidl('10^{25} erg cm^{-3} in 1 day')
	
	

	T9=T/1d9
	rhofac = 0.5
	eanal = 4d21 * T9 * rhofac^(11.0/3.0) * exp(-3.3*rhofac^(1.0/3.0)/T9)
	oplot, T, eanal, col=250



	if keyword_set(ps) then begin
     	device,/close
     	set_plot,'x'
  	endif
	
end


pro ponsperna

	readcol, 'ponsperna.dat', Lq, dL,ULflag, name,B, format=('D,D,D,A,D')

	print, name,B
	plot,B*1d14, dL*Lq*1d33, /xlog,/ylog, psym=4, charsize=1.5, yrange=[1d35,3d36],ystyle=1,$
	xtitle=('B (G)'), xrange=[3d13,8d15],xstyle=1
	
	
	for i=0,n_elements(name)-1 do begin

		xyouts, 1d14*B[i]*1.1, dL[i]*Lq[i]*1d33, name[i],charsize=1.1
	endfor
	
	
	B=14.0+3.0*0.5*dindgen(3)
	L=36.3-0.9*(B-14.0)
	oplot, 10^B,10^L
	oplot, [1d13,1d14],[2d36,2d36]
	
end


pro fcontour

	print, 'Reading out/grid (this will take a while)'
	readcol, 'out/grid', ya, Ta, Fa, rhoa, epsa,cva, format=('D,D,D,D,D,D')
	print, 'Done..'
	
	plot, [1d19,1d32], [1d7,1d10], /xlog,/ylog,/nodata,xtitle=textoidl('Pressure (cgs)'),$
				ytitle=textoidl('Temperature (K)'), xstyle=1, ystyle=1, charsize=1.6
	
	grav = 2.28d14

	yvals = ya[where(Fa eq 17.0)]

	yy = []
	TT = []

	openw, lun, 'out/fcontour_14.dat',/get_lun, width=500

	foreach y, yvals do begin
				
		F = Fa[where(ya eq y)]
		eps = epsa[where(ya eq y)]
		T = Ta[where(ya eq y)]
		rho = rhoa[where(ya eq y)]
		cv = cva[where(ya eq y)]
		
		ind = where((eps - F) lt 0.0, count)
		if (count lt n_elements(F) and count gt 0 and y gt 10.0) then begin
			j = count-1
			interp = (F[j]-eps[j])/(F[j]-eps[j] - (F[j+1]-eps[j+1]))			
			interp=1.0-interp
			Fout = F[j]*interp + (1.0-interp)*F[j+1]
			epsout = eps[j]*interp + (1.0-interp)*eps[j+1]
			Tout = T[j]*interp + (1.0-interp)*T[j+1]
			rhoout = rho[j]*interp + (1.0-interp)*rho[j+1]
			cvout = alog10(cv[j])*interp + (1.0-interp)*alog10(cv[j+1])			
			printf, lun, j, y, Fout, epsout, Tout, rhoout, (10^y*10^cvout*10^Tout)/10^Fout
			yy=[yy,y]
			TT=[TT,Tout]
		endif
				
	endforeach
		
	free_lun, lun
		
	plot, 10^yy, 10^TT, /xlog,/ylog

	
end

	
	



pro neut4, ps=ps
	
	if keyword_set(ps) then open_ps, 'ps/neut4.ps'

	readcol, 'out/neut4.dat', rho, e1,e2,e3,e4,e5,e6,e7, format=('D,D,D,D,D,D,D,D')

	plot, rho, e1, /xlog,/ylog, xtitle=textoidl('\rho (g cm^{-3})'), yrange=[1d18,1d23],$
	 			ystyle=1, xrange=[1d8,1d14], xstyle=1, charsize=1.6, $
				ytitle=textoidl('Q_\nu (erg cm^{-3} s^{-1})')

;	oplot, rho, e2, linestyle=1
;	oplot, rho, e3, linestyle=2
;	oplot, rho, e4, linestyle=3
	oplot, rho, e6, linestyle=2

;	oplot, rho, e7, col=250

	xyouts, 3d11, 3d19,textoidl('B=10^{14} G'),col=80
	xyouts, 2d12, 8d20,textoidl('B=10^{15} G')



	readcol, 'out/neut5.dat', rho, e1,e2,e3,e4,e5,e6,e7, format=('D,D,D,D,D,D,D,D')

	oplot, rho, e1, linestyle=0, col=80
;	oplot, rho, e2, linestyle=1
;	oplot, rho, e3, linestyle=2
;	oplot, rho, e4, linestyle=3
	oplot, rho, e6, linestyle=2, col=80



	if keyword_set(ps) then close_ps

end



pro econd

	readcol, 'gon_out/initial_condition', rho, sigma, nel,format=('X,X,X,D,X,X,X,X,X,X,X,X,X,D,D')
	
;	plot, rho, sigma, /xlog,/ylog
	plot, nel, sigma, /xlog,/ylog
	
end


