PRO PK_FIT, FREQ, CDFT, CCOMP, FPEAK, B_SIGMA, cos_par, $
            thresh=thresh, nclip=nclip, nscan=nscan,    $
	    plot=plot

; This procedure estimates the parameters of the cosine function 
; (cos_par = [amplitude, frequency, phase]) for the highest peak in the 
; CLEANed Fourier Transform CDFT nearest the input frequency FPEAK.  These
; estimates are made from the CLEAN components, CCOMP.  CDFT and CCOMP
; are computed on the frequency grid FREQ, and CDFT has been convolved
; with a Gaussian restoring beam with sigma = B_SIGMA.   One-sigma 
; uncertainties for the frequency are estimated by means of the "post-mortem"
; analysis described by Schwarzenberg-Czerny (1991, MNRAS, 253, 198).
; These errors should be treated cautiously, especially for weak peaks,
; though they are likely more reliable than the usual alternatives; 
; see Schwarzenberg-Czerny for illuminating discussion.
;
; Why is this routine so involved?  Although the frequency and phase of a
; cosine signal can be determined quite accurately from the CLEANed DFT 
; itself, the amplitude cannot be recovered precisely from it.  This is 
; particularly true in cases where the peak falls between frequency grid 
; points, since the convolution of the CLEAN components with the Gaussian
; restoring beam smears the amplitude of individual components (which 
; otherwise add linearly) together.  This routine circumvents these 
; difficulties, at the price of having to retain the CLEAN components.
;
; INPUTS: 
;  FREQ    - frequency vector
;  CDFT    - CLEANed Discrete Fourier Transform   [complex]
;  CCOMP   - individual CLEAN components          [complex]
;  FPEAK   - approximate location of the frequency peak
;            An interval corresponding to NSCAN*B_SIGMA (usually: 3*B_SIGMA) 
;            will be scanned in order to locate the highest peak in the 
;            vicinity of FPEAK.
;  B_SIGMA - sigma associated with Gaussian restoring beam, in index units
;             
; OUTPUTS:
;  cos_par - fltarr(4), the values of the cosine parameters associated 
;            with the input frequency estimated from DFCT
;              cos_par(0) = amp_est = semi-amplitude of the cosine component
;              cos_par(1) = f_est   = frequency of the cosine component
;              cos_par(2) = ph_est  = phase of the cosine component [-pi, pi]
;              cos_par(3) = f_sig   = estimated 1-sigma uncertainty in f_sig
;
; OPTIONAL KEYWORDS: [Default values are in square brackets]
;  thresh = [0.0] - This keyword can be used to set an *amplitude* threshold
;                   below which any peaks are deemed to be unreliable or 
;                   spurious. The threshold is specified in units of the 
;                   continuum, so that a typical value might be 0.01 (i.e., 
;                   1% of F_c).  The cosine parameters will still be returned 
;                   for peaks below this threshold, but the estimated 
;                   frequency will be made negative in order to indicate 
;                   "below threshold".  The Schwarzenberg-Czerny "post-mortem" 
;                   procedure to estimate 1-sigma uncertainties will *not*
;                   be performed.  
;
;  nscan  = [3]   - the number of "beam sigmas" to scan on either side of 
;                   the input fpeak when searching for the local maximum
; 
;  nclip  = [2]   - the number of iterative 3-sigma clippings to perform when
;                   estimating the background noise level around the peak in 
;                   the power spectrum (required to estimate f_sig) 
;
;  plot           - if set, a diagnostic plot will be produced, provided 
;                   the peak in question is above the specified threshold
;
; OTHER ROUTINES CALLED:
;  meanvar, pk_gfit
;
; HISTORY;
;  Apr. 1996: Written and tested  [A. W. Fullerton, Uni-Sternwarte Muenchen]
;  July 1996: Added error trap in case total(abs(CDFT))=0.    [AWF, USM] 
;=============================================================================
    twopi=6.2831853d0
    if (n_params(0) lt 5) then begin
       print,' '
       print,' PK_FIT, FREQ, CDFT, CCOMP, FPEAK, B_SIGMA, cos_par, $ '
       print,'         thresh = [0.00], nclip = [2],  nscan=[3],     '
       print,'         /plot                                         '
       print,' NOTE:   cos_par = [amp_est, f_est, ph_est, f_sig]     '
       print,' '
       return
    endif
    if not keyword_set(thresh) then thresh=0.0
    if not keyword_set(nclip) then nclip=2
    if not keyword_set(nscan) then nscan=3
    dfreq=freq(1)                            ; frequency increment
    b_sig=b_sigma*dfreq                      ; beam sigma in frequency units

; Check that the peak can be isolated...
    a=where((freq ge fpeak-dfreq/2.) and (freq le fpeak+dfreq/2.))
    i_peak=a(0)                              ; index of peak
    if (i_peak eq -1) then begin
        print,' '
	print,'** ERROR in PK_FIT: cannot isolate specified frequency peak'
	print,'** Returning -->'
        print,' '
	cos_par=[-1.,-1.,-1.,-1.]
        return
    endif 

; Define the frequency range to scan for the local maximum. 
    freq_max=fpeak + nscan*b_sig
    freq_min=fpeak - nscan*b_sig
    a_freq=where((freq ge freq_min) and (freq le freq_max))       

; Scan for peak in CDFT...
    fvec=freq(a_freq)
    cvec=abs(cdft(a_freq))
    imax=where(cvec eq max(cvec))
    n_fpeak=fvec(imax(0))

; Recenter on new peak...consider only +/- 3*b_sig for estimation...
    freq_max=n_fpeak + 3.0*b_sig 
    freq_min=n_fpeak - 3.0*b_sig        
    a_freq=where((freq ge freq_min) and (freq le freq_max))       

; Now work with the CLEAN components directly...
    fvec=freq(a_freq)
    cvec=ccomp(a_freq)
    r_cvec=float(cvec)      ; real part of CLEAN components within +/-3 sigma
    i_cvec=imaginary(cvec)  ; imaginary part...                 |  of maximum
    w=abs(cvec)                                             ; weights 
    norm=total(w)                                           ; sum of weights
    if (norm ne 0.) then begin
        f_est=total(w*fvec)/norm                            ; frequency    
        amp_est=2.0*sqrt(total(r_cvec)^2 + total(i_cvec)^2) ; amplitude	
        ph_est=atan(total(i_cvec),total(r_cvec))            ; phase	
    endif else begin
        f_est=0.                                            ; frequency
	amp_est=0.                                          ; amplitude
	ph_est=-0.                                          ; phase
   endelse

; ESTIMATE UNCERTAINTY IN F_EST
; Use the method described by Schwarzenberg-Czerny (1991) to estimate
; the 1-sigma uncertainty from the height of the POWER SPECTRUM peak
; above the background noise.  Strictly speaking, this approach is only 
; valid for large, "high S/N" peaks, but we persist for all peaks above
; the specified threshold, with the caveat "let the user beware"!  
; The "sigma-clip" routine to estimate the noise level near the peak in 
; question is a bit of a kludge, particularly in regions where multiple peaks 
; or incompletely CLEANed features from the window function remain.  In 
; these circumstances, the uncertainty will be overestimated.
;
; OUTLINE OF STEPS
;   a) Estimate the background noise level, N; iterate NCLIP times
;   b) Model the POWER peak that WOULD be observed if a component
;      with amplitude f_est were convolved with a Gaussian beam
;      of standard deviation b_sig.  This approach is preferable to
;      Gaussian fitting because: (a) accurate Gaussian fitting requires 
;      f_est to fall on a frequency grid point; and (b) Gaussian fitting
;      will fail for weak peaks, multiple peaks, or peaks with complex
;      backgrounds (due, e.g., to incomplete removal of window structure).
;   c) Invert the Gaussian model for the POWER peak to determine its full 
;      width at the PMAX-N level.  This corresponds approximately to the
;      1-sigma error on f_est.
;
; NOTE: 
; The Gaussian parameters for a model of the *amplitude* peak are:
;     A(x) = a0 * exp { -0.5 * [ (x - a1) / a2 ] }
; where  a[0, 1, 2] = gpar[0, 1, 2] = [amplitude, location, std.dev]
; The *power* peak associated with this will also be a Gaussian:
;     P(x) = p0 * exp { -0.5 * [ (x - p1) / p2 ] } 
;          = A^2
; where p0 = a0^2, p1 = a1, and p2 = a2/sqrt(2).
; Hence, our model for the POWER PEAK, convolved with the restoring beam, is:
;       gpar=[amp_est^2, f_est, b_sig/sqrt(2)]
;----------------------------------------------------------------------------

    if ((amp_est ge thresh) and (f_est ne 0.)) then begin

;   Gaussian parameters for model of POWER peak
       gpar=[amp_est^2, f_est, b_sig/sqrt(2)] 

;   Now estimate the noise in the background around the peak.
;   Sample background in beamwidth to right (higher frequencies) of peak
       bplus_min = f_est + 3.0*b_sig       
       bplus_max = f_est + 9.0*b_sig
       aplus = where((freq ge bplus_min)  and (freq le bplus_max))

;   Sample background in beamwidth to left (lower frequencies) of peak
       bminus_min= f_est - 9.0*b_sig
       bminus_max= f_est - 3.0*b_sig
       aminus= where((freq ge bminus_min) and (freq le bminus_max))

;   Use background vector to estimate noise; iterate NCLIP times
;   Remember to DOUBLE amplitude estimates to allow for negative frequencies!
       if (aminus(0) ne -1) then begin
          backvec=2.0*[abs(cdft(aminus)),abs(cdft(aplus))]
       endif else begin
          backvec=2.0*[abs(cdft(aplus))]
       endelse
       meanvar,backvec,n,back,v
       for i=1,nclip do begin
           a=where(abs(backvec-back) le 3.*sqrt(v))
 	   backvec=backvec(a)
	   meanvar,backvec,n,back,v
       endfor
       back=back^2

; Invert the Gaussian model to find the frequency width of the POWER peak at
; the PMAX - BACK leve.  This is the 1-sigma confidence limit on f_est.
       nsr=back/gpar(0)             ; noise to signal power ratio
       if (nsr le 1.) then begin
       
          f_sig=2.0*gpar(2)*sqrt(-2.0*alog(1.0-back/gpar(0)))

	  if (f_sig lt 0.) then begin
	     f_est=-f_est
	     f_sig=0.
	  endif
	  
       endif else begin              ; something has gone horribly wrong
       
           f_est=-f_est              ; flag this case with a negative amp_est!
	   amp_est=-amp_est
	   f_sig=0.
	   
       endelse

; OPTION to produce a diagnostic plot
       if keyword_set(plot) then begin
          data_min  = f_est - 5.0*b_sig
          data_max  = f_est + 5.0*b_sig
          adata = where((freq ge data_min) and (freq le data_max)) 
          xvec=freq(adata)
          xx=findgen(201)/200*(max(xvec)-min(xvec))+min(xvec)
	  pk_gfit,xx,gpar,dfit
          ymax=1.2*max(dfit)
	  ymin=-0.1*ymax
          aplot=[aminus,adata,aplus]  
          plot,freq(aplot),(2.*abs(cdft(aplot)))^2, psym=1, $
                 xstyle=1,                                  $
		 xticklen=-0.02,                            $
    	         xtitle='!6Frequency',                      $
		 ystyle=1,yrange=[ymin,ymax],               $
		 yticklen=-0.01,                            $
                 ytitle='!6 Power',                         $
                 charsize=1.4
	    oplot,xx,dfit
	    oplot,freq(aplot),(2.0*abs(ccomp(aplot)))^2,linestyle=3
	    oplot,freq(aplot),freq(aplot)*0.+back,linestyle=2
	    oplot,freq(aminus),freq(aminus)*0.+0.2*ymax
	    oplot,freq(aplus),freq(aplus)*0.+0.2*ymax	    
	    oplot,[f_est,f_est],[ymin,max(dfit)],linestyle=0
	    oplot,[f_est,f_est]-f_sig,[ymin,max(dfit)],linestyle=1
	    oplot,[f_est,f_est]+f_sig,[ymin,max(dfit)],linestyle=1	    
	    lab0='    !6PK_FIT:'
	    lab1='    !6Freq = ' +strtrim(string(f_est),2) +' !9+!6 '
	    lab1=lab1+strtrim(string(f_sig),2)
	    lab2='    !6Back = '+strtrim(string(back),2)
	    xyouts,min(freq(aplot)),1.15*max(dfit),lab0,size=1.2
	    xyouts,min(freq(aplot)),1.10*max(dfit),lab1,size=1.2
	    xyouts,min(freq(aplot)),1.05*max(dfit),lab2,size=1.2
	    
	    ans=''
            read,'-> Press <CR> to continue...',ans
	    
         endif

    endif else begin             ; negative frequency indicates        
                                 ; "below threshold" 
       f_sig=0.0                     
       f_est=-f_est               
       
   endelse
   
   cos_par=[amp_est, f_est, ph_est, f_sig]

   return
   end

