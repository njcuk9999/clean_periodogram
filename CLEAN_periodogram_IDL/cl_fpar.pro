PRO CL_FPAR,TVEC,fpar,nfreq,silent=silent

; This procedure uses the input vector TVEC (which is the independent 
; variable associated with a time series of measurements) to compute the
; *default* parameters that govern the frequency grid that will subsequently
; be used by DFOURT to calculate the discrete Fourier Transform of the 
; the time series.  The defaults are described below.
;
; INPUTS:
;   TVEC   - time vector for data set under consideration.
;            *MUST INCREASE STRICTLY MONOTONICALLY*
;
; OUTPUTS:
;   fpar   - fltarr(3) - the default values governing the frequency grid
;            for the calculation of the discrete Fourier Transform 
;              fpar(0) = frequency resolution       [default: 1/T]
;              fpar(1) = maximum frequency          [default: 1/min(dt)]
;              fpar(2) = points per restoring beam  [default: 4]
;            Note the the frequency resolution element fpar(0) is oversampled
;            by PPB = fpar(2) to ensure accurate determination of the location
;            of peaks in the Fourier Transform.
;
;   nfreq  - the number of frequencies in the default frequency grid
; 
; KEYWORDS:
;   silent - if set, the output message listing the values of fpar will
;            be suppressed
;
; OTHER ROUTINES CALLED:
;   None.
;
; HISTORY:
;  Sept. 1994:  Written by A. W. Fullerton, Uni-Sternwarte Muenchen
;  April 1996:  Upgraded error checking and documentation  [AWF, USM]
;=============================================================================
    if (n_params(0) lt 1) then begin
       print,' '
       print,' CL_FPAR, TVEC, fpar, nfreq, /silent   '
       print,' NOTE:    fpar = [df, fmax, ppb]       '
       print,' '
       return
    end

; Check that the time vector is monotonically increasing
    nt=n_elements(tvec)          ; number of elements in TVEC    
    check=findgen(nt)            ; test 
    order=sort(tvec)             ; indices sorted into ascending order
    diff=order-check             ; should be zero if monotonically increasing
    a=where(diff ne 0.)          
    if (a(0) ne -1) then begin   
       print,'** ERROR in CL_FPAR: TVEC is not monotonically increasing'
       print,'** Returning with negative fpar  -->'
       fpar=[-1,-1,-1]
       return
       print,' '
    endif   

; Frequency resolution: default is 1./(total time spanned)
    dt=max(tvec)-min(tvec)
    if (dt eq 0.) then begin
       print,'** ERROR in CL_FPAR: total time spanned is 0. '
       print,'** Returning with negative fpar -->'
       fpar=[-1,-1,-1]
       return
       print,' '
    endif   
    fres=1./dt

; Maximum frequency: default is 1./(2.*[minimum time separation])
    tlo=tvec(0:nt-2)
    thi=tvec(1:nt-1)
    tmin=min(thi-tlo)
    if (tmin le 0.) then begin
       print,'** ERROR in CL_FPAR: minimum time different is le 0. '
       print,'** Returning with negative fpar -->'
       fpar=[-1,-1,-1]
       return
       print,' '
    endif       
    fmax=1./(2.*tmin)
    
; Points per beam: default is 4
    ppb=4.
    
; Total number of frequencies
    nfreq=long(fmax/(fres/ppb))+1

    if not keyword_set(silent) then begin 
       print,'  '
       print,'** CL_FPAR: default frequency parameters are - '
       print,'   df    = '+strtrim(string(fres),2)
       print,'   fmax  = '+strtrim(string(fmax),2)
       print,'   ppb   = '+strtrim(string(ppb),2)
       print,'   nfreq = '+strtrim(string(nfreq),2)
       print,'  '
    endif

    fpar=[fres, fmax, ppb]      ; output vector

    return
    end

