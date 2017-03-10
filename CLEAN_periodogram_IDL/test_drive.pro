PRO TEST_DRIVE
;
; This routine demonstrates the use of the routines DFOURT, CLEAN, and PK_FIT.
;
; INPUTS:
;  None.
;
; OUTPUTS:
;  None.  Previously prepared test data with two known cosinusoidal signals
;  (but no noise) are analyzed by means of the CLEAN routines, and the 
;  estimated parameters associated with the signals (i.e., amplitude, 
;  frequency, phase) are compared with the input values.
;
; RESTRICTIONS:
;  The routine assumes that the current plotting device is X-windows.
;  Graphical output will be suppressed if this is not the case.
;
; OTHER ROUTINES CALLED:
;  cl_beam, clean, dfourt, rdfil, pk_fit
;
; HISTORY:
;  April 1996: Written [A. W. Fullerton, Uni-Sternwarte Muenchen]
;=============================================================================
; Read in the test data.  This file was originally prepared by Roberts et al., 
; and distributed with their FORTRAN version of the CLEAN routines that are
; the basis of this IDL version.
    cospar1_in=[1.0, 31.,  0.25]            ; [amp, freq, phase] for signal 1
    cospar2_in=[0.4, 57., -0.25]            ; [amp, freq, phase] for signal 2
    rdfil,'test.dat',2,buff,/silent         ; read the data file
    tvec=reform(buff(0,*))                  ; extract time vector  
    dvec=reform(buff(1,*))                  ; extract data vector

; Test that current plotting device is "X"; if so, plot the data.
    if (strupcase(!d.name) ne 'X') then gplot=0 else gplot=1
    if (gplot eq 1) then begin
        window,0
        plot,tvec,dvec,psym=-1,                           $ 
             /xstyle, xtitle='Time',                      $
	     ytitle='Data',                               $
	     title='Input Time Series: Sum of 2 Cosines', $
   	     charsize=1.4
    endif
    
; Compute the "dirty" discrete Fourier transform.  Use the default frequency 
; parameters to describe the frequency grid.  The defaults are selected by 
; entering zeroes.
    fpar=[0.,0.,0.]                         ; fpar = [df, fmax, ppb]
    dfourt,tvec,dvec,fpar,freq,wfn,dft      ; compute window function and DFT

; Clean the DFT.  For this demonstration, use a gain of 0.5 and continue
; for 100 iterations.
    gain=0.5                                ; define "gain" per iteration
    ncl=100                                 ; define number of iterations
    clean,freq,wfn,dft,gain,ncl,cdft,ccom   ; clean the time series

; Plot the CLEANed amplitude spectrum.  The factor of 2 allows for
; the "mirror image" of the DFT at negative frequencies.
    if (gplot eq 1) then begin
       window,1
       plot,freq(0:n_elements(cdft)-1),2.0*abs(cdft),  $
            /xstyle, xtitle='Frequency',               $
	    ytitle='Amplitude',                        $
	    title='CLEANed Amplitude Spectrum',        $
            charsize=1.4
    endif

; The next few lines show how to estimate paramters associated with any 
; peaks in the CLEANed Fourier transform.  The CLEAN components themselves
; are the key inputs to the routine pk_fit.
     f1=30.                                 ; *approximate* frequency of peak 1
     f2=58.                                 ; *approximate* frequency of peak 2
     beam=cl_beam(wfn,b_sig)                ; get the sigma of the "beam"
     pk_fit,freq,cdft,ccom,f1,b_sig,cospar1 ; get cosine parameters for peak 1
     pk_fit,freq,cdft,ccom,f2,b_sig,cospar2 ; get cosine parameters for peak 2

; Compare input results with those estimated from CLEANed Fourier transform.
; Remember that this is an ideal case: uniform, well-sampled data and no noise!
     fmt='(2X, A11, 6X, 3(F8.4,8X))'
     print,'  '
     print,' ================================================================'
     print,'                    Amplitude      Frequency        Phase / Pi   '
     print,' ---------------------------------------------------------------- '
     print,'  Signal 1 - '
     print,'Input:     ',cospar1_in(0:1),cospar1_in(2),format=fmt   
     print,'Recovered: ',cospar1(0:1),cospar1(2)/(2.*!pi),format=fmt
     print,' '
     print,'  Signal 2 - '
     print,'Input:     ',cospar2_in(0:1),cospar2_in(2),format=fmt   
     print,'Recovered: ',cospar2(0:1),cospar2(2)/(2.*!pi),format=fmt     
     print,' ---------------------------------------------------------------- '
     print,'  '     
     print,'  '     

     return
     end
