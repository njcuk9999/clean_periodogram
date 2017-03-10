PRO DFOURT, TIME, DATA, FPAR, freq, wfn, dft
;
; This procedure computes the "dirty" discrete Fourier Transform, DFT,
; for a 1-D time series, DATA, which is sampled at arbitrarily spaced
; time intervals given by the TIME vector.  The frequency grid, FREQ, on 
; which the spectral window function WFN and DFT are computed, is controlled
; by the input vector FPAR, which is described below. Note that this
; implementation is completely general, and therefore slow, since it cannot
; make use of the timing enhancements of the FFT.
;
; This IDL implementation of the DFT is based on a suite of FORTRAN
; routines developed by Roberts et al.  For more information concerning
; the algorithm, please refer to:
;   Roberts, D.H., Lehar, J., & Dreher, J. W. 1987, AJ, 93, 968
;  "Time Series Analysis with CLEAN. I. Derivation of a Spectrum"
;
; INPUT PARAMETERS:
;    TIME  - input time (independent) vector
;    DATA  - input dependent vector  
;    FPAR  - fltarr(3), vector controlling output frequency grid
;            ENTER 0 TO SELECT DEFAULT VALUES
;               fpar(0) -> frequency increment for the FT  [default: 1/T]
;               fpar(1) -> max. frequency in the FT        [default: 1/min(dt)]
;               fpar(2) -> points per beam                 [default: 4]
;
; OUTPUT PARAMETERS:
;    freq  - frequency vector
;    wfn   - spectral window function           [complex]
;    dft   - "dirty" discrete Fourier transform [complex]
;
; OTHER ROUTINES CALLED:
;    cl_fpar, meanvar
;
; HISTORY: 
; Jan. 91:  translated to IDL                       [A. W. Fullerton, BRI]
;           source: FORTRAN code distributed by Roberts et al.
; Apr. 96:  vectorized code for efficiency          [AWF, USM]
;           reorganized I/O list, cleaned up documentation 
;           incorporated subroutine cl_dfour directly
;============================================================================
    twopi=6.2831853D0       
    
    if (n_params(0) lt 3) then begin
        print,' '
        print,' DFOURT, TIME, DATA, FPAR, freq, wfn, dft '
        print,' '
        return
    endif

; Find the mean time and mean data value. Subtract these "DC" terms from
; their respective vectors.
    meanvar,time,ndat,tmean
    tvec=time-tmean
    meanvar,data,ndat,dmean
    dvec=data-dmean
    
; Compute default frequency grid values, and substitute them if the
; values supplied by the user are 0.
    cl_fpar,tvec,fpar_def,/silent
    a=where(fpar_def eq -1.)
    if (a(0) ne -1) then begin
       print,' '
       print,' ** ERROR in DFOURT: something wrong with the time vector'
       print,' '
       return
    endif
    a=where(float(fpar) eq 0.)
    if (a(0) ne -1) then fpar(a)=fpar_def(a)

    dfreq=fpar(0)/fpar(2)           ; size of frequency increment
    mfreq=long(fpar(1)/dfreq) + 1   ; number of frequencies

; Dimension the output arrays. Compute the frequency vector.
    freq = dindgen(2*mfreq)*dfreq   ; frequency vector
    wfn  = complexarr(2*mfreq)      ; spectral window function
    dft  = complexarr(mfreq)        ; "dirty" discrete Fourier Transform

; Calculate the Fourier transform frequency by frequency. The DFT is 
; normalized to have the mean value of the data at zero frequency.
    for i=0L, mfreq-1 do begin
        phase  = -twopi*freq(i)*tvec
	phvec  = complex(cos(phase),sin(phase))
	wfn(i) = total(phvec)/ndat
	dft(i) = total(dvec*phvec)/ndat
    endfor
    
    for i=mfreq,2*mfreq-1 do begin   ; Complete the spectral window function
        phase  = -twopi*freq(i)*tvec    
	phvec  = complex(cos(phase),sin(phase))
	wfn(i) = total(phvec)/ndat
    endfor
    
   return                            ; All done!
   end

