FUNCTION CL_BEAM, WFN, b_sigma 
; 
; This function returns the "CLEAN beam", which will subsequently be 
; convolved with the CLEAN components to obtain the CLEANed Fourier 
; Transform.  The beam is modelled by a Gaussian function truncated 
; at +/- 5*b_sigma, where b_sigma is the standard deviation determined
; from the HWHM of the amplitude of the peak at 0-frequency in the spectral 
; window function WFN.  Since the time vector was symmetrized by subtraction 
; of the mean value (so that t=0 corresponds to the middle point of the 
; time series) prior to the calculation of WFN by the routine DFOURT, the 
; beam will be a purely real entity.  For details, refer to  Roberts et al. 
; 1987, AJ, 93, 968.
;
; INPUTS:
;  WFN     - the spectral window function computed by DFOURT [complex]
;
; OUTPUTS:
;  beam    - the Gaussian restoring beam, returned as a COMPLEX entity 
;            although the imaginary parts are all 0.  
;
; NOTE: in contrast to the "beam" computed by the FORTRAN FILLB routine 
;       by Roberts et al., this beam is "two-sided"; i.e., it runs from
;       (-5*b_sigma, +5*b_sigma), not just from (0, +5*b_sigma).  The
;       two-sided beam can be used directly by the intrinsic IDL procedure
;       CONVOLV, which is much more efficient than the 
;
; OPTIONAL OUTPUTS:
;  b_sigma - the sigma of the Gaussian restoring beam, in "index units"
;            of the frequency vector.
;
; OTHER ROUTINES CALLED:
;  None.
;
; HISTORY:
;  Apr. 1996:  Written by A. W. Fullerton [Uni-Sternwarte Muenchen] 
;              Based on FORTRAN routines of Roberts et al.
;              Concatenation of routines CL_HWHM and CL_FILLB
;=============================================================================
; First, estimate the half width at half max of the spectral window function.  
; Assume that the maximum is located at the 0th element.
    hmax=0.5*abs(wfn(0))
    a=where(abs(wfn) le hmax)
    i2=a(0)
    w2=abs(wfn(i2))
    if (w2 lt hmax) then begin           ; linearly interpolate to get a
       i1=a(0)-1                         ;          more accurate estimate
       w1=abs(wfn(i1))
       q=(hmax - w1)/(w2-w1)             ; interpolating function 
       hwidth=i1 + q*(i2-i1)             ; HWHM of WFN (frequency index units)
    endif else begin
       hwidth=i2                         ; HWHM of WFN (frequency index units)
    endelse

; Now fill the restoring beam with a Gaussian characterized by HWHM = hwidth
; (in frequency index units).  The beam will be a purely real function, 
; since the time vector was symmetrized (tvec = tvec - tmean) prior to the
; calculation of WFN by DFOURT.
    b_sigma= hwidth/sqrt(2.*alog(2.))    ; beam sigma in index units
    const=1./(2.*b_sigma*b_sigma)        ; Gaussian normalization constant
    n_beam=fix(5*b_sigma)+2              ; size of restoring beam
    x=findgen(n_beam)                    ;     [truncated at 5*b_sigma]  
    y=x
    y(0)=complex(exp(-const*x*x),0.)     ; "one-sided" beam, purely real

    mb=n_beam-1                          ; max index in y
    beam=complexarr(2*n_beam -1)         ; define beam 
    beam(0)=reverse(y)                   ; -tive freq: no need to take conj 
    beam(mb+1:*)=y(1:mb)                 ; +tive freq
    
    return, beam
    end
    
