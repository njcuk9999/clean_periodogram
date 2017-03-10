PRO CL_SUBCMP, WFN, DFT, CCOMP, L
;                                                                              
; This procedure computes and removes the influence of CLEAN component CCOMP 
; at frequency index L on the "dirty" DFT.  It proceeds as follows:
;  a) The spectral window function WFN is aligned with the CLEAN component 
;     in "index space".
;  b) The realigned WFN is scaled by CCOMP for positive frequencies
;     and CONJ(CCOMP) for negative frequencies.
;  c) The realigned, rescaled version of the WFN is subtracted from the DFT.
; The realigned, rescaled version of the WFN contains all the frequency
; couplings that arise in DFT due to the interaction of the CLEAN component 
; CCOMP with the spectral window function WFN.  Thus, the influence of 
; this component can be removed from DFT by a simple subtraction.  However,
; since only a fraction (given by the gain parameter) of the full CCOMP at 
; a given frequency is in general removed at a time, the process usually has 
; to be iterated.
;
; INPUTS:  
;  WFN   - COMPLEXARR(2*NFREQ) - the spectral window function
;  DFT   - COMPLEXARR(  NFREQ) - the dirty Discrete Fourier Transform
;  CCOMP - COMPLEX             - the CLEAN component to remove
;  L     - INTEGER             - the frequency index of the CLEAN Component
;
; OUTPUTS:
;  DFT is **modified** via subtraction of the the frequency couplings 
;  associated with the CLEAN component CCOMP at frequency index L.
;
; OTHER ROUTINES CALLED:
;  None.
;
; HISTORY:
;  Jan. 91: translated from FORTRAN code by Roberts et al.   [AWF, BRI]
;  Apr. 96: vectorized and upgraded documentation            [AWF, USM]
;
; NOTE: [AWF, April 22, 1996] 
; This version of CL_SUBCMP does not make the basic mathematical operations 
; as clear as the previous version, but is nonetheless preferable because it 
; avoids do loops and extra calls to functions.  For completeness, and also 
; to help illuminate the present version of the code, here is a stipped 
; down version of the old procedure:
;     PRO CL_SUBCMP, WFN, DIRTY, L, CCOMP
;         cneg=conj(ccomp)
;         for i=0L,n_elements(dirty)-1 do begin
;             dirty(i)=dirty(i) - ccomp*ncl_cval(wfn,i-l) - cneg*wfn(i+l)
;         endfor
;         return
;         end
;
; where CL_CVAL is a function with the following characteristics:
;     FUNCTION NCL_CVAL, ARRAY, I                                   
;     This function returns CVAL, the value of the complex ARRAY at the index 
;     location I, subject to the following rules:
;     (a) if  0 < I < n_elements(array)-1, CVAL = ARRAY(I)
;     (b) if      I > n_elements(array)-1, CVAL = 0.0
;     (c) if  I < 0                      , CVAL = CONJ(ARRAY(I))
;     Item (c) implies that ARRAY is Hermitian, since ARRAY(I) = CONJ(ARR(-I)).
;           ii=abs(i)
;           if (ii le n_elements(array)-1) then begin ; i is defined
;	    if (i ge 0) then cval=array(ii) $         ; i > 0, take array(i)
;	                else cval=conj(array(ii))     ; i < 0, take conjugate
;           endif else begin                          ; i is undefined
;              cval=complex(0.,0.)         
;           endelse
;           return, cval
;           end
;=============================================================================
    czero = complex(0.,0.)               ; define once and for all
    nwind = n_elements(wfn)              ; max element = nwind - 1 
    ndirt = n_elements(dft)              ; max element = ndirt - 1
    
; Compute the effect of +L component  
    cplus = complexarr(ndirt)            
    index = lindgen(ndirt) - l           ; index for wfn, shifted to +l comp
    a=where(index lt 0.)                          ; take conj(wfn(index(a))) 
    b=where((index ge 0.) and (index le nwind-1)) ; take wfn(index(b)) 
    c=where(index ge nwind)                       ; set to czero 
    if (a(0) ne -1) then cplus(a) = conj(wfn(abs(index(a))))    
    if (b(0) ne -1) then cplus(b) = wfn(index(b))
    if (c(0) ne -1) then cplus(c) = replicate(czero,n_elements(c))
    
; Compute effect of -L component
    cminus = complexarr(ndirt)
    index  = lindgen(ndirt) + l          ; index for wfn, shifted to -l comp 
    a=where(index lt 0.)                          ; take conj(wfn(index(a))) 
    b=where((index ge 0.) and (index le nwind-1)) ; take wfn(index(b)) 
    c=where(index ge nwind)                       ; set to czero     
    if (a(0) ne -1) then cminus(a) = conj(wfn(abs(index(a))))
    if (b(0) ne -1) then cminus(b) = wfn(index(b))
    if (c(0) ne -1) then cminus(c) = replicate(czero,n_elements(c))    

; Return the realigned, rescaled window function.  The effect of
; CCOMP can subsequently be removed from DFT by subtracting "subcmp"
    dft = dft - ccomp*cplus - conj(ccomp)*cminus
    
    return
    end
