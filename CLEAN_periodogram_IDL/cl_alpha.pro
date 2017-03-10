FUNCTION CL_ALPHA, WFN, DFT, L
;
; This function returns an estimate for the component ALPHA, which 
; produces the DFT at frequency index L via the relation:
;
;       DFT(L) = ALPHA*WFN(0) + CONJ(ALPHA) * WFN(2L)
;
; ALPHA is then given by:
; 
;                  DFT(L) - CONJ(DFT(L))*WFN(2L)
;       ALPHA   =  --------------------------------          
;                      1  - ABS( WFN(2L) )^2
;
; See Section III b) [especially equation (24)] of Roberts et al. 
; (1987, AJ, 93, 968).
;
;
; INPUTS: 
;    WFN  - complexarr(2*NFREQ), the spectral window function
;    DFT  - complexarr(  NFREQ), the "dirty" discrete Fourier Transform
;    L    - integer            , the frequency index of the desired component
;
; OUTPUTS:
;   The amplitude alpha is returned by the function.
;
; OTHER ROUTINES CALLED:
;   None.
;   
; HISTORY:
;  Jan. 91: translated for FORTRAN code by Roberts et al.  [AWF, Bartol]
;  Apr. 96: recoded for efficiency and added documentation [AWF, USM]
;=============================================================================
    err=1.0D-4                         ; allowed error in WNORM
                                                                              
; Find the (L, -L) components which produce DFT(L) through WFN
    win2l = wfn(2*l)                   ; (L,-L) interference              
    wnorm = 1.0 - abs(win2l)^2

; Trap to avoid singularities 
    if (wnorm lt err) then alpha = 0.5*dft(l)   $
                      else alpha = (dft(l) - win2l*conj(dft(l)))/wnorm

    return, alpha                      ; return with estimate of alpha
    end
