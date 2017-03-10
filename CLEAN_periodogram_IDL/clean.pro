PRO CLEAN, FREQ, WFN, DFT, GAIN, NCL, cdft, ccomp, rdft
;
; This procedure deconvolves the spectral window function WFN from the 
; "dirty" discrete Fourier Transform DFT by using a 1-D version of the 
; iterative CLEAN algorithm [Hoegbom, J. A. 1974, A&AS, 15, 417].  Through 
; successive subtractions of fractions of WFN, the couplings between physical 
; frequencies and their aliases or pseudo-aliases will be removed, while 
; preserving the information contained in the peak (i.e., the frequency, 
; amplitude, phase of the cosinusoidal component).
;
; The routine proceeds as follows:
;  1) During each iteration, the WFN is centered on the peak that currently 
;     has the  maximum amplitude.  
;  2) A fraction of the amplitude of this peak, specified by the GAIN 
;     parameter, is entered into the complex array of "CLEAN components", 
;     CCOMP, at the location corresponding to the peak.  Eventually, the 
;     entire amplitude of the peak will be restored to this CLEAN component.
;  3) The WFN is scaled by the current CCOMP, and subtracted from the input 
;     DFT to form a residual Fourier spectrum, RDFT.  
;  4) The process is repeated, with the RDFT from the previous iteration used
;     as the input spectrum to the current iteration.  
;  5) After NCL iterations, CCOMP is convolved with the Gaussian "beam",
;     truncated at 5*b_sigma.  The standard deviation, b_sigma, is determined
;     from the HWHM of the amplitude of the primary (0-frequency) peak of 
;     WFN.  In this suite of routines, the beam is a purely real function 
;     because the time vector is symmetrized [i.e., tvec = tvec - tmean] 
;     prior to the computation of WFN and DFT by the procedure DFOURT.
;  6) The residual Fourier transform is added to this convolution to produce 
;     the CLEANed discrete Fourier Transform, CDFT.
;
; Since the Fourier Transforms of real data are Hermitian, only the 
; non-negative frequencies need to be dealt with explicitly.  The negative 
; frequency elements may be recovered by using the function CL_CVAL, which 
; returns the complex conjugate for negative frequencies, and zero for 
; frequencies outside the defined range.  However, for essentially all 
; practical purposes, the  negative component can be accounted for by 
; doubling the amplitude of  positive component determined directly from CDFT.
;
; This IDL implementation of the CLEAN algorithm is based on a suite of
; FORTRAN routines developed by Roberts et al.  For details of the 
; concerning algorithm and practical advice regarding its use, please 
; refer to:  Roberts, D.H., Lehar, J., & Dreher, J. W. 1987, AJ, 93, 968
;            "Time Series Analysis with CLEAN. I. Derivation of a Spectrum"
; 
; INPUT PARAMETERS: 
;     FREQ  - frequency vector
;     WFN   - spectral window function                      [complex]
;     DFT   - "dirty" discrete Fourier transform            [complex] 
;     GAIN  - fraction of window function to subtract per iteration 
;             (0 < gain < 2)
;     NCL   - number of CLEAN iterations to perform    
;
; NOTE: The input parameters FREQ, WFN, and DFT will normally be computed
;       by the routine DFOURT.
;
; OUTPUT PARAMETERS:
;     CDFT  - CLEANed discrete Fourier Transform            [complex]
;     CCOMP - individual CLEAN components                   [complex]
;     RDFT  - residual transform after NCL CLEAN iterations [complex]
;
; OTHER ROUTINES CALLED:
;     cl_alpha, cl_beam, cl_subcmp
;
; HISTORY: 
; Jan. 91: translated to IDL                              [A.W.Fullerton, BRI]
;          Source:  FORTRAN code distributed by Roberts et al.
; Apr. 96: reorganized code for more efficient execution  [AWF, USM]
;          replaced slow cl_convolv with much faster intrinsic IDL function
;          reorganized I/O list, improved documentation
; COMMENT: The current implementation of the convolution is very efficient, 
; because it makes use of an intrinsic IDL function.  However, there are 
; small disagreements between this version and the Roberts et al. code for 
; the first MB indices of CDFT.  These difference are due to a minor bug in 
; the Roberts et al. code, which incorrectly includes data with *negative* 
; indices in the convolution summation for indices 0 < i < MB.  The ultimate
; culprit is the routine CVAL, which returns the complex conjugate of 
; C when given a negative index. In this context, "0" should be returned.
;;=============================================================================
; Error checking...
    if (n_params(0) lt 5) then begin
       print,' ' 
       print,' CLEAN, FREQ, WFN, DFT, GAIN, NCL, cdft, ccomp, rdft '
       print,' '
       return
    endif

; Set up arrays
    rdft =dft                             ; first residual DFT = dirty DFT
    ccomp=dft*0.0                         ; array of CLEAN components
    beam=cl_beam(wfn)                     ; compute Gaussan restoring beam

; CLEAN loop: iterate NCL times
    for icl=1,ncl do begin                ; CLEAN the residual spectrum

        pk=where(abs(rdft) eq max(abs(rdft))) ; find current maximum in RDFT
	l=pk(0)                               ; index of current maximum
	cc=gain*cl_alpha(wfn,rdft,l)          ; estimate CLEAN component
        cl_subcmp,wfn,rdft,cc,l               ; shift & scale WFN; subtract CC
	ccomp(l)=ccomp(l) + cc                ; store component in CCOMP
	
    endfor

; Convolve CCOMP with the BEAM to produce the CLEAN FT.  
    mb=(n_elements(beam)-1)/2                         ; elements per half beam
    pad=replicate(complex(0.,0.),mb)                  ; define padding
    input=[pad,ccomp,pad]                             ; pad the data
    cdft=shift(convol(input, beam, center=0), -mb)    ; convolve + recenter
    cdft=cdft(mb:n_elements(input)-mb-1)              ; strip padding  

; Add the final residual vector to preserve noise distribution
    cdft = cdft + rdft                                ; CLEANed DFT
					 
    return                                            ; All done!
    end