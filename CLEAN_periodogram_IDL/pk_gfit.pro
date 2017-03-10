PRO PK_GFIT, X, GPAR, f, pder
;
; This procedure computes the sum of a Gaussian and a polynomial background:
;
;   F(x) = GPAR(0) * exp(-z^2/2.) + GPAR(3) + GPAR(4)*x + ...
;   where 
;   z(x) = [x - GPAR(1)] / GPAR(2)
; This procedure is intended to be used with CURVEFIT.
;
; INPUTS:
;       X = values of independent variable
;    GPAR = parameters of the function
;            GPAR(0) = amplitude of Gaussian component
;            GPAR(1) = location of Gaussian component
;            GPAR(2) = sigma of Gaussian component
;            GPAR(3) = mean value of background
;            GPAR(4) = slope of linear component of background
;            GPAR(5) = coefficient for the quadratic term
;            GPAR(6) = coefficient for the cubic term... etc.
; There must be at least 3 elements in GPAR.  The number of additional
; terms determines the order of the polynomial that is fit to the background.
;    
; OUTPUTS:
;       f = value of function for X
;    pder = fltarr(n_elements(x),n_elements(gpar)) - optional output array 
;           containing the partial derivative required for LSq fitting with 
;           CURVEFIT.  pder(i,j) = deriviative at ith point wrt jth parameter.
;
; OTHER ROUTINES CALLLED:
;  None.
;
; HISTORY:
; Apr. '96: Modified from GAUSS_FUNCT in IDL User's Library [AWF, USM]
;=============================================================================
    on_error,2                     ; return to caller if an error occurs

    if (n_params(0) lt 2) then begin
       print,' '
       print,' PK_GFIT, X, GPAR, f, [pder]
       print,' '
       return
    endif

    if (gpar(2) eq 0.) then begin
        print,' ** ERROR in PK_GFIT: gpar(2) = sigma is zero '
	print,' ** RETURNING -->'
	print,' '
	return
    endif

    npar=n_elements(gpar)

; Compute function
    z = (x-gpar(1))/gpar(2)
    ez=exp(-z^2/2.)*(abs(z) le 7.)     ; Gaussian part, ignore small terms
    f = gpar(0)*ez 

; Compute partial derivatives
    pder=fltarr(n_elements(x), npar)
    pder(0,0)=ez
    pder(0,1)=gpar(0)*ez*z/gpar(2)
    pder(0,2)=pder(*,1)*z

; Add polynomial part
    if (npar gt 3) then begin 
        fp=poly(x, gpar(3:npar-1))
        f= f + fp 
	for i=3,npar-1 do pder(*,i) = x^(i-3)
    endif

    return
    end
