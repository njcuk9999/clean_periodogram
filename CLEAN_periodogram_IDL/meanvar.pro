PRO MEANVAR, VECTOR, n, mean, var

; Procedure to evaluate the mean and the variance of an input vector.
;
; INPUTS:
;   VECTOR - vector of n elements, of any type except string
; 
; OUTPUTS:
;   n      - the number of elements in the input vector
;   mean   - the mean value of the input vector
;   var    - the variance about the mean of the input vector
;
; OTHER ROUTINES CALLED:
;   None.
;
; HISTORY:
;   June 1986:  Written by A. W. Fullerton [DDO, U. Toronto]
;   April  96:  Added error trapping.  [AWF, USM]
;=========================================================================
    if (n_params(0) lt 1) then begin
       print,' '
       print,' MEANVAR, VECTOR, n, mean, var '
       print,' '
       return
    endif
    
    n = n_elements(vector) 
    
    if (n eq 1) then begin
        mean = total(vector)/n
	var=0.
	print,' '
	print,' ** WARNING - MEANVAR - only 1 element in the vector! '
	print,' '
    endif else begin
        mean = total(vector)/n
        var  = total((vector - mean)^2)/(n-1)
    endelse

    return
    end