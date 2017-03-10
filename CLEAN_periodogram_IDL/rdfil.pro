PRO RDFIL,FILIN,NDIM,buff,text,silent=silent
;
; Utility to read the contents of a ASCII disk file into a 1-D vector or 
; a 2-D array in the IDL environment. The disk file must have the "standard 
; format" (described below).
;
; INPUTS:
;   FIL  - string variable that specifies name of disk file
;   NDIM - integer that specifies dimensionality of data
;          (1--> vector; 2--> array)
;
; OUTPUTS:
;   buff - output structure containing data from file
;            if NDIM=1, buff=DBLARR(NROWS)
;            if NDIM=2, buff=DBLARR(NCOLS,NROWS)
;
;   text - optional strarr containing the header line and any text that 
;          trails the data in the file [maximum of 100 lines]
;
; RESTRICTIONS:
;   The disk file must have the following "standard" structure:
;   line 1: Title line (up to 80 characters)
;   line 2: Dimensions: NROW  --> if NDIM=1 (i.e., vector)
;                             OR
;                       NCOL, NROW  --> if NDIM=2 (i.e., array)
;   lines 3 -> NROW+3: Data:   D_i or D_ij             
;   lines NROW+4 on  : text comments 
;
; MODIFICATION HISTORY:
;   1. Jan. 1987: Written by A. W. Fullerton          [DDO, U. Toronto] 
;   2. June 1990: tested under IDL V2.0 (VAX VMS)     [AWF, Bartol]
;   3. Sept.1991: added optional output of title      [AWF, Bartol]
;   4. Apr. 1996: added ability to read trailing text [AWF, USM]
;============================================================================
    line=''                          ; dummy variable for reading text lines
    text=strarr(100)                 ; text array
    il=0                             ; counter for no. entries in text
    ncol=0L
    nrow=0L

    if (n_params(0) lt 3) then begin
        print,' '
	print,' RDFIL, FILIN, NDIM, buff, [text], /silent '
        print,' '
	return
    endif
    ndim=fix(ndim)

    get_lun,unit
    openr,unit,filin
    readf,unit,line
    text(il)=line

    if not keyword_set(silent) then begin
       print,line
       print,' '
    endif
    
    case ndim of
            1: begin
	         readf,unit,nrow
		 buff=dblarr(nrow)
		 readf,unit,buff
	       end
	    2: begin
	         readf,unit,ncol,nrow
		 buff=dblarr(ncol,nrow)
		 readf,unit,buff
	       end
	 else: begin
	         print,' '
	         print,' ** ERROR in RDFIL: ndim must be 1 or 2'
		 print,' '
                 return
	       end
    endcase

    if not eof(unit) then begin
       while not eof(unit) do begin
             il=il+1
             readf,unit,line
	     text(il)=line
       endwhile
    endif
    text=text(0:il)

    close,unit
    free_lun,unit
    
    return
    end

