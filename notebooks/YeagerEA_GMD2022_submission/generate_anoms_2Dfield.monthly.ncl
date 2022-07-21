;===============================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
; ================================================;
begin

;; Primary settings
  writeanom = True
  writedrift = False

  dir = "/glade/work/xianwu/analysis_output/smyle_output/postprocess/"

  timeavg = "mon"
  
  season = "Nov"
;  season = "Feb"
;  season = "May"
;  season = "Aug"
  
  styr1 =  1970              ; first forecast year 1
  styr2 =  2019               ; last forecast year 1
  nlag  = 24
  nmem  = 20
  
  if(season .eq. "Nov")then
    monnum = 11
  end if
  
  if(season .eq. "Feb")then
    monnum = 2
  end if
  
  if(season .eq. "May")then
    monnum = 5
  end if

  if(season .eq. "Aug")then
    monnum = 8
  end if

  var = "TS"
;  var = "Z200"
;  var = "UBOT"
;  var = "VBOT"
;  var = "FSNS"
;  var = "FLNS"
;  var = "LHFLX"
;  var = "SHFLX"
;  var = "PSL"
;  var = "PREC"
;  var = "TREFHT"
  



  fin = dir+"SMYLE."+season+"."+var+"."+timeavg+"mean"

  fout = fin + ".anom.nc"
  fin = fin + ".nc"

  f       = addfile(fin,"r")
  dpdata    = f->$var$(S|:,L|0:nlag-1,M|:,lat|:,lon|:) 
  delete(f)
  clim = dpdata(0,:,0,:,:)
  istyr = 0             ; keep track of ensembles to include in climatology
  do ilead=0,nlag-1 ; loop over lead months
    imon = mod((ilead + monnum-1),12)         ; month of year (0-11)
    if (imon.eq.0) then
      istyr = istyr+1
    end if
    st0 = 2-istyr       ; start with "1972-mon" ensemble when istyr=0
    st1 = st0+47        ; end with "2015-mon" ensemble when istyr=0
    print("ilead = "+ilead+". average over ensembles: "+(1970+st0)+"0"+monnum+" to "+(1970+st1)+"0"+monnum)
    clim(ilead,:,:) = (/ dim_avg_n_Wrap(dim_avg_n_Wrap(dpdata(st0:st1,ilead,:,:,:),0),0) /)
  end do

;;  Convert DP fields into anomalies from the above climatology
  dpdata_an = dpdata
  printVarSummary(dpdata_an)

  tmp = conform(dpdata,clim,(/1,3,4/))
  dpdata_an = (/ dpdata - tmp /)

; ============================================================================
; Output
; ============================================================================
  system("\rm -f "+fout)
  setfileoption("nc","Format","LargeFile")
  f = addfile (fout, "c")
    setfileoption(f,"DefineMode",True)
    gAtt = True
    datestamp = systemfunc("date")
    gAtt@history     = "created by X. Wu on "+datestamp
    gAtt@script      = "/glade/u/home/xianwu/ncl/CESM2-SMYLE/postprocess/generate_anoms_2Dfield.monthly.ncl"
    setfileoption(f,"DefineMode",False)
    fileattdef(f,gAtt)
  
  f@climatology = "1972-2019, computed separately for each lag"
  printVarSummary(dpdata_an)
  printVarSummary(clim)
  f->anom=dpdata_an
;  f->clim=clim
end
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
