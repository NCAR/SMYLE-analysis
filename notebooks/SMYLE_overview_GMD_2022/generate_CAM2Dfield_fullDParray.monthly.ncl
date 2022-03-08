;===============================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
; ================================================;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Purpose:  Generate a complete DP array file, dimensioned
;; (Start year,Lag,Member,nlon,nlat), for any given 2D atmospheric field  
;; (e.g., TS).  Output is raw (not bias-corrected).
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
begin

;; Primary settings
  outdir = "/glade/work/xianwu/analysis_output/smyle_output/postprocess/"
  smyledir = "/glade/campaign/cesm/development/espwg/SMYLE/archive/"

;  season = "Nov"
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
  
  nlat = 192
  nlon = 288 
  
  varstr   = "TS"
;  varstr   = "PSL"
;  varstr   = "UBOT"
;  varstr   = "VBOT"
;  varstr = "TREFHT"
;  varstr = "FSNS"
;  varstr = "FLNS"
;  varstr = "LHFLX"
;  varstr = "SHFLX"  
;   varstr = "PREC"

  basevarstr = varstr
  if (.not.ismissing(str_index_of_substr(varstr,"PD",-1))) then
    basevarstr = "PD"
  end if
  
  if (.not.ismissing(str_index_of_substr(varstr,"UBOT",-1))) then
    basevarstr = "U"
  end if

  if (.not.ismissing(str_index_of_substr(varstr,"VBOT",-1))) then
    basevarstr = "V"
  end if

  if(.not.ismissing(str_index_of_substr(varstr,"PREC",-1))) then
    basevarstr = "PRECC"
    indir = smyledir +"PRECC"+"/"
    basevarstr2 = "PRECL"
    indir2 = smyledir +"PRECL"+"/"
  end if
  
  fout = outdir + "SMYLE."+season+"."+varstr+".monmean.nc"
  
;; Start work
  nstyr = styr2 - styr1 + 1
  startyear = ispan(styr1,styr2,1)
    startyear!0 = "S"
    startyear@long_name = "ensemble"
    startyear@description = "historical year corresponding to forecast year 1"
    startyear@example = "S=1970 for forecasts initialized on "+season+" 1 1970"
  lag = ispan(1,nlag,1)
    lag!0 = "L"
    lag@long_name = "forecast month"
  mem = ispan(1,nmem,1)
    mem!0 = "M"
    mem@long_name = "ensemble member"
  field = new((/nstyr,nlag,nmem,nlat,nlon/),"float",1.e30)
     field!0 = "S"
     field&S = startyear
     field!1 = "L"
     field&L = lag
     field!2 = "M"
     field&M = mem
     field!3 = "lat"
     field!4 = "lon"

  do iyr = 0,nstyr-1
      styr1str = sprinti("%0.4i", iyr+styr1)
      styr2str = sprinti("%0.4i", iyr+styr1+2)
      monnum2 = monnum-1
      monnumstr = sprinti("%0.2i", monnum)
      monnumstr2 = sprinti("%0.2i", monnum2)
    do im = 0,nmem-1
      memstr = sprinti("%0.3i", im+1)
      run_name= "b.e21.BSMYLE.f09_g17." + styr1str + "-"+monnumstr+"." +memstr
      indir = smyledir+run_name
      f = addfile(indir+"/atm/proc/tseries/month_1/"+run_name+".cam.h0."+basevarstr+"."+styr1str+monnumstr+"-"+styr2str+monnumstr2+".nc", "r")
      
      print("processing : "+run_name)  

      if (varstr.eq."PSL") then
         data = f->PSL(0:nlag-1,:,:)
      end if
      if (varstr.eq."PRECC") then
         data = f->PRECC(0:nlag-1,:,:)
      end if
      if (varstr.eq."PRECL") then
         data = f->PRECL(0:nlag-1,:,:)
      end if
      if (varstr.eq."TS") then
         data = f->TS(0:nlag-1,:,:)
      end if
      if (varstr.eq."TREFHT") then
         data = f->TREFHT(0:nlag-1,:,:)
      end if
      if (varstr.eq."UBOT")then
         data = f->U(0:nlag-1,29,:,:)
      end if
      if (varstr.eq."VBOT")then
         data = f->V(0:nlag-1,29,:,:)
      end if
      if (varstr.eq."FSNS")then
         data = f->FSNS(0:nlag-1,:,:)
      end if
      if (varstr.eq."FLNS")then
         data = f->FLNS(0:nlag-1,:,:)
      end if 
      if (varstr.eq."LHFLX")then
         data = f->LHFLX(0:nlag-1,:,:)
      end if  
      if (varstr.eq."SHFLX")then
         data = f->SHFLX(0:nlag-1,:,:)
      end if        
      if (varstr.eq."PREC") then
        data = f->PRECC(0:nlag-1,:,:)
        f2 = addfile(indir+"/atm/proc/tseries/month_1/"+run_name+".cam.h0."+basevarstr2+"."+styr1str+monnumstr+"-"+styr2str+monnumstr2+".nc", "r")
        data = (/ data + f2->PRECL(0:nlag-1,:,:) /)
      end if   
    
      field(iyr,:,im,:,:) = (/ data /)
    end do
  end do

  field@units = data@units
  field@long_name = data@long_name

  printVarSummary(field)
  
  system("rm -f "+fout)
  setfileoption("nc","Format","LargeFile")
  f = addfile(fout,"c")
  setfileoption(f,"DefineMode",True)
       dimNames = (/"S","L","M","lat","lon"/)
       dimSizes = (/nstyr,nlag,nmem,nlat,nlon/)  
       dimUnlim = (/False,False,False,False,False /)
     filedimdef(f,dimNames,dimSizes,dimUnlim)  	; hard-coded dims
       gAtt = True
       datestamp = systemfunc("date")
       gAtt@history     = "created by X. Wu on "+datestamp
       gAtt@script      = "/glade/u/home/xianwu/ncl/CESM2-SMYLE/postprocess/generate_CAM2Dfield_fullDParray.monthly.ncl"
     field&lon = lonGlobeF (nlon, "lon", "longitude", "degrees_east")
     field&lat = latGlobeF (nlat, "lat", "latitude", "degrees_north")
     
     fileattdef(f,gAtt)           			; add global attributes
     filevardef(f,"S" ,typeof(field&S),getvardims(field&S))
     filevardef(f,"L" ,typeof(field&L),getvardims(field&L))
     filevardef(f,"M" ,typeof(field&M),getvardims(field&M))
     filevardef(f,varstr,typeof(field),getvardims(field))
     filevarattdef(f,"S",field&S)
     filevarattdef(f,"L",field&L)
     filevarattdef(f,"M",field&M)
     filevarattdef(f,varstr,field)
  setfileoption(f,"DefineMode",False)
  f->S= (/ field&S /)
  f->L= (/ field&L /)
  f->M= (/ field&M /)
  f->$varstr$= (/ field /)

end
