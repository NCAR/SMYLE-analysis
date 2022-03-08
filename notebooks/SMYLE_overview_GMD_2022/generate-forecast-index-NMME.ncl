;---generate an index (e.g.,Nino-3.4) for NMME forecasts

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl" 

begin

  dir = "/glade/work/xianwu/analysis_output/smyle_output/output/"
  
  ystr = 1982
  yend = 2016
  ny = yend-ystr+1

  print(ny)
  var = "sst"
  ground = "ocean"
  indvar = "Nino-34"
    
  if(indvar .eq. "Nino-34")then
    latmin = -5.
    latmax = 5.
    lonmin = 190.
    lonmax = 240.
  end if
  
;---NMME Forecasts

  modelname = "CMC1-CanCM3"
;  modelname = "CMC2-CanCM4"
;  modelname = "COLA-RSMAS-CCSM4"
;  modelname = "GFDL-CM2p1-aer04"
;  modelname = "GFDL-CM2p5-FLOR-A06"
;  modelname = "GFDL-CM2p5-FLOR-B01"
;  modelname = "NASA-GMAO-062012"
;  modelname = "NCEP-CFSv2"

  f= addfile("/glade/work/xianwu/data/model/nmme/"+var+"-"+modelname+".nc","r")
  print(f)
  x = f->sst(S|:,M|:,L|:,{Y|-5.:5.},{X|190.:240.})    ;196001(0)   198001(240)   198201(264)
  printVarSummary(x)
  n = wgt_areaave_Wrap( x, 1., 1., 0)
  printVarSummary(n)
  init = floattointeger(n&S(0))
  print(init)
  n := n(264-init:264-init+ny*12-1,:,:)
  printVarSummary(n)
  
  dM = dimsizes(x&M)
  dL = dimsizes(x&L)
  print(dM)
  print(dL)
  n := reshape(n(S|:,L|:,M|:),(/ny,12,dL,dM/))    ;12,35,12,10
  printVarSummary(n)

  if( modelname .eq. "NCEP-CFSv2" .or. modelname .eq. "COLA-RSMAS-CCSM4")then    ;climatology definition according to Barnston et al. 2019
    nclim1 = dim_avg_n_Wrap(dim_avg_n_Wrap(n(1982-1982:1998-1982,:,:,:),0),2)
    nclim2 = dim_avg_n_Wrap(dim_avg_n_Wrap(n(1999-1982:2016-1982,:,:,:),0),2)
    printVarSummary(nclim1)
    n(1982-1982:1998-1982,:,:,:) = n(1982-1982:1998-1982,:,:,:) - conform(n(1982-1982:1998-1982,:,:,:),nclim1,(/1,2/))
    n(1999-1982:2016-1982,:,:,:) = n(1999-1982:2016-1982,:,:,:) - conform(n(1999-1982:2016-1982,:,:,:),nclim2,(/1,2/))
    printVarSummary(n)
  else 
    nclim = dim_avg_n_Wrap(dim_avg_n_Wrap(n(:,:,:,:),0),2)
    printVarSummary(nclim)
    n = n - conform(n,nclim,(/1,2/))
    printVarSummary(n)
  end if

  
  n!0 = "S"
  n&S = fspan(ystr,yend,ny)
  n!1 = "initmonth"
  n&initmonth = ispan(1,12,1)
  n!2 = "L"
  n&L = x&L
  n!3 = "M"
  n&M = x&M

  if (any(isnan_ieee(n))) then
      value = 1.e30
      replace_ieeenan (n, value, 0)
      n@_FillValue = value
  end if
  
  printVarSummary(n)
  
  system("rm -f " + dir + get_script_prefix_name() + "-"+ indvar+"-"+modelname+".nc")
  f = addfile(dir + get_script_prefix_name() + "-"+ indvar+"-"+modelname+".nc", "c")
  f->n = n
end

  

