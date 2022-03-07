;---Make plot for timeseries of an index (e.g., Nino-3.4) for a selected season

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl" 

begin
  diro = "/glade/work/xianwu/data/obs/hadisst/"
  dirf = "/glade/work/xianwu/analysis_output/smyle_output/postprocess/"
  dirom = "/glade/work/xianwu/data/obs/mergedata/"
 
  ystro = 1970
  yendo = 2021
  nyo = yendo-ystro+1
 
  ystr = 1970
  yend = 2021
  ny = yend - ystr +1
   
  detrend = "linear"

  season = "DJF"

  var = "TS"
  ground = "ocean"

  indvar = "Nino-34"

  if(indvar .eq. "Nino-34")then
    latmin = -5.
    latmax = 5.
    lonmin = 190.
    lonmax = 240.
  end if

  fmonmax = 45

;---Observation
  if(var .eq. "TS")then
    f = addfile(dirom + "merge-ts-sst-obsdata-BEST-HadISST-1900-2021.nc", "r")
    if(indvar .eq. "Nino-34")then
      t = f->ts(time|(ystro-1900)*12:(yendo-1900+1)*12-1,{lat|latmin:latmax},{lon|lonmin:lonmax}) 
      n = wgt_areaave_Wrap( t, 1., 1., 0)
      delete( t)
    end if
  end if
  
  n = rmAnnCycle1D( n)
  
  if(detrend .eq. "quadratic")then
    n = (/ dtrend_quadratic_msg_n( n, False, False, 0) /)
  end if
  
  if(detrend .eq. "linear")then
    n = (/ dtrend_msg_n( n&time, n, False, False, 0) /)
  end if
   
  n1 = n
  n2 = n  
  
  n1 = runave(n, 3, 0)   ;three-month average for DJF 
  n2 = runave(n, 2, 0)   ;two-month average for DJ
  
  n = n1
  
  printVarSummary(n)
  printVarSummary(n2)

  if(season .eq. "DJF")then
    nosea = n(12::12)     ;1970DJF-
    nosea2 = n2(12::12)   ;1970DJ
    nosea&time = fspan(ystro,yendo-1,nyo-1) 
    nosea2&time = fspan(ystro,yendo-1,nyo-1) 
    nosea := array_append_record(nosea,nosea(0:2),0)
    nosea2 := array_append_record(nosea2,nosea2(0:2),0)
    printVarSummary(nosea)
    print(nyo)
    nosea(nyo-1:) = nosea@_FillValue
    nosea2(nyo-1:) = nosea2@_FillValue
  end if
  
  nosea&time = ispan(ystro,yendo+2,1)
  printVarSummary(nosea)

  delete(n)
  
;---Forecasts    
  startmon =(/"Feb","May","Aug","Nov"/)
  fmon = 24
  n = new((/4,ny,fmon,20/),"float")

  do i=0,3 
    f = addfile( dirf + "SMYLE."+startmon(i)+"."+var+".monmean.anom.nc", "r")
    x = f->anom(S|ystr-1970:2019-1970,L|0:fmon-1,M|:,lat|:,lon|:)
    x&lat = latGlobeF (192, "lat", "latitude", "degrees_north")
    x&lon = lonGlobeF (288, "lon", "longitude", "degrees_east")   
    if(var .eq. "TS" .or. var .eq. "PREC")then  
      if(ground .eq. "ocean")then
        f = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
        LSMASK = f->LSMASK 
        lsmask = landsea_mask( LSMASK, x&lat, x&lon)
        x(S|:,L|:,M|:,lat|:,lon|:) = mask( x(S|:,L|:,M|:,lat|:,lon|:), lsmask.le.0.5, True)  ;only take ocean value
        delete( lsmask)    
      end if
      if(ground .eq. "land")then
        f = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
        LSMASK = f->LSMASK 
        lsmask = landsea_mask( LSMASK, x&lat, x&lon)
        x(S|:,L|:,M|:,lat|:,lon|:) = mask( x(S|:,L|:,M|:,lat|:,lon|:), lsmask.eq.1, True)  ;only take land value
        delete( lsmask)    
      end if    
    end if
    
    if(i .eq. 0 .or. i .eq. 1)then
      f = addfile( dirf + "SMYLE."+startmon(i)+"."+var+".monmean-realtime.anom.nc", "r")
      xrealtime = f->anom(S|2020-2020:2021-2020,L|0:fmon-1,M|:,lat|:,lon|:)
    end if
    
    if(i .eq. 2 .or. i .eq. 3)then
      f = addfile( dirf + "SMYLE."+startmon(i)+"."+var+".monmean-realtime.anom.nc", "r")
      xrealtime = f->anom(S|2020-2020:2020-2020,L|0:fmon-1,M|:,lat|:,lon|:)
    end if
    xrealtime&lat = latGlobeF (192, "lat", "latitude", "degrees_north")
    xrealtime&lon = lonGlobeF (288, "lon", "longitude", "degrees_east")   
    if(var .eq. "TS" .or. var .eq. "PREC")then  
      if(ground .eq. "ocean")then
        f = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
        LSMASK = f->LSMASK 
        lsmask = landsea_mask( LSMASK, x&lat, x&lon)
        xrealtime(S|:,L|:,M|:,lat|:,lon|:) = mask( xrealtime(S|:,L|:,M|:,lat|:,lon|:), lsmask.le.0.5, True)  ;only take ocean value
        delete( lsmask)    
      end if
      if(ground .eq. "land")then
        f = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
        LSMASK = f->LSMASK 
        lsmask = landsea_mask( LSMASK, x&lat, x&lon)
        xrealtime(S|:,L|:,M|:,lat|:,lon|:) = mask( xrealtime(S|:,L|:,M|:,lat|:,lon|:), lsmask.eq.1, True)  ;only take land value
        delete( lsmask)    
      end if    
    end if    
    
    x := x(S|:,L|:,M|:,{lat|latmin:latmax},{lon|lonmin:lonmax})
    xrealtime := xrealtime(S|:,L|:,M|:,{lat|latmin:latmax},{lon|lonmin:lonmax})
    n(i,0:ny-1-2,:,:) = wgt_areaave_Wrap( x, 1., 1., 0)  
    if(i .eq. 0 .or. i .eq. 1)then
      n(i,ny-1-1:ny-1,:,:) = wgt_areaave_Wrap( xrealtime, 1., 1., 0)  
    end if
    if(i .eq. 2 .or. i .eq. 3)then
      n(i,ny-1-1:ny-1-1,:,:) = wgt_areaave_Wrap( xrealtime, 1., 1., 0)  
    end if     
    delete(x)
    delete(xrealtime)
  end do
  
  n!0 = "startmon"
  n!1 = "S"
  n!2 = "L"
  n!3 = "M"
  n&S = ispan(ystr,yend,1)
  printVarSummary(n) 
  
  nm = dim_avg_n_Wrap(n,3)
  
  printVarSummary(nm)

  if(detrend .eq. "linear")then
    nm_dt = dtrend_msg_n(nm&S, nm, True, False, 1)
    nm_t = nm-nm_dt  ;trend based on ensemble mean
  end if

  ntemp = n-conform(n,nm_t,(/0,1,2/))  
  n = runave_n_Wrap(ntemp,3,0,2)  ;detrend all ensemble members


  if(season .eq. "DJF")then
    n(0,:,23,:) = (ntemp(0,:,22,:)+ntemp(0,:,23,:))/2.  ;DJ for 22-month lead
    n_new = new((/8,ny+1,20/),"float")  ;1970 DJF-
    
    printVarSummary(n)
    printVarSummary(n_new)
    print(ny)
    n_new(0,1:ny,:) = n(0,:,23,:) ;Feb1970(0) - 22, D(22)J(23)F1971
    n_new(1,1:ny,:) = n(1,:,20,:) ;May1970
    n_new(2,1:ny,:) = n(2,:,17,:) ;Aug1970
    n_new(3,1:ny,:) = n(3,:,14,:) ;Nov1970
    n_new(4,:ny-1,:) = n(0,:,11,:) ;Feb1970 - 10, D(10)J(11)F1971
    n_new(5,:ny-1,:) = n(1,:,8,:) ;May1970
    n_new(6,:ny-1,:) = n(2,:,5,:) ;Aug1970
    n_new(7,:ny-1,:) = n(3,:,2,:) ;Nov1970
  end if
  
  nfseam = dim_avg_n_Wrap(n_new,2)
  nfseas = dim_stddev_n_Wrap(n_new,2)
  
  printVarSummary(nosea)
  printVarSummary(nfseam)
  
;---Figures
  wks  = gsn_open_wks( "ps",  get_script_prefix_name()+"-"+indvar)  
  
;---Time Series
  res = True
  res@gsnDraw = False
  res@gsnFrame = False
  res@gsnRightString = ""
  res@gsnStringFontHeightF = 0.015
  res@gsnLeftString = ""  

  res@trXMinF = 1969.5
  res@trXMaxF = 2022.5
  res@gsnXRefLineColor = "Gray50"
  res@gsnXRefLineDashPattern = 2
  res@gsnYRefLineColor = "Gray50"
  res@gsnYRefLineDashPattern = 2
  res@gsnYRefLine = 0
  res@tmXTOn = False
  res@tmYROn = False
  res@tmXBLabelFontHeightF = 0.012
  res@tmYLLabelFontHeightF = 0.012
  res@tmXBMinorThicknessF = 1.
  res@tmXBLabelJust = "CenterCenter"
  res@tiXAxisString = ""
  res@tiYAxisString = ""
  res@tiYAxisFontHeightF = 0.014
  res@xyDashPattern = "0"
  res@xyLineThicknessF = "2."  
  res@xyMarkLineMode = "MarkLines"
  res@xyMarker = 16
  res@xyMarkerSizeF = 0.005
  
;---res2 
  res2 = res
  res2@xyLineColor = "transparent"
  delete(res2@gsnXRefLine)
  delete(res2@xyMarkLineMode)

;---Text
  restext = True
  
  restext@txFontHeightF =  0.012
  
  col = (/"red","blue","orange","forest green"/)
  
  plot = new(4, "graphic") 
  plot2 = new(4, "graphic") 

  if(indvar .eq. "Nino-34")then
    res@trYMinF = -3.5
    res@trYMaxF = 3.5
    res@tmYLPrecision = 1
  end if
  
  res@vpWidthF = 0.4
  res@vpHeightF = 0.13  
  
  restext@txJust = "TopCenter"
  ntilde  = "n~H-13V2F35~D~FV-2H3~"
  Nino = "Ni" + ntilde + "o-3.4 "
  gsn_text_ndc( wks, season+" "+Nino +"Index", 0.5, 0.99, restext)

  if(season .eq. "DJF")then
    forecastname = (/"Feb Forecasts (22-mon lead)","May Forecasts (19-mon lead)","Aug Forecasts (16-mon lead)","Nov Forecasts (13-mon lead)",\
                   "Feb Forecasts (10-mon lead)","May Forecasts (7-mon lead)","Aug Forecasts (4-mon lead)","Nov Forecasts (1-mon lead)"/)
  end if

  restext@txJust = "TopLeft"  
  res@vpXF = 0.53
  do i=0,3
    res@vpYF = 0.38+(res@vpHeightF+0.06)*i
    if(i .eq. 0)then
      corr = escorc(nosea2(0:ny),nfseam(i,:))
      rmse = dim_rmsd(nosea2(0:ny),nfseam(i,:))/ dim_stddev_n_Wrap(nosea2(0:ny), 0)
    else
      corr = escorc(nosea(0:ny),nfseam(i,:))
      rmse = dim_rmsd(nosea(0:ny),nfseam(i,:))/ dim_stddev_n_Wrap(nosea(0:ny), 0)
    end if
    gsn_text_ndc( wks, forecastname(i)+", ACC="+sprintf("%3.2f", corr)+", nRMSE="+sprintf("%3.2f", rmse), res@vpXF, res@vpYF+0.02, restext)
    res@xyLineColor = "black"
    res@xyMarkLineMode = "MarkLines"
    res@xyMarkerSizeF = 0.003
    if(i .eq. 0)then
      plot(i) = gsn_csm_xy( wks, nosea2&time, nosea2, res)     ;DJ
    else
      plot(i) = gsn_csm_xy( wks, nosea&time, nosea, res)       ;DJF
    end if  
    res@xyLineColor = col(i)
    delete(res@xyMarkLineMode)
    plotf = gsn_csm_xy( wks, nosea&time(0:ny), nfseam(i,:), res)     
    overlay(plot(i), plotf)
    
    res2@xyLineColor = col(i)
    res2@xyLineThicknessF = "0"
    res2@gsnXYFillColors = col(i)
    res2@gsnXYFillOpacities = 0.1
    plots = gsn_csm_xy(wks, nosea&time(0:ny), (/nfseam(i,:)+nfseas(i,:), nfseam(i,:)-nfseas(i,:)/), res2)
    overlay(plot(i), plots)   
  end do
  
  res@vpXF = 0.07
  do i=0,3
    res@vpYF = 0.38+(res@vpHeightF+0.06)*i
    corr = escorc(nosea(0:ny),nfseam(i+4,:))
    rmse = dim_rmsd(nosea(0:ny),nfseam(i+4,:))/ dim_stddev_n_Wrap(nosea(0:ny), 0)
    gsn_text_ndc( wks, forecastname(i+4)+", ACC="+sprintf("%3.2f", corr)+", nRMSE="+sprintf("%3.2f", rmse), res@vpXF, res@vpYF+0.02, restext)
    res@xyLineColor = "black"
    res@xyMarkLineMode = "MarkLines"
    res@xyMarkerSizeF = 0.003
    plot2(i) = gsn_csm_xy( wks, nosea&time, nosea, res)     
    res@xyLineColor = col(i)
    delete(res@xyMarkLineMode)
    plotf = gsn_csm_xy( wks, nosea&time(0:ny), nfseam(i+4,:), res)     
    overlay(plot2(i), plotf)
    
    res2@xyLineColor = col(i)
    res2@xyLineThicknessF = "0"
    res2@gsnXYFillColors = col(i)
    res2@gsnXYFillOpacities = 0.1
    plots = gsn_csm_xy(wks, nosea&time(0:ny), (/nfseam(i+4,:)+nfseas(i+4,:), nfseam(i+4,:)-nfseas(i+4,:)/), res2)
    overlay(plot2(i), plots) 
  end do
  

;  drawNDCGrid( wks)
  draw( wks)
  frame( wks)
end
