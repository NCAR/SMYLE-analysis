;---ACC and RMSE for any index

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin
  diro = "/glade/work/xianwu/data/obs/"
  dirf = "/glade/work/xianwu/analysis_output/smyle_output/postprocess/"
  dirfo = "/glade/work/xianwu/analysis_output/smyle_output/output/"

  figure = "pdf"
  figure = "ps"
  
  modelname = (/"CMC1-CanCM3", "CMC2-CanCM4", "COLA-RSMAS-CCSM4","GFDL-CM2p1-aer04",\
                "GFDL-CM2p5-FLOR-A06","GFDL-CM2p5-FLOR-B01","NASA-GMAO-062012","NCEP-CFSv2"/)

  
  ystro = 1970
  yendo = 2021
  
  fmon = 24
  
  ystr = 1970
  yend = 2019
  ny = yend-ystr+1
  
;  detrend = "quadratic"
;  detrend = "nondetrend"
  detrend = "linear"
  
  define = "12yr"  

  var = "TS"
  indvar = "Nino-34"
  ground = "ocean"

  if(indvar .eq. "Nino-34")
    latmin = -5.
    latmax = 5.
    lonmin = 190.
    lonmax = 240.
  end if
  
;---Observation
  if(var .eq. "TS")then
    if(indvar .eq. "Nino-34")then
      f = addfile( diro + "hadisst/HadISST_sst.nc", "r")
      t = f->sst(time|(ystro-1870)*12:(yendo-1870+1)*12-1,{latitude|-5.:5.},{longitude|-170.:-120.})
      n = wgt_areaave_Wrap( t, 1., 1., 0)
      delete( t)
    end if
  end if
  printVarSummary(n)
  n = rmAnnCycle1D( n)
  
  if(detrend .eq. "quadratic")then
    n = (/ dtrend_quadratic_msg_n( n, False, False, 0) /)
  end if
  
  if(detrend .eq. "linear")then
    n = (/ dtrend_msg_n( n&time, n, False, False, 0) /)
  end if
    
  printVarSummary(n)
  no = new((/4,ny,fmon/),"float")
  do i=0,ny-1
    no(0,i,:) = n(time|i*12+1:i*12+1+fmon-1)   ;1970 Feb - 
    no(1,i,:) = n(time|i*12+4:i*12+4+fmon-1)   ;1970 May - 
    no(2,i,:) = n(time|i*12+7:i*12+7+fmon-1)   ;1970 Aug - 
    no(3,i,:) = n(time|i*12+10:i*12+10+fmon-1)   ;1970 Nov - 
  end do  
 
  delete(n)
  
;---SMYLE Forecasts    
  startmon =(/"Feb","May","Aug","Nov"/)
  fmon = 24
  n = new((/4,ny,fmon,20/),"float")

  do i=0,3 
    f = addfile( dirf + "SMYLE."+startmon(i)+"."+var+".monmean.anom.nc", "r")
    x = f->anom(S|ystr-1970:yend-1970,L|0:fmon-1,M|:,lat|:,lon|:)
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
    
    x := x(S|:,L|:,M|:,{lat|latmin:latmax},{lon|lonmin:lonmax})
    n(i,0:ny-1,:,:) = wgt_areaave_Wrap( x, 1., 1., 0)  
    delete(x)
  end do
  
  n!0 = "startmon"
  n!1 = "S"
  n!2 = "L"
  n!3 = "M"
  n&S = ispan(ystr,yend,1)
  printVarSummary(n) 
  
  nm = dim_avg_n_Wrap(n,3)

;---NMME forecasts
  nnmme = new((/dimsizes(modelname)+3,ny,4,fmon/),"float")
  do i=0,dimsizes(modelname)-1
    f = addfile( dirfo + "/generate-forecast-index-NMME-"+indvar+"-"+modelname(i)+".nc", "r")
    ntemp = f->n
    dtemp = dimsizes(ntemp)
    nnmme(i,1982-ystr:2016-ystr,0,:dtemp(2)-1) = dim_avg_n_Wrap(ntemp(:,1,:,:),2)  ;select 4 seasons
    nnmme(i,1982-ystr:2016-ystr,1,:dtemp(2)-1) = dim_avg_n_Wrap(ntemp(:,4,:,:),2)  ;select 4 seasons
    nnmme(i,1982-ystr:2016-ystr,2,:dtemp(2)-1) = dim_avg_n_Wrap(ntemp(:,7,:,:),2)  ;select 4 seasons
    nnmme(i,1982-ystr:2016-ystr,3,:dtemp(2)-1) = dim_avg_n_Wrap(ntemp(:,10,:,:),2)  ;select 4 seasons
 
    delete(ntemp)
    delete(dtemp)
  end do

  printVarSummary(nnmme)

;---Combine SMYLE, NMME individual models, NMME MME together
  nnmme(dimsizes(modelname),:,:,:) = dim_avg_n_Wrap(nnmme(0:dimsizes(modelname)-1,:,:,:),0) 
  nnmme(dimsizes(modelname)+1,:,:,:) = nm(S|:,startmon|:,L|:)   ;1970-2019
  nnmme(dimsizes(modelname)+2,1982-ystr:2016-ystr,:,:) = nm(S|1982-ystr:2016-ystr,startmon|:,L|:)  ;1982-2016

  if(detrend .eq. "linear")then
    nnmme = dtrend_msg_n(nnmme&S, nnmme, True, False, 1)
  end if

;---3-month running mean for seasonal values
  no3 = runave_n_Wrap(no, 3, 0, 2)
  nnmme3 = runave_n_Wrap(nnmme, 3, 0, 3)

;---Correlation SKill
   
  r = nnmme(:,0:2,:,:)
  r_sig = nnmme(:,0:2,:,:)
  rp = r
  rp_sig = r_sig
  
  r3 = nnmme(:,0:2,:,:)
  r_sig3 = nnmme(:,0:2,:,:)
  rp3 = r3
  rp_sig3 = r_sig3


  do k=0,dimsizes(modelname)-1+3
    do i=0,fmon-1
      do j=0,3
      r(k,0,j,i) = escorc_n(no(j,:,i), nnmme(k,:,j,i), 0, 0)  ;All years
      r3(k,0,j,i) = escorc_n(no3(j,:,i), nnmme3(k,:,j,i), 0, 0)  ;All years
      Nr0 = ny
      prob0 = rtest(r(k,0,j,i),Nr0,0)
      r_sig(k,0,j,i) = mask(r(k,0,j,i), prob0 .gt. 0.05 , False)  ;<0.05 significant      

      prob03 = rtest(r3(k,0,j,i),Nr0,0)
      r_sig3(k,0,j,i) = mask(r3(k,0,j,i), prob0 .gt. 0.05 , False)  ;<0.05 significant
      end do
    end do
  end do

  e = nnmme(:,0:2,:,:)
  ep = e
  s = e
  
  e3 = nnmme3(:,0:2,:,:)
  ep3 = e3
  s3 = e3
  
  temps = no
  temps = 0.
  
  do k=0,dimsizes(modelname)-1+3
    do i=0,fmon-1
      do j=0,3
        e(k,0,j,i) = dim_rmsd_n_Wrap(no(j,:,i), nnmme(k,:,j,i), 0)  ;All years
        s(k,0,j,i) = dim_stddev_n_Wrap(no(j,:,i), 0)  ;All years  RMS temps=0 
        
        e3(k,0,j,i) = dim_rmsd_n_Wrap(no3(j,:,i), nnmme3(k,:,j,i), 0)  ;All years
        s3(k,0,j,i) = dim_stddev_n_Wrap(no3(j,:,i), 0)  ;All years  RMS temps=0 
      end do
    end do
  end do
  e = e/s
  e3 = e3/s3
  
  printVarSummary(r)
  printVarSummary(e)

  r&L = ispan(1,fmon,1)
  e&L = ispan(1,fmon,1)

  r3&L = ispan(1,fmon,1)
  e3&L = ispan(1,fmon,1)

  rp&L = ispan(1,fmon,1)
  ep&L = ispan(1,fmon,1)

  s&L = ispan(1,fmon,1)
  s&L = ispan(1,fmon,1)

  printVarSummary(r)
  printVarSummary(e)
  rall = dim_avg_n_Wrap(r, 2)
  eall = dim_avg_n_Wrap(e, 2)
  r_sig_all = dim_avg_n_Wrap(r, 2)
  
  rnew = new((/11,3,5,24/),"float")
  enew = new((/11,3,5,24/),"float")
  r_sig_new = new((/11,3,5,24/),"float")
  
  rnew(:,:,0:3,:) = r3
  enew(:,:,0:3,:) = e3
  r_sig_new(:,:,0:3,:) = r_sig3
  
  
  rnew(:,:,4,:) = rall
  enew(:,:,4,:) = eall
  r_sig_new(:,:,4,:) = r_sig_all
  
  r := rnew
  e := enew
  r_sig := r_sig_new
  
;---Figures
;  
  wks  = gsn_open_wks( figure, get_script_prefix_name() + "-" + indvar + "-"+ detrend) 
  gsn_merge_colormaps( wks, "amwg_blueyellowred",(/"forest green","purple","black"/))
  cmap = gsn_retrieve_colormap( wks)
  mapcolor = ispan( 2, 20, 1)  
;---Time Series
  res = True
  res@gsnDraw = False
  res@gsnFrame = False
  res@gsnRightString = ""
  res@gsnStringFontHeightF = 0.015
  res@gsnLeftString = ""


  res@tmXBMode = "Explicit"
  res@tmXBLabelFontHeightF = 0.01
  res@tmXBValues = ispan( 3, 22, 3) 
;  res@tmXBLabels = (/"1","4","7","10","13", "16","19","22","25", "28","31","34"/)
  res@trXMinF = 0.5
  res@trXMaxF = 24.5
  
  res@tmYLLabelFontHeightF = 0.012
  res@tmYLPrecision = 2
  res@tmYLMode = "Manual"
  res@tmYLTickSpacingF =0.4

  res@gsnXRefLine = ispan( 3, 44, 3)
  res@gsnXRefLineColor = "Gray50"
  res@gsnXRefLineDashPattern = 2
;  res@gsnYRefLine = (/-0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.2,1.4/)
;  res@gsnYRefLineColor = "Gray50"

  res@tmXTOn = False
  res@tmYROn = False
  res@tmXBMinorThicknessF = 2.
  res@tmXBMode = "Explicit" 
  res@tmXBLabelJust = "CenterCenter"
  res@tiXAxisString = ""
  res@tiYAxisString = ""
  res@tiYAxisFontHeightF = 0.014

  res@xyLineThicknessF = "2"
  res@xyMarkLineMode = "MarkLines"
  res@xyMarkerSizeF = 0.006

  
;---res3 
  res3 = res
  delete(res3@gsnXRefLine)
  delete(res3@xyMarkLineMode)


;---Text
  restext = True
  restext@txFontHeightF = 0.012
  

  plot = new(12, "graphic")
  plotmme = new(12, "graphic") 
  res@vpXF = 0.12
  
  res@vpWidthF = 0.35
  res@vpHeightF = 0.14
  
  Ntilde  = "N~H-15V6F35~D~FV-6H3~"
  ntilde  = "n~H-13V2F35~D~FV-2H3~"
  ElNino = "El Ni" + ntilde  +  "o"
  LaNina = "La Ni" + ntilde  +  "a"  
  
  typenum = 0
  type = (/"All years","ENSO growth (observed ENSO peak at Dec~S~1~N~)","ENSO decay (observed ENSO peak at Dec~S~1~N~)"/)
  restext@txJust = "TopCenter"
  gsn_text_ndc( wks, "Correlation Skill", res@vpXF + res@vpWidthF/2., 0.985, restext)
  gsn_text_ndc( wks, "Normalized RMSE", 0.53 + res@vpWidthF/2., 0.985, restext)
  restext@txJust = "TopLeft"

  restext@txFontHeightF = 0.011

;---Correlation
   
  monthtitle_l = (/"(a) Feb Initialized","(b) May Initialized","(c) Aug Initialized","(d) Nov Initialized","(e) All Initialization Average"/)
  monthtitle_r = (/"(f) Feb Initialized","(g) May Initialized","(h) Aug Initialized","(i) Nov Initialized","(j) All Initialization Average"/)


;          labels = (/"1(MAM)","4(JJA)","7(SON)","10(DJF)","13(MAM)","16(JJA)","19(SON)",\
;                     "1(JJA)","4(SON)","7(DJF)","10(MAM)","13(JJA)","16(SON)","19(DJF)",\
;                     "1(SON)","4(DJF)","7(MAM)","10(JJA)","13(SON)","16(DJF)","19(MAM)",\
;                     "1(DJF)","4(MAM)","7(JJA)","10(SON)","12(DJF)","16(MAM)","19(JJA)",\
;                      "1","4","7","10","12","16","19"/)
          labels = (/"MAM","JJA","SON","DJF","MAM","JJA","SON",\
                     "JJA","SON","DJF","MAM","JJA","SON","DJF",\
                     "SON","DJF","MAM","JJA","SON","DJF","MAM",\
                     "DJF","MAM","JJA","SON","DJF","MAM","JJA",\
                      "2","5","8","11","14","17","20","23"/)

  tmin = (/0,0,4/)
  tmax = (/3,3,7/)  
  res@trYMinF = 0.1
  res@trYMaxF = 1.
  res@gsnYRefLine = (/0.,0.5/)
  res@gsnYRefLineDashPatterns = (/0,0/)  
  do i=0,4
    if(i .le. 3)then
      mstr = 2
      mend = 20
      step = 3
    else
      mstr = 0
      mend = 23
      step = 1
    end if
    res@vpYF = 0.95-(res@vpHeightF+0.05)*i
    gsn_text_ndc( wks, monthtitle_l(i), res@vpXF+0.005, res@vpYF+0.019, restext)
    res@tmXBLabels = labels(i*7:i*7+6)
    res@xyLineColor = "black"
    res@xyMarkerColors = "black"
    res@xyDashPatterns = 0
    res@xyMarker = 4
    plot(i)=gsn_csm_xy(wks,r&L(mstr:mend:step),r(dimsizes(modelname)+1,0,i,mstr:mend:step),res)    ;1970-2019
    res@xyMarker = 16
    plotoverlay=gsn_csm_xy(wks,r&L(mstr:mend:step),r_sig(dimsizes(modelname)+1,0,i,mstr:mend:step),res)
    overlay(plot(i), plotoverlay)


    res@xyLineColor = "blue"
    res@xyMarkerColors = "blue"
    res@xyDashPatterns = 16
    res@xyMarker = 4
    plotoverlay=gsn_csm_xy(wks,r&L(mstr:mend:step),r(dimsizes(modelname)+2,0,i,mstr:mend:step),res)    ;1970-2019
    overlay(plot(i), plotoverlay)
    res@xyMarker = 16    
    plotoverlay=gsn_csm_xy(wks,r&L(mstr:mend:step),r_sig(dimsizes(modelname)+2,0,i,mstr:mend:step),res)
    overlay(plot(i), plotoverlay)

    res@vpYF = 0.95-(res@vpHeightF+0.05)*i
    res@xyLineColor = "red"
    res@xyMarkerColors = "red"
    res@xyDashPatterns = 16
    res@xyMarker = 4    
    plotmme(i)=gsn_csm_xy(wks,r&L(mstr:mend:step),r(dimsizes(modelname),0,i,mstr:mend:step),res)
    res@xyMarker = 16
    plotmmesig=gsn_csm_xy(wks,r&L(mstr:mend:step),r_sig(dimsizes(modelname),0,i,mstr:mend:step),res)
    overlay(plotmme(i), plotmmesig)   
    
    res3@xyLineColor = "red"
    res3@xyLineThicknessF = 0.1
    res3@xyLineOpacities = 0
    res3@gsnXYFillColors = "red"
    res@xyDashPatterns = 16
    res3@gsnXYFillOpacities = 0.1
    xxx1 := dim_min_n_Wrap(r(0:dimsizes(modelname)-1,0,i,mstr:mend:step),0)
    xxx2 := dim_max_n_Wrap(r(0:dimsizes(modelname)-1,0,i,mstr:mend:step),0)
    if(i .eq. 4)then
      print(xxx1)
      print(xxx2)
    end if
    plots = gsn_csm_xy(wks, r&L(mstr:mend:step), (/xxx1, xxx2/), res3)
    overlay(plotmme(i), plots)   
  end do
  
  restext@txFontCoblor = "black"
  gsn_text_ndc( wks, "SMYLE (1970-2019)", res@vpXF+0.01, 0.865, restext)
  restext@txFontColor = "blue"
  gsn_text_ndc( wks, "SMYLE (1982-2016)", res@vpXF+0.01, 0.845, restext)
  restext@txFontColor = "red"
  gsn_text_ndc( wks, "NMME (1982-2016)", res@vpXF+0.01, 0.825, restext)
  restext@txFontColor = "black"

  
  res@vpXF = 0.53
  res@trYMinF = 0.2
  res@trYMaxF = 1.3
  res@gsnYRefLine = (/0.,3.6/)
  res@gsnYRefLineDashPatterns = (/0,0/)  
  do i=0,4
    if(i .le. 3)then
      mstr = 2
      mend = 20
      step = 3
    else
      mstr = 0
      mend = 23
      step = 1
    end if
    res@vpYF = 0.95-(res@vpHeightF+0.05)*i
    restext@txFontColor = "black"
    gsn_text_ndc( wks, monthtitle_r(i), res@vpXF+0.005, res@vpYF+0.019, restext)
    res@tmXBLabels = labels(i*7:i*7+6)
    
    res@xyLineColor = "black"
    res@xyMarkerColors = "black"
    res@xyDashPatterns = 0
    plot(i+5)=gsn_csm_xy(wks,e&L(mstr:mend:step),e(dimsizes(modelname)+1,0,i,mstr:mend:step),res)    ;1970-2019


    res@xyLineColor = "blue"
    res@xyMarkerColors = "blue"
    res@xyDashPatterns = 16
    plotoverlay=gsn_csm_xy(wks,e&L(mstr:mend:step),e(dimsizes(modelname)+2,0,i,mstr:mend:step),res)    ;1970-2019
    overlay(plot(i+5), plotoverlay)

    res@vpYF = 0.95-(res@vpHeightF+0.05)*i
    res@xyLineColor = "red"
    res@xyMarkerColors = "red"   
    res@xyDashPatterns = 16
    plotmme(i+5)=gsn_csm_xy(wks,e&L(mstr:mend:step),e(dimsizes(modelname),0,i,mstr:mend:step),res)
    res@xyMarker = 16

    
    res3@xyLineColor = "red"
    res3@xyLineThicknessF = 0.
    res3@gsnXYFillColors = "red"
    res3@gsnXYFillOpacities = 0.1
    res@xyDashPatterns = 16
    xxx1 := dim_min_n_Wrap(e(0:dimsizes(modelname)-1,0,i,mstr:mend:step),0)
    xxx2 := dim_max_n_Wrap(e(0:dimsizes(modelname)-1,0,i,mstr:mend:step),0)

    plots = gsn_csm_xy(wks, e&L(mstr:mend:step), (/xxx1, xxx2/), res3)
    overlay(plotmme(i+5), plots)   
  end do

  draw( wks)
  frame( wks)

end 
