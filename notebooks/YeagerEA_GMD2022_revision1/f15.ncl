;========================================================
; This file reads in a traj file produced by the TC tracker
; and produces a lat-lon plot of trajectories color-coded by intensity
;========================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "/glade/u/home/zarzycki/ncl/functions/surfaceWindCorrFactor.ncl"
load "/glade/u/home/zarzycki/ncl/functions/getTrajectoriesGeneric.ncl"
load "/glade/u/home/zarzycki/ncl/functions/mask_tc.ncl"

;========================================================
begin


maskspacing=1.0
tcMask = mask_tc(0,-999,-999,maskspacing,360)


allyears = ispan(1970,2019,1)
nyears = dimsizes(allyears)

basin = (/"NA","EP","NI","NWP","SI","SP"/)
nbasin = dimsizes(basin)


fout = addfile("ACE_correlation_6basin.nc","r")
ACE = fout->ACE 
ACE_sh = fout->ACE_sh 
corr_nh = fout->corr_nh
corr_sh = fout->corr_sh


datas = dim_standardize_n(ACE(:,:,:), 1, 1)  ; standardize as to years
datas!0 = "lead"
datas!1 = "year"
datas!2 = "basin"

datas_sh = datas
datas_sh = (/dim_standardize_n(ACE_sh(:,:,:),1,1)/)
printVarSummary(datas)
printVarSummary(datas_sh)

print(corr_sh)



psig = 0.1                       ; test significance level                     
do i = 0, 6  ; lead time order
  do ib = 0, 3   ; nbasin
    df = num(.not.ismissing(ACE(i, :,ib))) - 2
    t    = corr_nh(i,ib)*sqrt(df/(1-corr_nh(i,ib)^2))  
    p    = student_t(t, df)
    print("t="+t+"  p="+p)            ; t=2.02755  p=0.0732238
   
    if (p.le.psig) then
      print("lead "+i+",basin= "+ib+", df = "+df+", r="+corr_nh(i, ib)+" is significant at the 90% level")
    else
      print("lead "+i+", basin="+ib+", df = "+df+", r="+corr_nh(i, ib)+" is NOT significant at the 90% level")
    end if
  end do
end do




do i = 0, 6  ; lead time order
  do ib = 4,5
    df = num(.not.ismissing(ACE_sh(i, :,ib))) - 2
    t    = corr_sh(i,ib)*sqrt(df/(1-corr_sh(i,ib)^2))  
    p    = student_t(t, df)
    print("t="+t+"  p="+p)            ; t=2.02755  p=0.0732238
   
    if (p.le.psig) then
      print("lead "+i+",basin= "+ib+", df = "+df+", r="+corr_sh(i, ib)+" is significant at the 90% level")
    else
      print("lead "+i+", basin="+ib+", df = "+df+", r="+corr_sh(i, ib)+" is NOT significant at the 90% level")
    end if
  end do
end do





;--------------------------------------------------------
print("Beginning plot sequence")

 wks = gsn_open_wks("pdf","forecast_interanual_at_ep_nwp_sp_3line")          
 gsn_define_colormap(wks, "BlGrYeOrReVi200")

     labels = (/"BestTrack","1-mon" ,"4-mon","7-mon"," 10-mon", " 13-mon", " 16-mon", " 19-mon"/)
     lres = True
     lres@vpWidthF = 0.8
     lres@vpHeightF = 0.2
     lres@vpXF = .1
     lres@vpYF = 0.8
     lres@trXMinF = 1.
     lres@trXMaxF = nyears
     lres@trYMinF = -4.
     lres@trYMaxF = 5.
     lres@gsnDraw = False
     lres@gsnFrame = False
     lres@gsnYRefLine = 0
     lres@gsnYRefLineDashPattern = 1
     lres@gsnYRefLineColor = "gray"

     lres@tmXTOn = False
     lres@tmXBMode          = "Explicit"         ; explicit labels
     lres@tmXBValues        = ispan(1,nyears,5); (/1,2,3,4,5,6/)
     lres@tmXBLabels = ispan(1970,2019,5)
     lres@tmXBLabelFontHeightF = 0.0205
     lres@tmXTLabelFontHeightF = 0.0205
     lres@tmYLLabelFontHeightF = 0.0225
     lres@tmXMajorGrid = True
     lres@tmYMajorGrid = True
     lres@tmXMajorGridLineDashPattern =1
     lres@tmYMajorGridLineDashPattern =1
     lres@tmXMajorGridLineColor = "grey"
     lres@tmYMajorGridLineColor = "grey"


     lres@xyDashPatterns = 0
     lres@xyLineThicknesses = (/5,3,3,3,3,3,3,3/)
     lres@xyMarkLineMode  =  "MarkLines"
     lres@xyMarkers = (/16,16,16,16,16,16,16,16/)
     colors =  (/"black","red","orange","salmon","palegreen4", "olivedrab", "steelblue","deepskyblue4"/)
;     lres@xyMarkerSizes = (/0.02,0.02, 0.02, 0.02/)



     plot = new(4,graphic)

;************************************************
; set legend resources for simple_legend_ndc
;************************************************
    genres                           = True
    pstart = 24
    pdelta = 18
    genres@XPosPercent               = pstart                      ; orientation on page
    genres@YPosPercent               = 10
    genres@ItemSpacePercent          = 3
    genres@LineLabelWhiteSpacePercent = 1
    textres                          = True
    textres@lgLabels                 = labels(0)
    textres@lgPerimOn                = False                   ; no perimeter
    textres@lgItemCount              = 2
    lineres                          = True
    lineres@lgLineLabelFontHeightF   = 0.015                   ; font height
    lineres@LineLengthPercent        = 7
    lineres@lgDashIndexes            = 0 
    lineres@lgLineColors             = colors(0)
    lineres@lgLineThicknesses        =3.5

     simple_legend_ndc(wks, genres, lineres, textres)

    delete(textres@lgLabels)
    delete(lineres@lgLineColors)
    genres@XPosPercent               = pstart + pdelta    
    genres@YPosPercent               = 10
    textres@lgLabels                 = labels(1:2)
    lineres@lgLineColors             = colors(1:2)

     simple_legend_ndc(wks, genres, lineres, textres)

    genres@XPosPercent               =  pstart + 2*pdelta                     ; orientation on page
    genres@YPosPercent               = 10
    textres@lgLabels                 = (/labels(5),labels(7)/)
    lineres@lgLineColors             = (/colors(5),colors(7)/)

     simple_legend_ndc(wks, genres, lineres, textres)



       lres@gsnCenterString = basin(0)
       ; for NA, only plot obs, month-1 and month-19
       selmon = (/0, 1, 7/)
       lres@xyLineColors = colors(selmon)
       lres@xyMarkerColors = lres@xyLineColors
       plot(0) = gsn_csm_xy(wks,ispan(1,nyears,1), datas(selmon,:,0),lres)
       lres@pmLegendDisplayMode = "Never"
       lres@gsnCenterString = basin(1)

       ; for EP 
       selmon = (/0,1,2/)
       lres@xyLineColors = colors(selmon)
       lres@xyMarkerColors = lres@xyLineColors
       plot(1) = gsn_csm_xy(wks,ispan(1,nyears,1), datas(selmon,:,1),lres)


       ; for NWP
       selmon = (/0,1,5/)
       lres@xyLineColors = colors(selmon)
       lres@xyMarkerColors = lres@xyLineColors
       lres@gsnCenterString = basin(3)
       plot(2) = gsn_csm_xy(wks,ispan(1,nyears,1), datas(selmon,:,3),lres)

       ; for sp
       selmon = (/0,1,2/)
       lres@xyLineColors = colors(selmon)
       lres@xyMarkerColors = lres@xyLineColors
       lres@gsnCenterString = basin(5)
       plot(3) = gsn_csm_xy(wks,ispan(1,nyears,1), datas_sh(selmon,:,5),lres)


;        lres@trYMaxF = 10
;        lres@tiMainString = basin(2)
;        plot(2) = gsn_csm_xy(wks,ispan(1,18,1), data(0:2,:,2),lres)

;        lres@trYMaxF = 50
;        lres@tiMainString = basin(3)
;        plot(3) = gsn_csm_xy(wks,ispan(1,18,1), data(0:2,:,3),lres)
     
;        lres@trYMaxF = 120
;        lres@tiMainString = basin(5)
;        plot(5) = gsn_csm_xy(wks,ispan(1,18,1), data(0:2,:,5),lres)


     resp = True
     resp@gsnFrame           = False
     resp@gsnFrame=False
     resp@gsnPanelCenter  = False
     resp@gsnPanelFigureStrings = (/"a)", "b)","c)","d)"/)
     resp@amJust = "TopLeft"
     ; resp@gsnPanelYF = (/0.95,0.6, 0.25/)
     ; resp@gsnPanelXF = (/0.1, 0.1, 0.1/)
     resp@gsnPanelBottom = 0.1
     gsn_panel(wks,plot,(/4,1/),resp)





     txres3               = True             
     txres3@txAngleF      = 90.          ; Rotate string clockwise 
     txres3@txFontHeightF = 0.02
     gsn_text_ndc(wks,"standardized basin-scale TC ACE",0.12,0.55,txres3)

     frame(wks)




end



