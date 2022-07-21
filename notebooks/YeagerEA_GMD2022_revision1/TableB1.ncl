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


tcnum = new((/nyears,nbasin/),float) ; lead year basin
tcnum = 0.
tcnum!0 = "year"
tcnum!1 = "basin"



tcnum&year = ispan(1970, 2019,1)
tcnum&basin = ispan(0, 5, 1)

tcnum_sh  = tcnum ; winter tcnum


; now observations
  files = systemfunc("ls ../../S2S/trimbasin_traj_files/trim_obs_traj_{1970..2019}*.nc")
  f = addfiles(files, "r")
  ListSetType(f, "join")  ; year, nstorm, nstep

  lat = f[:]->lat
  lon = f[:]->lon
  wind = f[:]->wind
  month = f[:]->month

  do iy = 1,nyears-1  ; 1971, 2019
    year = allyears(iy)
    nstorms =num(wind(iy,0,:).gt.0)
    print(nstorms)

    do istorm = 0, nstorms-1  ; loop over storms
      ;NH TC season
      orimon = month(iy,0,istorm)
      ibasin  = 999
      orilat = lat(iy,0,istorm)
      orilon = min((/lon(iy,0,istorm), 359.0/))
      tcmask = tcMask({orilat},{orilon})
            
           if (tcmask.eq.3) then  ; central pacific counts for nwp
              ibasin = 3
            else if (tcmask.eq.5) then ; indian ocean
              ibasin = 2
            else if (tcmask.eq.6) then  ; southern hemisphere
              ibasin = 4
            else if (tcmask.eq.7) then
              ibasin = 5
            else
              ibasin  = tcmask-1
            end if
          end if
          end if
        end if

        if (ibasin.ge.0.and.ibasin.le.5) then

          if (orimon.ge.6.and.orimon.le.11) then
            tcnum(iy,ibasin) = tcnum(iy,ibasin) +  1
          else if (orimon.eq.12) then
            tcnum_sh(iy, ibasin) = tcnum_sh( iy, ibasin) + 1
          else if (orimon.le.5) then
            tcnum_sh(iy-1, ibasin) = tcnum_sh(iy-1, ibasin) + 1
          end if
          end if
          end if

;        delete(inds)
       end if
    end do  ; end do strom
   end do ; end do year
delete(f)
delete(lat)
delete(lon)
delete(wind)
delete(month)
delete(files)

;-------------------


;--------------
; forecast data
f= addfile("./leadmonth_5_tc_count_1970-2018_6basin.nc","r")
tc_num = f->tc_num   ; lead, year, ens, basin
tc_num_sh = f->tc_num_sh
printVarSummary(tc_num)

tc_num_mean_nh = dim_avg_n(tc_num,2)  ; ens mean: lead, year, basin
tc_num_mean_sh = dim_avg_n(tc_num_sh, 2)
printVarSummary(tc_num_mean_nh)

print(tc_num_mean_nh(0,:,0))


tc_num_mean_nh = where(tc_num_mean_nh.eq.0, tc_num_mean_nh@_FillValue, tc_num_mean_nh)
tc_num_mean_sh = where(tc_num_mean_sh.eq.0, tc_num_mean_nh@_FillValue, tc_num_mean_sh)

do ilead = 0, 6
  print("ilead = "+ilead)
  print(tc_num_mean_sh(ilead,:,4:5))
end do

delete(f)

;--------------------------------------------------------
print("Beginning plot sequence")

 wks = gsn_open_wks("x11","forecast_interanual_count_5")          ; send graphics to PNG file


;      labels = (/" BestTrack"," 1:MAY" ," 4:FEB"," 7:NOV"," 10:AUG", " 13:MAY", " 16:FEB", " 19:NOV"/)

;      lres = True
;      lres@vpWidthF = 0.9
;      lres@vpHeightF = 0.2
;      lres@vpXF = .15
;      lres@trXMinF = 1.
;      lres@trXMaxF = 48.
;      lres@trYMinF = -4.
;      lres@trYMaxF = 5.
;      lres@gsnDraw = False
;      lres@gsnFrame = False
;      lres@gsnYRefLine = 0
;      lres@gsnYRefLineDashPattern = 1
;      lres@gsnYRefLineColor = "gray"

;      lres@tmXTOn = False
;      lres@tmXBMode          = "Explicit"         ; explicit labels
; ;     lres@tmXBValues        = ispan(0,48,2); (/1,2,3,4,5,6/)
; ;     lres@tmXBLabels = ispan(1970,2018,2)
;      lres@tmXBLabelFontHeightF = 0.0205
;      lres@tmXTLabelFontHeightF = 0.0205
;      lres@tmYLLabelFontHeightF = 0.0225


;      lres@xyDashPatterns = 0
;      lres@xyLineThicknessF = 3.
;      lres@xyMarkLineMode  =  "MarkLines"
;      lres@xyMarkers = (/16,16,16,16,16,16,16,16/)
;      lres@xyLineColors = (/"black","lightpink3","salmon","orange","palegreen4", "olivedrab", "steelblue","deepskyblue4"/)
;      lres@xyMarkerColors = lres@xyLineColors


     data_nh = new((/8,49,nbasin/),float)  ; data_nhset, year, nbasin  ; 1972-2018
     data_sh = data_nh
     data_nh(0,:,:) = (/tcnum(0:48,0:nbasin-1)/)
     data_nh(1:7,:,:) = (/tc_num_mean_nh(:,0:48,:)/)

     data_sh(0,:,:) = (/tcnum_sh(0:48,0:nbasin-1)/)
     data_sh(1:7,:,:) = (/tc_num_mean_sh(:,0:48,:)/)
     
; standardize
     datas = dim_standardize_n(data_nh, 1, 1)
     datas_sh = dim_standardize_n(data_sh, 1, 1)

     print(datas_sh)

     corr_nh = data_nh(1:7,0,:)
     corr_nh = 0.
     corr_sh = corr_nh

do  id =1, 7
  corr_nh(id-1,:) = escorc_n(datas(0,:,:), datas(id,:,:), 0,0)
  corr_sh(id-1,:) = escorc_n(datas_sh(0,:,:), datas_sh(id,:,:),0,0)
end do
print(corr_nh)
print(corr_sh)



df = 47
psig = 0.1                       ; test significance level                     
do i = 0, 6  ; lead time order
  do ib = 0,5   ; nbasin
    t    = corr_nh(i,ib)*sqrt(df/(1-corr_nh(i,ib)^2))  
    p    = student_t(t, df)
    print("t="+t+"  p="+p)            ; t=2.02755  p=0.0732238
   
    if (p.le.psig) then
      print("lead "+i+",basin= "+ib+", r="+corr_nh(i, ib)+" is significant at the 90% level")
    else
      print("lead "+i+", basin="+ib+", r="+corr_nh(i, ib)+" is NOT significant at the 90% level")
    end if
  end do
end do

nyears = 49


     labels = (/"BestTrack","1-mon: MAY" ,"4-mon: FEB","7-mon: NOV"," 10-mon: AUG", " 13-mon: MAY", " 16-mon: FEB", " 19-mon: NOV"/)
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
     print(nyears)
     print(ispan(1971,2019,5))
     lres@tmXBLabels = ispan(1971,2019,5)
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
     lres@xyLineColors = (/"black","lightpink3","salmon","orange","palegreen4", "olivedrab", "steelblue","deepskyblue4"/)
     lres@xyMarkerColors = lres@xyLineColors
;     lres@xyMarkerSizes = (/0.02,0.02, 0.02, 0.02/)



     plot = new(3,graphic)

;************************************************
; set legend resources for simple_legend_ndc
;************************************************
    genres                           = True
    genres@XPosPercent               = 10                      ; orientation on page
    genres@YPosPercent               = 10
    genres@ItemSpacePercent          = 3
    genres@LineLabelWhiteSpacePercent = 1
    textres                          = True
    textres@lgLabels                 = labels(0:1)
    textres@lgPerimOn                = False                   ; no perimeter
    textres@lgItemCount              = 2
    lineres                          = True
    lineres@lgLineLabelFontHeightF   = 0.015                   ; font height
    lineres@LineLengthPercent        = 7
    lineres@lgDashIndexes            = 0 
    lineres@lgLineColors             = (/"black","lightpink3"/);,"salmon","orange"/);,"palegreen4", "olivedrab", "steelblue","deepskyblue4"/)
    lineres@lgLineThicknesses        =3.5

     simple_legend_ndc(wks, genres, lineres, textres)

    genres@XPosPercent               = 32                      ; orientation on page
    genres@YPosPercent               = 10
    textres@lgLabels                 = labels(2:3)
    lineres@lgLineColors             = (/"salmon","orange"/);,"palegreen4", "olivedrab", "steelblue","deepskyblue4"/)

     simple_legend_ndc(wks, genres, lineres, textres)

    genres@XPosPercent               =  54                     ; orientation on page
    genres@YPosPercent               = 10
    textres@lgLabels                 = labels(4:5)
    lineres@lgLineColors             = (/"palegreen4", "olivedrab"/);, "steelblue","deepskyblue4"/)

     simple_legend_ndc(wks, genres, lineres, textres)

    genres@XPosPercent               = 76                      ; orientation on page
    genres@YPosPercent               = 10
    textres@lgLabels                 = labels(6:7)
    lineres@lgLineColors             = (/ "steelblue","deepskyblue4"/)

     simple_legend_ndc(wks, genres, lineres, textres)


       lres@gsnCenterString = basin(0)
       plot(0) = gsn_csm_xy(wks,ispan(1,nyears,1), datas(:,:,0),lres)

       lres@pmLegendDisplayMode = "Never"
       lres@gsnCenterString = basin(1)

       plot(1) = gsn_csm_xy(wks,ispan(1,nyears,1), datas(:,:,1),lres)

       lres@gsnCenterString = basin(3)
       plot(2) = gsn_csm_xy(wks,ispan(1,nyears,1), datas(:,:,3),lres)


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
     resp@gsnPanelFigureStrings = (/"a)", "b)","c)"/)
     resp@amJust = "TopLeft"
     ; resp@gsnPanelYF = (/0.95,0.6, 0.25/)
     ; resp@gsnPanelXF = (/0.1, 0.1, 0.1/)
     resp@gsnPanelBottom = 0.1
     gsn_panel(wks,plot,(/3,1/),resp)





     txres3               = True             
     txres3@txAngleF      = 90.          ; Rotate string clockwise 
     txres3@txFontHeightF = 0.02
     gsn_text_ndc(wks,"standardized basin-scale TC number",0.03,0.55,txres3)

     frame(wks)








;      plot = new(3,graphic)

;        lres@pmLegendDisplayMode = "Always"
;        lres@pmLegendSide           = "Top"               ; Change location of 
;        lres@pmLegendParallelPosF   = .88                 ; move units right
;        lres@pmLegendOrthogonalPosF = -0.53                ; move units down

;        lres@pmLegendWidthF         = 0.12                ; Change width and
;        lres@pmLegendHeightF         = 0.08                ; Change width and
;        lres@lgLabelFontHeightF     = .02                 ; change font height
;        lres@lgPerimOn = False
;        lres@xyExplicitLegendLabels = labels(0:3)         ; explicit labels
;        lres@lgItemOrder  = (/2, 1,0 /) 

;        lres@gsnCenterString = basin(0)
;        plot(0) = gsn_csm_xy(wks,ispan(0,48,1), datas(:,:,0),lres)

;        lres@pmLegendDisplayMode = "Never"
;      do ib = 1,4
;        lres@gsnCenterString = basin(ib)
;        plot(ib) = gsn_csm_xy(wks,ispan(0,48,1), datas(:,:,ib),lres)
;      end do

; ;        lres@trYMaxF = 10
; ;        lres@tiMainString = basin(2)
; ;        plot(2) = gsn_csm_xy(wks,ispan(1,18,1), data(0:2,:,2),lres)

; ;        lres@trYMaxF = 50
; ;        lres@tiMainString = basin(3)
; ;        plot(3) = gsn_csm_xy(wks,ispan(1,18,1), data(0:2,:,3),lres)
     
; ;        lres@trYMaxF = 120
; ;        lres@tiMainString = basin(5)
; ;        plot(5) = gsn_csm_xy(wks,ispan(1,18,1), data(0:2,:,5),lres)


;      resp = True
;      resp@gsnFrame=False
;      resp@gsnPanelCenter  = False
;      resp@gsnPanelFigureStrings = (/"a)", "b)","c)","d)","e)"/)
;      resp@amJust = "TopLeft"
;      gsn_panel(wks,plot,(/5,1/),resp)

;      txres3               = True             
;      txres3@txAngleF      = 90.          ; Rotate string clockwise 
;      txres3@txFontHeightF = 0.02
;      gsn_text_ndc(wks,"standardized annual TC number",0.15,0.5,txres3)
;      frame(wks)





end



