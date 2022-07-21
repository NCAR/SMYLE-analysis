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

ACE = new((/8,nyears,nbasin/),float) ; lead year basin
ACE = 0.
ACE!1 = "year"
ACE!2 = "basin"
ACE!0 = "lead"  ; obs,1,4,7,10,13,16,19

ACE&lead = (/0,1,4,7,10,13,16,19/)
ACE&year = ispan(1970, 2019,1)
ACE&basin = ispan(0, 5, 1)

ACE_sh  = ACE ; winter ACE



; now observations
  files = systemfunc("ls ../../S2S/trimbasin_traj_files/trim_obs_traj_{1970..2019}*.nc")
  f = addfiles(files, "r")
  ListSetType(f, "join")  ; year, nstorm, nstep

  lat = f[:]->lat
  lon = f[:]->lon
  wind = f[:]->wind
  wind = wind / 0.51
  month = f[:]->month


  do iy = 1,nyears-1  ; 1971, 2019
    year = allyears(iy)
    nstorms =num(wind(iy,0,:).gt.0)
    print(nstorms)

    do istorm = 0, nstorms-1  ; loop over storms
      ;NH TC season
      orimon = month(iy,0,istorm)
      ace0 = tofloat(sum(wind(iy,:,istorm)^2) / 10000.)

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
            ACE(0,iy,ibasin) = ACE(0,iy,ibasin) +  ace0  ; curent year summer
          else if (orimon.eq.12) then
            ACE_sh(0, iy, ibasin) = ACE_sh(0, iy, ibasin) + ace0   ;current year winter
          else if (orimon.le.5) then
            ACE_sh(0, iy-1, ibasin) = ACE_sh(0, iy-1, ibasin) + ace0 ; last year winter
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



;--------------------------------------
; now SMYLE forecast for NH

do iw = 1, 19, 3
;  print("forecast lead "+iw)
  ; lead 1, May; 
  ; lead 4, Feb
  ; lead 7, Nov
  ; lead 10, Aug
  ; lead 13, May
  ; lead 16, Feb
  ; lead 19, Nov

  leadmonth = mod(iw * 3 +2, 12)  ; 5, 2,  11, 8, 5, 2, 11, 8
;  print("leadmonth = "+leadmonth)

  if (leadmonth.ge.10) then
    leadmonth_str = tostring(leadmonth)
  else
    leadmonth_str = "0"+tostring(leadmonth)
  end if


  do iy = 0, nyears-2  ; 1970-2018
    year = allyears(iy)

        if (iw.le.4) then   ;     year = year   , 1 and 4 months lead JJASOM, the current year
          yearnew = year
        else if (iw.ge.7.and.iw.le.16) then
          yearnew = year + 1         ; lead =7, 10, 13, 16 JJASON, the next year
        else ; eq.19
          yearnew = year + 2         ; leand = 19, JJASON, two years later
        end if
        end if

     if (yearnew.le.2019) then

      do ens= 1, 19
      
        if (ens.lt.10) then
          ENSNUM="0"+ens
        else
          ENSNUM =ens
        end if
     
       f = addfile("../traj_files/trimmed/trajectories.txt.SMYLE."+year+"."+leadmonth_str+".0"+ENSNUM+"_trim2.nc","r" ) ; this is for week 2, use the n-1th fi
       traj_data = f->traj_data ;storm, var,steps
    ;   printVarSummary(traj_data)
       lon = traj_data(:,2,:)  
       lat = traj_data(:,3,:)
       wind = traj_data(:,5,:)
       wind = wind / 0.51  ; m/s to kt
       month = traj_data(:,8,:)
       yyyy = traj_data(:,7,:)

       nstorms =num(wind(:,0).gt.0)
       print(nstorms)

       do istorm = 0, nstorms-1  ; loop over storms
         ace0 = tofloat(sum(wind(istorm,:)^2) / 10000.)

         ibasin  = 999
         orilat = lat(istorm,0)
         orilon = min((/lon(istorm,0), 359.0/))
         orimon  = month(istorm, 0)
         oriyear = yyyy(istorm,0)
         tcmask = tcMask({orilat},{orilon})

         if (orimon.ge.6.and.orimon.le.11.and.oriyear.eq.yearnew) then             ; for NH summer

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
          ACE({iw},{yearnew},{ibasin}) = ACE({iw},{yearnew},{ibasin}) + ace0
        end if
    
       end if  ; end if NH summer season for the current lead
 end do  ; end do storm
delete(f)
delete(traj_data)
delete(lat)
delete(lon)
delete(wind)
delete(month)
delete(yyyy)
end do ; end do ensemble member

end if ; end if yearnew.le.2019

end do ; end do year

end do ; end do lead 


;--------------------------------------
; now SMYLE forecast for SH

do iw = 1, 19, 3
  print("forecast lead "+iw)
 
  leadmonth = mod(24 - iw,12)

  print("leadmonth = "+leadmonth)

  if (leadmonth.ge.10) then
    leadmonth_str = tostring(leadmonth)
  else
    leadmonth_str = "0"+tostring(leadmonth)
  end if


  do iy = 0, nyears-2  ; 1970-2018
    year = allyears(iy)
    print("year = "+year)
      if (iw.lt.13) then
        yearnew = year
      else
        yearnew = year + 1
      end if

      print("yearnew = "+yearnew)

     if (yearnew.le.2019) then

      do ens= 1, 19
      
        if (ens.lt.10) then
          ENSNUM="0"+ens
        else
          ENSNUM =ens
        end if
     
       f = addfile("../traj_files/trimmed/trajectories.txt.SMYLE."+year+"."+leadmonth_str+".0"+ENSNUM+"_trim2.nc","r" ) ; this is for week 2, use the n-1th fi
       traj_data = f->traj_data ;storm, var,steps
    
       lon = traj_data(:,2,:)  
       lat = traj_data(:,3,:)
       wind = traj_data(:,5,:)
       wind = wind / 0.51  ; m/s to kt
       month = traj_data(:,8,:)
       yyyy = traj_data(:,7,:)

       nstorms =num(wind(:,0).gt.0)
     ;  print(nstorms)

       do istorm = 0, nstorms-1  ; loop over storms
         ace0 = tofloat(sum(wind(istorm,:)^2) / 10000.)
         print("ace0 = "+ace0)
         ibasin  = 999
         orilat = lat(istorm,0)
         orilon = min((/lon(istorm,0), 359.0/))
         orimon  = month(istorm, 0)
         oriyear = yyyy(istorm,0)
         tcmask = tcMask({orilat},{orilon})
        
        if ((orimon.eq.12.and.oriyear.eq.yearnew).or.(orimon.le.5.and.oriyear.eq.(yearnew+1))) then


           if (tcmask.eq.3) then  ; central pacific counts for nwp
              ibasin = 3
            else if (tcmask.eq.5) then ; indian ocean
              ibasin = 2
            else if (tcmask.eq.6) then  ; southern hemisphere
              ibasin = 4
              print("found SI")
            else if (tcmask.eq.7) then
              ibasin = 5
              print("found SP")
            else
              ibasin  = tcmask-1
            end if
           end if
           end if
           end if



        if (ibasin.ge.0.and.ibasin.le.5.and..not.ismissing(ace0)) then
          ACE_sh({iw},{yearnew},{ibasin}) = ACE_sh({iw},{yearnew},{ibasin}) + ace0

        end if
    
       end if  ; end if SH summer season for the current lead
 end do  ; end do storm
delete(f)
delete(traj_data)
delete(lat)
delete(lon)
delete(wind)
delete(month)
delete(yyyy)
end do ; end do ensemble member

end if ; end if yearnew.le.2019

end do ; end do year

end do ; end do lead 

;-----------------------


ACE(1:7,:,:) = ACE(1:7,:,:)  / 19. ; get ensemble mea
ACE_sh(1:7,:,:) = ACE_sh(1:7,:,:)  / 19. ; get ensemble mea



;--------------------------------------------------------
print("Beginning plot sequence")

 wks = gsn_open_wks("x11","forecast_interanual_at_nwp")          ; send graphics to PNG file
 gsn_define_colormap(wks, "BlGrYeOrReVi200")

     labels = (/" BestTrack"," 1:MAY" ," 4:FEB"," 7:NOV"," 10:AUG", " 13:MAY", " 16:FEB", " 19:NOV"/)
     lres = True
     lres@vpWidthF = 0.9
     lres@vpHeightF = 0.2
     lres@vpXF = .2
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
     lres@tmXBLabels = ispan(1972,2018,5)
     lres@tmXBLabelFontHeightF = 0.0205
     lres@tmXTLabelFontHeightF = 0.0205
     lres@tmYLLabelFontHeightF = 0.0225


     lres@xyDashPatterns = 0
     lres@xyLineThicknessF = 3.
     lres@xyMarkLineMode  =  "MarkLines"
     lres@xyMarkers = (/16,16,16,16,16,16,16,16/)
     lres@xyLineColors = (/"black","lightpink3","salmon","orange","palegreen4", "olivedrab", "steelblue","deepskyblue4"/)
     lres@xyMarkerColors = lres@xyLineColors
;     lres@xyMarkerSizes = (/0.02,0.02, 0.02, 0.02/)


     ACE =where(ACE.eq.0, ACE@_FillValue, ACE)
     ACE_sh =where(ACE_sh.eq.0, ACE@_FillValue, ACE_sh)

     datas = dim_standardize_n(ACE(:,:,:), 1, 1)  ; standardize as to years
     datas_sh = dim_standardize_n(ACE_sh(:,:,:),1,1)
     printVarSummary(datas)



     corr_nh = datas(1:7,0,:)
     corr_nh = 0.
     corr_sh = corr_nh

     do  id =1, 7
       corr_nh(id-1,:) = escorc_n(datas(0,:,:), datas(id,:,:), 0,0)
       corr_sh(id-1,:) = escorc_n(datas_sh(0,:,:), datas_sh(id,:,:),0,0)
     end do
     print(corr_nh)
     print(corr_sh)


df = nyears-2
psig = 0.1                       ; test significance level                     
do i = 0, 6  ; lead time order
  do ib = 0, 4   ; nbasin
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


fout = addfile("ACE_correlation_6basin.nc","c")
fout->ACE = ACE
fout->ACE_sh = ACE_sh
fout->corr_nh = corr_nh
fout->corr_sh = corr_sh


     plot = new(2,graphic)

       lres@pmLegendDisplayMode = "Always"
       lres@pmLegendSide           = "Top"               ; Change location of 
       lres@pmLegendParallelPosF   = .88                 ; move units right
       lres@pmLegendOrthogonalPosF = -0.42                ; move units down

       lres@pmLegendWidthF         = 0.1               ; Change width and
       lres@pmLegendHeightF         = 0.18                ; Change width and
       lres@lgLabelFontHeightF     = .015                 ; change font height
       lres@lgPerimOn = True
       lres@lgPerimFill ="SolidFill"
       lres@lgPerimFillColor = "White"
       lres@xyExplicitLegendLabels = labels         ; explicit labels
       lres@lgItemOrder  = (/8, 7, 6, 5, 4, 3, 2, 1,0 /) 

       lres@gsnCenterString = basin(0)
       plot(0) = gsn_csm_xy(wks,ispan(1,nyears,1), datas(:,:,0),lres)

       lres@pmLegendDisplayMode = "Never"
       lres@gsnCenterString = basin(3)
       plot(1) = gsn_csm_xy(wks,ispan(1,nyears,1), datas(:,:,3),lres)

;       lres@gsnCenterString = basin(3)
 ;      plot(2) = gsn_csm_xy(wks,ispan(1,nyears,1), datas(:,:,3),lres)


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
     resp@gsnFrame=False
     resp@gsnPanelCenter  = False
     resp@gsnPanelFigureStrings = (/"a)", "b)"/)
     resp@amJust = "TopLeft"
     gsn_panel(wks,plot,(/2,1/),resp)

     txres3               = True             
     txres3@txAngleF      = 90.          ; Rotate string clockwise 
     txres3@txFontHeightF = 0.02
;     gsn_text_ndc(wks,"standardized annual TC number",0.02,0.5,txres3)
     frame(wks)





end



