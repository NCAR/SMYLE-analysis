;;=================== Description ==================================
; This script reads in the observed best track data and data from the SMYLE hindcasts,
; calculate track density,
; and compute regression to nino index
;
; NOTE: SH DJFMAM TC track for YEAR starts from december of YEAR and ends
; at May in YEAR+1, it should be SH tracks for YEAR+1. djfmam nino is for YEAR
; so when calculating correlation, we use nino(is+1:iend) and track(is:iend-1)
;==================================================================

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "/glade/u/home/zarzycki/ncl/functions/trackDensity.ncl"
load "/glade/u/home/zarzycki/ncl/functions/surfaceWindCorrFactor.ncl"
load "/glade/u/home/zarzycki/ncl/functions/getTrajectoriesGeneric.ncl"

begin

gridsize=5.0  ; track density grids -- global map dimensino 19 * 36

allyears = ispan(1970,2019,1)
nyears = dimsizes(allyears)
ntime = nyears

count =new((/7, ntime,37, 72/),float) ; lead, years, nlat, nlon
count!0 = "lead"
count!1 = "time"
count!2 = "lat"
count!3 = "lon"
countobs = count(0,:,:,:)
count&lead = (/1,4,7,10,13,16,19/)
count&time = allyears

count_sh = count  ; winter season
countobs_sh = countobs



;--------------------------------
; Track density in the observations
;-------------------------------

;   files = systemfunc("ls ../../../S2S/trimbasin_traj_files/trim_obs_traj_{1970..2019}*.nc")
;   print(files)
;   f = addfiles(files, "r")
;   ListSetType(f, "join")  ; year, nstorm, nstep

;   lat = f[:]->lat
;   lon = f[:]->lon
;   wind = f[:]->wind
;   month = f[:]->month

;   ; create a year coordinate
;   yyyy = month
;   yyyy = 0

;   do iy = 0, nyears-1
;     year =allyears(iy)
;     print("obs year "+year)
;     yyyy(iy,:,:) = year
;   end do
;   printVarSummary(yyyy)

;   lat1d = ndtooned(lat)
;   lon1d = ndtooned(lon)
;   wind1d = ndtooned(wind)
;   mon1d = ndtooned(month)
;   year1d = ndtooned(yyyy)


;   inds= ind(wind1d.ge.18)
;   lat1d_ = lat1d(inds)
;   lon1d_ = lon1d(inds)
;   mon1d_ = mon1d(inds)
;   year1d_ = year1d(inds)
;   print(year1d_(0:10))
;   delete(inds)

;   do iy = 1,nyears-1
;         year =allyears(iy)
;         monind = ind(mon1d_.ge.6.and.mon1d_.le.11.and.year1d_.eq.year)  ; same year JJASON  
;         countobs(iy,:,:) =track_density(gridsize,0.0,lat1d_(monind),lon1d_(monind),False)

;         delete(monind)
;         print("year+1 = "+(year+1))
;         monind = ind((mon1d_.eq.12.and.year1d_.eq.year).or.(mon1d_.le.5.and.year1d_.eq.year+1)) ; same year D and next year JFMAM
;         countobs_sh(iy,:,:) = track_density(gridsize,0.0,lat1d_(monind),lon1d_(monind),False)
;         delete(monind)


;   end do
; delete(f)
; delete(files)
; delete(lat)
; delete(lon)
; delete(wind)
; delete(month)
; delete(yyyy)
; delete(lat1d)
; delete(lon1d)
; delete(mon1d)
; delete(wind1d)
; delete(year1d)
; delete(lat1d_)
; delete(lon1d_)
; delete(mon1d_)
; delete(year1d_)

; printVarSummary(countobs)
; count&lat = countobs&lat
; count&lon = countobs&lon

; count_sh&lat = countobs&lat
; count_sh&lon = countobs&lon



; ;--------------------------------
; ; Track density in SMYLE
; ;-------------------------------

; ;----------------for NH TC season JJASON-----------
; do iw = 1, 19, 3
;   print("forecast lead "+iw)
;   ; lead 1, May; 
;   ; lead 4, Feb
;   ; lead 7, Nov
;   ; lead 10, Aug
;   ; lead 13, May
;   ; lead 16, Feb
;   ; lead 19, Nov

;  leadmonth = mod(iw * 3 +2, 12)  ; 5, 2,  11, 8, 5, 2, 11, 8
;   print("leadmonth = "+leadmonth)

;   if (leadmonth.ge.10) then
;     leadmonth_str = tostring(leadmonth)
;   else
;     leadmonth_str = "0"+tostring(leadmonth)
;   end if


;   ;------ read in data --------------
;   do iy = 0, nyears-2  ; 1970 - 2018
;     year = allyears(iy)
;     print("year "+year)

;     files=systemfunc("ls ../../traj_files/trimmed/trajectories.txt.SMYLE."+year+"."+leadmonth_str+".*trim2.nc") ; each file has 24 months of forecast
;     f  = addfiles(files, "r")
;     ListSetType(f, "cat") 
;     traj_data  = f[:]->traj_data   
;     printVarSummary(traj_data)

;     lon = traj_data(:,2,:)
;     lat = traj_data(:,3,:)
;     yyyy = traj_data(:,7,:)
;     mon = traj_data(:,8,:)
;     year1d = ndtooned(yyyy)
;     lon1d = ndtooned(lon)
;     lat1d = ndtooned(lat)
;     mon1d = ndtooned(mon)

;     dims = dimsizes(lon1d)
    
;     if (iw.le.4) then   ;    1 and 4 months lead for JJASOM, use the current year
;       yearnew = year
;     else if (iw.ge.7.and.iw.le.16) then
;       yearnew = year + 1         ; lead =7, 10, 13, 16 JJASON, the next year
;     else ; eq.19
;       yearnew = year + 2         ; leand = 19, JJASON, two years later
;     end if
;     end if

;     print("lead = "+iw+",forecast year for JJASON= "+year+", now ="+yearnew)
;     if (yearnew.le.allyears(nyears-1)) then
;       monind = ind(mon1d.ge.6.and.mon1d.le.11.and.year1d.eq.yearnew)  ; JJASON
;       count({iw},{yearnew},{:},{:}) = track_density(gridsize,0.0,lat1d(monind),lon1d(monind),False)
;       delete(monind)
;     end if


; delete(f)
; delete(files)
; delete(lat)
; delete(lon)
; delete(mon)
; delete(yyyy)
; delete(dims)
; delete(lon1d)
; delete(lat1d)
; delete(mon1d)
; delete(year1d)
; delete(traj_data)

; end do
; end do


; ;----------------for SH TC season DJFMAM-----------

; do iw = 1, 19, 3
;   print("forecast lead "+iw)
;   ; lead 1, Nov; 
;   ; lead 4, Aug
;   ; lead 7, May
;   ; lead 10, Feb
;   ; lead 13, Nov
;   ; lead 16, Aug
;   ; lead 19, May

;   leadmonth = mod(24 - iw,12) ; 11, 8, 5, 2, ..
;   print("winter leadmonth = "+leadmonth)

;   if (leadmonth.ge.10) then
;     leadmonth_str = tostring(leadmonth)
;   else
;     leadmonth_str = "0"+tostring(leadmonth)
;   end if


;   ;------ read in data -----------
;   do iy = 0, nyears-2  ; 1970 - 2018
;     year = allyears(iy)
;     print("year "+year)

;     files=systemfunc("ls ../../traj_files/trimmed/trajectories.txt.SMYLE."+year+"."+leadmonth_str+".*trim2.nc") ; each file has 24 months of forecast
;     f  = addfiles(files, "r")
;     ListSetType(f, "cat") 
;     traj_data  = f[:]->traj_data   
;     printVarSummary(traj_data)

;     lon = traj_data(:,2,:)
;     lat = traj_data(:,3,:)
;     yyyy = traj_data(:,7,:)
;     mon = traj_data(:,8,:)
;     year1d = ndtooned(yyyy)
;     lon1d = ndtooned(lon)
;     lat1d = ndtooned(lat)
;     mon1d = ndtooned(mon)

;     dims = dimsizes(lon1d)
    
;     if (iw.lt.13) then
;       yearnew = year
;     else
;       yearnew = year + 1
;     end if
    
;     print("winter lead = "+iw+",forecast year for DJF= "+yearnew+","+(yearnew+1))
;     if (yearnew.le.allyears(nyears-1)) then
;       monind = ind((mon1d.eq.12.and.year1d.eq.yearnew).or.(mon1d.le.5.and.year1d.eq.(yearnew+1)))   ; D of forecast year and JFMAM of the next year
;       count_sh({iw},{yearnew},{:},{:}) = track_density(gridsize,0.0,lat1d(monind),lon1d(monind),False)
;       delete(monind)
;     end if

; delete(f)
; delete(files)
; delete(lat)
; delete(lon)
; delete(mon)
; delete(yyyy)
; delete(dims)
; delete(lon1d)
; delete(lat1d)
; delete(mon1d)
; delete(year1d)
; delete(traj_data)

; end do
; end do


; count = count / 19. ; ensemble mean 
; count_sh  = count_sh / 19.


fout = addfile("tc_track_density.nc","r")
count = fout->count
count_sh =  fout->count_sh
countobs =  fout->countobs 
countobs_sh = fout->countobs_sh 

; fout->count = count
; fout->count_sh = count_sh
; fout->countobs = countobs
; fout->countobs_sh = countobs_sh
delete(fout)



;--------------------------------
; calculate nino regression
;----------------------------------

fn0 = addfile("../../nino34_series_1970-2019_ann.nc","r")
nino = fn0->nino34


fn = addfile("../../nino34_series_1970-2019_jjason.nc","r")
ninojja = fn->nino34

fn2 =addfile("../../nino34_series_1970-2019_djfmam.nc","r")
ninodjf = fn2->nino34
 

; get correlation with nino index
;----------------------------------

;fout = addfile("enso_regression_output.nc","c")

;-------------------------------
; for observations,JJASON

rcobs   = regCoef_n(ninojja(2:48), countobs(2:48,:,:), 0, 0) ; rc(nlat,mlon)  ; 1972-2018, 47 years
rcobs!0   = "lat"    ; name dimensions
rcobs!1   = "lon"
rcobs&lat = countobs&lat    ; assign coordinate values to named dimensions
rcobs&lon = countobs&lon
printVarSummary(rcobs)
dims = dimsizes(rcobs)

tval = rcobs@tval
yinterp = rcobs@yintercept

tval2d = onedtond(tval, dims)
tval2d@_FillValue = rcobs@_FillValue

rcobs_tval = tval2d

; fout->rcobs = rcobs
; fout->rcobs_tval = tval2d

delete(tval2d)
delete(tval)
delete(yinterp)
delete(dims)

;-------------------------------
; for observations,DJFMAM
rcobs_sh   = regCoef_n(ninodjf(1:48), countobs_sh(0:47,:,:), 0, 0) ; rc(nlat,mlon)  ; 1971-2018, 48 years
rcobs_sh!0   = "lat"    ; name dimensions
rcobs_sh!1   = "lon"
rcobs_sh&lat = countobs_sh&lat    ; assign coordinate values to named dimensions
rcobs_sh&lon = countobs_sh&lon
printVarSummary(rcobs_sh)
dims = dimsizes(rcobs_sh)

tval = rcobs_sh@tval
yinterp = rcobs_sh@yintercept

tval2d = onedtond(tval, dims)
tval2d@_FillValue = rcobs_sh@_FillValue
print(avg(abs(tval2d)))
print(max(abs(tval2d)))

rcobs_sh = rcobs_sh
rcobs_sh_tval = tval2d

delete(tval2d)
delete(tval)
delete(yinterp)
delete(dims)




;-----------------------------------
; for SMYLE JJASON

rc = regCoef_n(ninojja(2:48), count(:,2:48,:,:), 0, 1)
rc!1   = "lat"    ; name dimensions
rc!2   = "lon"
rc&lat = countobs&lat    ; assign coordinate values to named dimensions
rc&lon = countobs&lon
printVarSummary(rc)
dims = dimsizes(rc)

tval = rc@tval
yinterp = rc@yintercept



tval2d = onedtond(tval, dims)
tval2d@_FillValue = rc@_FillValue
tval2d!1   = "lat"    ; name dimensions
tval2d!2   = "lon"
tval2d&lat = countobs&lat    ; assign coordinate values to named dimensions
tval2d&lon = countobs&lon
print(max(abs(tval2d)))
print(avg(abs(tval2d)))

rc_tval = tval2d

delete(tval2d)
delete(tval)
delete(yinterp)
delete(dims)


;-----------------------------
; for SMYLE DJFMAM

rc_sh = regCoef_n(ninodjf(1:48), count_sh(:,0:47,:,:), 0, 1)
rc_sh!1   = "lat"    ; name dimensions
rc_sh!2   = "lon"
rc_sh&lat = countobs&lat    ; assign coordinate values to named dimensions
rc_sh&lon = countobs&lon
printVarSummary(rc_sh)
dims = dimsizes(rc_sh)

tval = rc_sh@tval
yinterp = rc_sh@yinterc_shept

tval2d = onedtond(tval, dims)
tval2d@_FillValue = rc@_FillValue
tval2d!1   = "lat"    ; name dimensions
tval2d!2   = "lon"
tval2d&lat = countobs&lat    ; assign coordinate values to named dimensions
tval2d&lon = countobs&lon
print(max(abs(tval2d)))
print(avg(abs(tval2d)))

rc_sh_tval = tval2d


rcobs = where(abs(rcobs_tval).ge.1.67, rcobs, 0)
rcobs_sh = where(abs(rcobs_sh_tval).ge.1.67, rcobs_sh, 0)

rc = where(abs(rc_tval).ge.1.67, rc, 0)
rc_sh = where(abs(rc_sh_tval).ge.1.67,rc_sh, 0)





;************************************************
; create plot
;************************************************

  wks = gsn_open_wks("pdf","track_density_nino_regression_smyle_season")
  plot = new(8,"graphic")
  plots = plot

  colorMap1="BlWhRd"
  colorMap2="BlueWhiteOrangeRed";"temp_diff_18lev"

  res                       = True     ; plot mods desired
  res@gsnDraw = False
  res@gsnFrame = False
  res@vpWidthF =0.7
  res@vpHeightF = 0.4
  res@gsnAddCyclic          = False    

  res@cnFillOn              = True     ; turn on color fill
  res@cnFillMode            = "RasterFill"       ; Raster Mode
  res@cnLinesOn             = False    ; turn of contour lines
  res@cnLineLabelsOn  = False

  res@cnLevelSelectionMode  = "ManualLevels"
  res@cnMinLevelValF        = -2.2
  res@cnMaxLevelValF        =2.2
  res@cnLevelSpacingF       = .1



  res@lbBoxEndCapStyle = "TriangleBothEnds"
  res@lbLabelStride         = 4            ; skip every other label
  res@lbBoxLinesOn          = False        ; turn off box between colors
  res@pmLabelBarOrthogonalPosF = 0.2

  res@mpMinLatF             = 0
  res@mpMaxLatF             = 50
  res@mpMinLonF             = 0
  res@mpMaxLonF             = 360
  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)*0.5
  res@mpGridAndLimbOn       = False  
  res@mpGridLineDashPattern = 2             ; Dashed lines
  res@mpGridLatSpacingF     = 5.0
  res@mpGridLonSpacingF     = 10.0
  res@tmYLLabelFontHeightF = 0.01
  res@mpFillOn  = False
  cutname = "> 0 m/s"

  res@cnFillPalette=colorMap2
  res@tiMainFontHeightF = 0.015
  res@lbLabelBarOn = False


  res2 = res
  res2@mpMinLatF             = -50
  res2@mpMaxLatF             = 0


  res@tiMainString = "BestTrack"
  res@tmXBOn = False
  res2@tmXTOn = False
  plot(0) = gsn_csm_contour_map(wks,rcobs,res)
  plots(0) = gsn_csm_contour_map(wks,rcobs_sh,res2)

 do i=0,6
  res@tiMainString = (i*3+1)+"-mon lead"
  plot(i+1) = gsn_csm_contour_map(wks, rc(i,:,:),res)
  plots(i+1) = gsn_csm_contour_map(wks, rc_sh(i,:,:),res2)
end do


;---------
; plot the basin boundaries
;----------------
  ; 1 = ATL  lat(0, 45), lon(257, 355) four corners: (257,22), (257, 45), (355,45),(355, 0), (295, 0)
  ; 2 = EPAC five corners (220,0), (220, 45), (257, 45)) + (295,0),(257, 22) 
  ; 3 = CPAC+ WPAC   (100, 0) (100, 45), (220 0) (220, 45)
  ; 5 = NIO       (30, 0) (30, 45), (100, 0) (100, 45)
  ; 6 = SI      (30,0), (30, -45) (135, 0) (135, -45)
  ; 6 = SP      (135,0), (135, -45) (290, 0) (290, -45)

lres = True
lres@gsLineColor = "olivedrab"
lres@gsLineThicknessF = 3.
;lres@gsLineDashPattern = 1

; AT
xpts = (/257., 257., 355., 355., 287., 257./)
ypts = (/22., 45., 45., 5., 5., 22./)
dum_at = new(5, graphic)
do i = 0, 4
  dum_at(i) = gsn_add_polyline(wks, plot(0), xpts(i:i+1), ypts(i:i+1), lres)
end do
delete(xpts)
delete(ypts)

; EP
xpts = (/220., 220., 257., 257., 287., 220./)
ypts = (/5., 45., 45., 22., 5., 5./)
dum_ep = new(5, graphic)
do i = 0, 4
  dum_ep(i) = gsn_add_polyline(wks, plot(0), xpts(i:i+1), ypts(i:i+1), lres)
end do
delete(xpts)
delete(ypts)


xpts = (/100., 100., 220., 220., 100./)
ypts = (/5., 45., 45., 5., 5./)
dum_nwp = new(4, graphic)
do i = 0, 3
  dum_nwp(i) = gsn_add_polyline(wks, plot(0), xpts(i:i+1), ypts(i:i+1), lres)
end do
delete(xpts)
delete(ypts)


xpts = (/30., 30., 100., 100., 30./)
ypts = (/5., 25., 25., 5., 5./)
dum_nio = new(4, graphic)
do i = 0, 3
  dum_nio(i) = gsn_add_polyline(wks, plot(0), xpts(i:i+1), ypts(i:i+1), lres)
end do
delete(xpts)
delete(ypts)


; xpts = (/30., 30., 200., 290., 290.,200., 30./)
; ypts = (/-45., -5., -5.,-5., -45., -45., -45./)
; dum_sh = new(6, graphic)
  xpts = (/30.,30., 135., 135., 30. /)
  ypts = (/-45., -5., -5., -45., -45/)
  dum_si = new(4, graphic)

do i = 0, 3
  dum_si(i) = gsn_add_polyline(wks, plots(0), xpts(i:i+1), ypts(i:i+1), lres)
end do


;--------
  xpts = (/135.,135., 290., 290., 135. /)
  ypts = (/-45., -5., -5., -45., -45/)
  dum_sp = new(4, graphic)

do i = 0, 3
  dum_sp(i) = gsn_add_polyline(wks, plots(0), xpts(i:i+1), ypts(i:i+1), lres)
end do


;------------
; add text
 txres               = True                     ; text mods desired
 txres@txFontHeightF = 0.02                     ; font smaller. default big
 txres@txFontColor = "black"
 dum_t1 = gsn_add_text(wks,plot(0),"NA",330.,10.0,txres) 
 dum_t2 = gsn_add_text(wks,plot(0),"EP",230.,40.0,txres) 
 dum_t3 = gsn_add_text(wks,plot(0),"NWP",115.,40.0,txres) 
 dum_t4 = gsn_add_text(wks,plot(0),"NI",40.,19.0,txres) 
 dum_t5 = gsn_add_text(wks,plots(0),"SI",42.,-10.0,txres) 
 dum_t6 = gsn_add_text(wks,plots(0),"SP",265.,-10.0,txres) 


print("before attaching plots")
res3 = True
res3@gsnAttachPlotsXAxis  = True
dum = new(8, graphic)
do i = 0, 7
  dum(i)= gsn_attach_plots(plot(i),(/plots(i)/),True,res3)
end do


;-----------------------------



resP = True
resP@amJust           = "TopLeft"

resP@gsnPanelLabelBar = True
resP@lbLabelFontHeightF = 0.008
resP@pmLabelBarOrthogonalPosF = -0.01
resP@lbBoxEndCapStyle = "TriangleBothEnds"
resP@lbLabelStride         = 4            ; skip every other label
;resP@pmLabelBarOrthogonalPosF = 0.2
;resP@gsnPanelMainString = ibstyr+" - "+ibenyr+" (count/5x5box/year)"     ; new resource added in NCL V6.4.0

resP@gsnPanelFigureStringsFontHeightF = 0.01
resP@gsnPanelFigureStrings= (/"a)","b)","c)","d)","d)","f)","g)", "h)","i)"/)
resP@gsnStringFont   = "pwritx_database"
resP@gsnPanelYWhiteSpacePercent = 5
gsn_panel(wks,plot,(/4,2/),resP)



end
