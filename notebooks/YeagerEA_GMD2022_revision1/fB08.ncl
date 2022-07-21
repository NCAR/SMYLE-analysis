;;=================== Description ==================================
; This script reads in either iBTracs data or data from the GFDL
; tracking code and plots a track density plot for a specific ocean
; basin
;==================================================================

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "/glade/u/home/zarzycki/ncl/functions/trackDensity.ncl"
load "/glade/u/home/zarzycki/ncl/functions/surfaceWindCorrFactor.ncl"
load "/glade/u/home/zarzycki/ncl/functions/getTrajectoriesGeneric.ncl"

begin

gridsize=5.0  ; global map dimensino 19 * 3
allyears = ispan(1970, 2018,1)
nyears = dimsizes(allyears)
count =new((/4, nyears, 37, 72/),float) ; week, years, nlat, nlon


do iw = 0, 3
    leadmonth = iw*3 + 2
    print("leadmonth = "+leadmonth)

    if (leadmonth.ge.10) then
      leadmonth_str = tostring(leadmonth)
    else
      leadmonth_str = "0"+tostring(leadmonth)
    end if

    do iy = 0, nyears-1
      year = allyears(iy)
      print("year "+year)
      
        files = systemfunc("ls ../traj_files/trimmed/trajectories.txt.SMYLE."+year+"."+leadmonth_str+".0*_trim2.nc" )  ;all ensemble members
        f = addfiles(files,"r")
        ListSetType(f, "cat")
        traj_data = f[:]->traj_data
        printVarSummary(traj_data)
        lon = traj_data(:,2,:)
        lat = traj_data(:,3,:)
        yyyy = traj_data(:,7,:)
        mm= tointeger(traj_data(:,8,:))

        
        lon1d = ndtooned(lon)
        lat1d = ndtooned(lat)
        yyyymm  = 100 * ndtooned(yyyy) + ndtooned(mm)
        dims = dimsizes(lon1d)
        
        yyyymm_c = year*100 + leadmonth +1
        print("current year month = "+yyyymm_c)
        yyyymm_c1= (year+1) * 100 + leadmonth
        print("next year month = "+yyyymm_c1)
        inds = ind(yyyymm.ge.yyyymm_c.and.yyyymm.le.yyyymm_c1)
        print(dimsizes(inds))

        count(iw, iy,:,: ) = track_density(gridsize,0.0,lat1d(inds),lon1d(inds),False)

       
  delete(inds)
  delete(dims)
  delete(f)
  delete(traj_data)
  delete(lat)
  delete(lon)
  delete(yyyy)
  delete(mm)
  delete(yyyymm)
  delete(lat1d)
  delete(lon1d)
  delete(files)

end do ; end year
end do ; end lead time

count_mean0 = dim_avg_n_Wrap(count, 0) ; mean lead
count_mean = dim_avg_n_Wrap(count_mean0, 0) / 19.   ;mean year




;----------------
; now observations

files = systemfunc("ls ../../S2S/trimbasin_traj_files/trim_obs_traj_*.nc")
nyears= dimsizes(files)
print(nyears)

f  = addfiles(files, "r")
ListSetType(f, "join")

lat = f[:]->lat(:,:,:)
lon = f[:]->lon(:,:,:) ; bs, ns
wind = f[:]->wind(:,:,:)
ndims = dimsizes(lon)
print(ndims)

lat1d = ndtooned(lat)
lon1d = ndtooned(lon)
wind1d = ndtooned(wind)

inds= ind(wind1d.ne.18)
lat1d_ = lat1d(inds)
lon1d_ = lon1d(inds)
print(min(lat1d_))
print(max(lat1d_))

countobs = track_density(gridsize,0.0,lat1d_,lon1d_,False)
countobs = countobs / nyears  ; 11 ensemble, 18 years

copy_VarCoords(countobs,count_mean)  

count_diff = countobs
count_diff = count_mean - countobs

;************************************************
; create plot
;************************************************

  wks = gsn_open_wks("pdf","track_density_SMYLE")
  plot = new(2,"graphic")

  colorMap1="WhiteYellowOrangeRed"


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
  res@cnMinLevelValF        = 1.
  res@cnMaxLevelValF        = 20.
  res@cnLevelSpacingF       = 1.

  res@lbLabelStride = 2
  res@lbBoxEndCapStyle = "TriangleBothEnds"
  res@pmLabelBarOrthogonalPosF = 0.2

  res@mpMinLatF             = -40
  res@mpMaxLatF             = 40
  res@mpMinLonF             = 0
  res@mpMaxLonF             = 360
  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)*0.5
  res@mpGridAndLimbOn       = False  
  res@mpGridLineDashPattern = 2             ; Dashed lines
  res@mpGridLatSpacingF     = 5.0
  res@mpGridLonSpacingF     = 10.0
  res@tmYLLabelFontHeightF = 0.01
  cutname = "> 0 m/s"

  res@cnFillPalette=colorMap1
  res@tiMainFontHeightF = 0.015
  res@tiMainString = "BestTrack"
 plot(0) = gsn_csm_contour_map_ce(wks,countobs,res)


 res@tiMainString = "SMYLE"
 plot(1) = gsn_csm_contour_map_ce(wks,count_mean,res)


resP = True
resP@amJust           = "TopLeft"


;resP@gsnPanelMainString = ibstyr+" - "+ibenyr+" (count/5x5box/year)"     ; new resource added in NCL V6.4.0

resP@gsnPanelFigureStringsFontHeightF = 0.02
resP@gsnPanelFigureStrings= (/"a)","b)","c)"/)
resP@gsnStringFont   = "pwritx_database"
resP@gsnPanelYWhiteSpacePercent = 5
gsn_panel(wks,plot,(/2, 1/),resP)





end
