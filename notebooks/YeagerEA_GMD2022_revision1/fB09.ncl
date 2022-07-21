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


cutwind=17
maskspacing=1.0
tcMask = mask_tc(0,-999,-999,maskspacing,360)

basin = (/"NA","EP","NI","NWP","SH"/)
nbasin = dimsizes(basin)

files = systemfunc("ls ../../S2S/tc_count_files/obs_tc_count_*.nc")
;files = systemfunc("ls ./tc_count_obs_persist_20*_new.nc")
f = addfiles(files,"r")
ListSetType(f, "join")
tc_num_obs = f[:]->tc_num   ; year, basin , month
printVarSummary(tc_num_obs)

tc_num_obs_mean = dim_avg_n(tc_num_obs,0)
printVarSummary(tc_num_obs_mean)
delete(f)
delete(files)

;-------------
;--------------
; forecast data
f= addfiles("leadmonth_tc_count_monthly_1970-2018.nc","r")
ListSetType(f,"join")
tc_num = f[:]->tc_num  ; leadmonth, year, ens, basin, month
printVarSummary(tc_num)

tmp = dim_avg_n_Wrap(tc_num,1) ; year mean
tc_num_mean_smyle = dim_avg_n_Wrap(tmp,1)  ; leadmonth, basin, month

printVarSummary(tc_num_mean_smyle)
delete(tmp)
delete(f)
delete(tc_num)
;----------------

;--------------------------------------------------------
print("Beginning plot sequence")

 wks = gsn_open_wks("pdf","forecast_seasonality_SMYLE")          ; send graphics to PNG file

     ; sres = True
     ; sres@vpWidthF = 0.7
     ; sres@vpHeightF = 0.5
     ; sres@vpXF = .15
     ; sres@trXMinF = 0.3
     ; sres@trXMaxF = 7.05
     ; sres@trYMinF = 0
     ; sres@trYMaxF = 100
     ; sres@gsnDraw = True
     ; sres@gsnFrame = False
     ; sres@gsnXYBarChart = True
     ; sres@gsnXYBarChartBarWidth = 0.15           ; change bar widths
     ; sres@tmXBMode          = "Explicit"         ; explicit labels
     ; sres@tmXBValues        = fspan(1,6.5,6); (/1,2,3,4,5,6/)
     ; sres@tmXBLabels = basin
     ; sres@tmXBLabelFontHeightF = 0.0205
     ; sres@tmXTLabelFontHeightF = 0.0205
     ; sres@tmYLLabelFontHeightF = 0.0225
     ; sres@tiMainFontHeightF = 0.025
     ; sres@tiMainFont = "helvetica"
     ; sres@tiMainString ="Annual mean TC number"
     ; sres@gsnRightString = ""
     ; sres@tiYAxisString = "TC number"
        

     ; colors = (/(/219.,200.,102./),(/165.,170,28./),(/92.,144.,14/),(/14.,87.,45/),(/48.,141.,158/),(/39., 61., 141/)/)
     ; colors = colors/255.


     ; sres@gsnXYBarChartColors = colors(0,:)
     ; plot1 = gsn_csm_xy(wks,fspan(.55,6.05,6),dim_sum_n(tc_num_obs_mean,1),sres); draw each time series
     ; sres@gsnXYBarChartColors = colors(1,:); seperately, not
     ; plot2 = gsn_csm_xy(wks,fspan(0.7,6.2,6), dim_sum_n(tc_num_mean1,1), sres)

     ; sres@gsnXYBarChartColors = colors(2,:)
     ; plot3 = gsn_csm_xy(wks,fspan(0.85,6.35,6), dim_sum_n(tc_num_mean2,1), sres)

     ; sres@gsnXYBarChartColors = colors(3,:); seperately, not
     ; plot4 = gsn_csm_xy(wks,fspan(1.,6.5,6), dim_sum_n(tc_num_mean3,1), sres)

     ; sres@gsnXYBarChartColors = colors(4,:); seperately, not
     ; plot5 = gsn_csm_xy(wks,fspan(1.15,6.65,6), dim_sum_n(tc_num_mean4,1), sres)

     ; sres@gsnXYBarChartColors = colors(5,:); seperately, not
     ; plot6 = gsn_csm_xy(wks,fspan(1.3,6.8,6), dim_sum_n(tc_num_mean5,1), sres)

     ; ; sres@gsnXYBarChartColors = (/"black"/); seperately, not
     ; ; plot7 = gsn_csm_xy(wks,fspan(1.45,6.45,6), dim_sum_n(tc_num_mean6,1), sres)

     ; lbres                    = True          ; labelbar only resources
     ; lbres@vpWidthF           = 0.25           ; labelbar width
     ; lbres@vpHeightF          = 0.08           ; labelbar height
     ; lbres@lbBoxMajorExtentF  = 0.36          ; puts space between color boxes
     ; lbres@lbFillColors       = colors(0:1:-1,:)
     ; lbres@lbMonoFillPattern  = True          ; Solid fill pattern
     ; lbres@lbLabelFontHeightF = 0.02         ; font height. default is small
     ; lbres@lbLabelJust        = "CenterLeft"  ; left justify labels
     ; lbres@lbPerimOn          = False
     ; lbres@lgPerimColor  = "white"
     ; labels = (/"BestTrack","Week_1" ,"Week_2","Week_3","Week_4","Week_5"/)

     ; gsn_labelbar_ndc(wks,2,labels(0:1:-1),0.12,0.23,lbres)
     ; lbres@lbFillColors       = colors(2:3:-1,:)
     ; gsn_labelbar_ndc(wks,2,labels(2:3:-1),0.37,0.23,lbres)
     ; lbres@lbFillColors       = colors(4:5:-1,:)
     ; gsn_labelbar_ndc(wks,2,labels(4:5:-1),0.62,0.23,lbres)

     ; frame(wks)      


     
     lres = True
     lres@vpWidthF = 0.7
     lres@vpHeightF = 0.4
     lres@vpXF = .15
     lres@trXMinF = 0.5
     lres@trXMaxF = 12.5
     lres@trYMinF = 0.
     lres@trYMaxF = 5.
     lres@gsnDraw = False
     lres@gsnFrame = False
     lres@tmXMajorGrid = True 
     lres@tmXMajorGridLineColor = "grey"
     lres@tmXMajorGridLineDashPattern = 1
     lres@tmYMajorGrid = True
     lres@tmYMajorGridLineColor = "grey"
     lres@tmYMajorGridLineDashPattern = 1

     lres@tmXBMode          = "Explicit"         ; explicit labels
     lres@tmXBValues        = ispan(1,12,1); (/1,2,3,4,5,6/)
     lres@tmXBLabels = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
     lres@tmXBLabelFontHeightF = 0.0205
     lres@tmXTLabelFontHeightF = 0.0205
     lres@tmYLLabelFontHeightF = 0.0225


     lres@xyDashPatterns = 0
     lres@xyLineThicknessF = 3.5
     lres@xyLineColors = (/"black","salmon","steelblue","olivedrab","red"/);,"purple"/)

     plot = new(5,graphic)


     data = new((/5,nbasin,12/),float)
     data(0,:,:) = (/tc_num_obs_mean(0:4,:)/)
     data(1:4,:,:) = (/tc_num_mean_smyle(0:3,:,:)/)


     do ib = 0,1
       lres@tiMainString = basin(ib)
       plot(ib) = gsn_csm_xy(wks,ispan(1,12,1), data(0:4,ib,:),lres)
     end do

       lres@trYMaxF = 1.8
       lres@tiMainString = basin(2)
       plot(2) = gsn_csm_xy(wks,ispan(1,12,1), data(0:4,2,:),lres)

       lres@trYMaxF = 8.
       lres@tiMainString = basin(3)
       plot(3) = gsn_csm_xy(wks,ispan(1,12,1), data(0:4,3,:),lres)
     
       lres@pmLegendDisplayMode = "Always"
       lres@pmLegendSide           = "Top"               ; Change location of 
       lres@pmLegendParallelPosF   = .82                 ; move units right
       lres@pmLegendOrthogonalPosF = -0.64           

       lres@pmLegendWidthF         = 0.15                ; Change width and
       lres@pmLegendHeightF        = 0.18                ; height of legend.
       lres@lgLabelFontHeightF     = .02                 ; change font height
       lres@xyExplicitLegendLabels = (/" BestTrack"," 1-mon"," 4-mon"," 7-mon"," 10-mon"/)               ; explicit labels
       lres@lgItemOrder  = (/4,3,2, 1,0 /) 
       lres@trYMaxF = 7.
       lres@tiMainString = basin(4)
       plot(4) = gsn_csm_xy(wks,ispan(1,12,1), data(0:4,4,:),lres)


     resp = True
     resp@gsnPanelCenter  = False

     gsn_panel(wks,plot,(/3,2/),resp)

     print(dim_sum_n(tc_num_obs,1))



delete(tc_num_obs)

delete(tc_num)
delete(data)



end



