;_____________________________________________________
;_____________________________________________________
;
; olr_lag_season.ncl
;
; This script produces seasonal time vs. longitude
; evolution of equatorial OLR anomalies lag-regressed
; onto an intraseasonally-filtered OLR time series
; which is averaged over an Indian Ocean (IO) box 
; region and a west Pacific Ocean (WP) box region. In 
; addition, the two leading propagating intraseasonal 
; OLR modes from an empirical orthogonal function 
; analysis over the tropical WP are overlaid to assess
; their fidelity in capturing eastward and westward 
; OLR propagation.
;
; Author:      	Alex O. Gonzalez
; Date:        	October 2018
; Affiliation: 	Department of Geological and Atmospheric 
;		Sciences, Iowa State Unvieristy
; Email:	agon@iastate.edu
;
;_____________________________________________________
;_____________________________________________________
;
load "ind_season_day.ncl"

begin
print_clock("Starting program")

dataset = "NOAA"
;folder  = "../"+mainfld+"/rain/"
yr0 = 1979
yrf = 2016

yrc0 = (/1979,2009,2014/)
yrcf = (/2016,2010,2015/)
nyr = dimsizes(yrc0)

mn0  = 11
day0 = 11
mnf  =  4
dayf = 20

; filtering cutoffs
ta = 100    ; upper cutoff period [days] 
tb =  25    ; lower cutoff period [days]

; type for all non-int vars - either float or double
type = "float"

nlag = 20
olr_val = -10.

; lats and lon boundaries
lats = -10.0
latn = 10.0
lonw = 45.
lone = 205.

; create daily lag vector
lag = ispan(-nlag,nlag,1)
lag!0 = "lag"
lag&lag = lag
lag@units = "days" 
lag@long_name = "lag (day)" 
nl = dimsizes(lag)

do mm=0,nyr-1
   common_yr0 = yrc0(mm)
   common_yrf = yrcf(mm)
   print_clock("Processing "+common_yr0+"-"+common_yrf+" NOAA OLR data")

   ; read data array from file
   fn = "olr_"+tostring(tb)+"-"+tostring(ta)+"day.nc"
   fin = addfile(fn,"r") 
   timei = fin->time
   timei@calendar = "standard"

   ; find indices of season (my built-in function)
   id_seas = ind_season_day(timei,day0,dayf,mn0,mnf,common_yr0,common_yrf)
   olrlag0 = fin->olr(id_seas,{lats:latn},{lonw:lone})

   ; load reconstructed olr 1+2
   ; read data array from file
   fn = "olr_reconstr_1-2_15d.nc"
   fin1 = addfile(fn,"r") 
   timer = fin1->time
   timer@calendar = "standard"

   ; find indices of season (my built-in function)
   id_seasr = ind_season_day(timer,day0,dayf,mn0,mnf,common_yr0,common_yrf)

   ; read data array from file
   fn = "olr_reconstr_3-4_15d.nc"
   fin2 = addfile(fn,"r")
   lon2 = fin2->lon({110:lone})
   nx2  = dimsizes(lon2)

;***********************************************************
;--------------------------Base Point Series---------------------------
; only look near equator
blats = -5
blatn =  5

; calculate rainfall anom averaged over Indian Ocean box (time series)
blonw = 75
blone = 85
o_boxi = dim_avg_n_Wrap(olrlag0(:,{blats:blatn},{blonw:blone}), (/1,2/))

; calculate rainfall anom averaged over west Pacific box (time series)
blonw = 165
blone = 175
o_boxw = dim_avg_n_Wrap(olrlag0(:,{blats:latn},{blonw:blone}), (/1,2/))
delete(olrlag0)

lat   = fin->lat({lats:latn})    ; latitude
lon   = fin->lon({lonw:lone})    ; longitude
ny = dimsizes(lat)
nx = dimsizes(lon)

; read in unfiltere olr anom
fn = "olr_anom.nc"
fin = addfile(fn,"r")

io_regress = new((/ny,nx,nl/),"float")
wp_regress = new((/ny,nx,nl/),"float")
io12_regress = new((/nx2,nl/),"float")
io34_regress = new((/nx2,nl/),"float")
wp12_regress = new((/nx2,nl/),"float")
wp34_regress = new((/nx2,nl/),"float")
do l =-nlag,nlag
   o = fin->olr(id_seas+l,{lats:latn},{lonw:lone})
   o_mjo = fin1->olr(id_seasr+l,{110:lone})
   o_wpi = fin2->olr(id_seasr+l,{110:lone})
   o_boxw_mjo = dim_avg_n_Wrap(fin1->olr(id_seasr,{blonw:blone}),1)
   o_boxw_wpi = dim_avg_n_Wrap(fin2->olr(id_seasr,{blonw:blone}),1)

   io_regress(:,:,l+nlag) = regCoef_n(o_boxi,o,0,0)*olr_val
   io12_regress(:,l+nlag) = regCoef_n(o_boxi,o_mjo,0,0)*olr_val
   io34_regress(:,l+nlag) = regCoef_n(o_boxi,o_wpi,0,0)*olr_val

   wp_regress(:,:,l+nlag) = regCoef_n(o_boxw,o,0,0)*olr_val
   wp12_regress(:,l+nlag) = regCoef_n(o_boxw,o_mjo,0,0)*olr_val
   wp34_regress(:,l+nlag) = regCoef_n(o_boxw,o_wpi,0,0)*olr_val
end do
delete([/o_boxi,o_boxw,o,o_mjo,o_wpi,o_boxw_wpi,o_boxw_mjo,id_seas,id_seasr/])

; create metadata
io_regress!0 = "lat"
io_regress!1 = "lon"
io_regress!2 = "lag"
io_regress&lat = lat
io_regress&lon = lon 
io_regress&lag = lag
copy_VarMeta(io_regress,wp_regress)

io12_regress!0 = "lon"
io12_regress!1 = "lag"
io12_regress&lon = lon2 
io12_regress&lag = lag
copy_VarMeta(io12_regress,io34_regress)

wp12_regress!0 = "lon"
wp12_regress!1 = "lag"
wp12_regress&lon = lon2
wp12_regress&lag = lag
copy_VarMeta(wp12_regress,wp34_regress)

if mm .eq. 0
	; create arrays for plotted fields
	bigri = new((/nl,nx,nyr/),type)
	bigrw = new((/nl,nx,nyr/),type)

        ; create arrays for plotted fields
	bigri12 = new((/nl,nx,nyr/),type)
	bigri34 = new((/nl,nx,nyr/),type)

	; create arrays for plotted fields
	bigrw12 = new((/nl,nx,nyr/),type)
	bigrw34 = new((/nl,nx,nyr/),type)

	; create metadata for plotted fields
	bigri!0 = "lag"
	bigri!1 = "lon"
	bigri!2 = "dataset"
	bigri&lon = lon
	bigri&lag = lag 
	;bigri&dataset = name(:ndatasets-1)
	copy_VarMeta(bigri,bigrw)
        copy_VarMeta(bigri,bigri12)
	copy_VarMeta(bigri,bigri34)
        copy_VarMeta(bigrw,bigrw12)
	copy_VarMeta(bigrw,bigrw34)
end if

; to make plots of day lag vs. lon, average over lats
; and reorder dims (required by ncl func gsn_csm_hov)
ri = dim_avg_n_Wrap(io_regress,0)
bigri(:,:,mm) = ri(lag|:,lon|:)
delete([/io_regress,ri/])

bigri12(:,{110:205},mm) = io12_regress(lag|:,lon|:)
bigri34(:,{110:205},mm) = io34_regress(lag|:,lon|:)

rw = dim_avg_n_Wrap(wp_regress,0)
bigrw(:,:,mm) = rw(lag|:,lon|:)
delete([/wp_regress,rw/])

bigrw12(:,{110:205},mm) = wp12_regress(lag|:,lon|:)
bigrw34(:,{110:205},mm) = wp34_regress(lag|:,lon|:)
end do

;bigri = where(bigri.lt.0,bigri,bigri@_FillValue)
printMinMax(bigri,0)
printMinMax(bigrw,0)
printMinMax(bigrw12,0)
printMinMax(bigrw34,0)
;exit()
;***********************************************************
; create panel plot
;***********************************************************
   
   ;polyline
   a = 5.0 ;m s^{-1}
   pi = 4.*atan(1.0)
   spd = 86400.
   mpdeg = 6.371e6*pi/180.0
   a = a*spd/mpdeg
   c = -20.0
   b = c - lonw/a 
   pline = lon({lonw:lone})/a + b
   resp = True
   resp@gsLineDashPattern = 5

   fn  = "olr_lag_season"
   wks = gsn_open_wks("eps",fn)

   res = True

   res@cnLevelSelectionMode = "ManualLevels" ; manually set the contour levels with the following 3 resources
   res@cnMinLevelValF       = -6.;-1.5           ; set the minimum contour level
   res@cnMaxLevelValF       =  6.;1.5            ; set the maximum contour level
   res@cnLevelSpacingF      = 1.	     ; set the interval between contours
   res@cnLineLabelsOn       = False          ; turn off contour line labels
   res@cnFillOn             = True
   res@cnLinesOn            = False
   res@cnInfoLabelOn        = False
   
   res@vpHeightF            = 0.4

   res@tmXBLabelFontHeightF = .02
   res@tmYLLabelFontHeightF = .02

   res@lbLabelBarOn = False ; turn off individual cb's

   res@gsnDraw  = False                  ; don't draw yet
   res@gsnFrame = False                  ; don't advance frame yet
   res@gsnContourNegLineDashPattern = 1  ; make negative values dashed
   res@gsnContourZeroLineThicknessF = 0. ; make zero line disappear
   res@tmXBLabelsOn          = False

   lon_plt = ispan(20,230,10)
   lon_plt_str = tostring(lon_plt)
   iE = ind(lon_plt.lt.180)
   iW = ind(lon_plt.gt.180)
   lon_plt_str(iE) = lon_plt_str(iE)+"E"
   lon_plt_str(iW) = lon_plt_str(iW)+"W"
   lon_plt_str(iW) = ispan(170,130,10)+"W"
   
   res@tmXBMode          = "Explicit"              ; explicit labels
   res@tmXBValues        = lon_plt(::2)            ; location of labels 
   res@tmXBLabels        = lon_plt_str(::2)        ; labels themselve
   res@tmXBMinorValues   = lon_plt                 ; location of labels 

;  Create a blue to red colormap   
   blrd = (/ (/49,54,149/), (/69,117,180/),\
             (/116,173,209/), (/171,217,233/), (/224,243,248/), (/255,255,255/), \
             (/254,224,144/), (/253,174,97/), (/244,109,67/), (/215,48,39/),\
             (/165,0,38/)/)

   cm = todouble(blrd)/256.                      ; Convert to RGB values.

   xi = ispan(-5,5,1)
   xo = fspan(-5,5,13)
   cms = linint1_n_Wrap(xi,cm,False,xo,0,0)

   res@cnFillPalette = cms

   let = (/"a","b","c"/)
   label = let+") "+tostring(yrc0)+"-"+tostring(yrcf)

   ores                      = True              ; set up a second resource list
   ores@gsnDraw              = False             ; do not draw the plot
   ores@gsnFrame             = False             ; do not advance the frame
   ores@cnLevelSelectionMode = "ManualLevels" ; manually set the contour levels with the following 3 resources
   ores@cnMinLevelValF       = -6.0           ; set the minimum contour level
   ores@cnMaxLevelValF       =  6.0           ; set the maximum contour level
   ores@cnLevelSpacingF      =  1.0	      ; set the interval between contours
   ores@gsnContourNegLineDashPattern = 1      ; make negative values dashed
   ores@gsnContourZeroLineThicknessF = 0.     ; make zero line disappear
   ores@cnLineLabelsOn       = False          ; turn off contour line labels
   ores@cnInfoLabelOn        = False

   plot   = new(nyr*2,graphic)                ; create a plot array
   plotm  = new(nyr*2,graphic)                ; create a plot array
   plotxy = new(nyr*2,graphic)
   do mm=0,nyr-1
      res@gsnLeftString = label(mm)
      if mm .eq. nyr-1
         res@tmXBLabelsOn = True
      end if
      if mm .eq. 0
         res@gsnRightString = "IO base region"
      else if mm .eq. 1
	 delete(res@gsnRightString)
      end if
      end if
      plot(mm) = gsn_csm_hov(wks,bigri(:,{lonw:lone},mm),res)
   end do
   res@tmXBLabelsOn = False
   res@tmYLLabelsOn = False
   ores@tmYLLabelsOn = False
   delete(bigrw&lag@long_name)
   delete(bigrw12&lag@long_name)
   delete(bigrw34&lag@long_name)

   let = (/"d","e","f"/)
   label = let+") "+tostring(yrc0)+"-"+tostring(yrcf)
   do mm=0,nyr-1
      res@gsnLeftString = label(mm)
      if mm .eq. nyr-1
         res@tmXBLabelsOn = True
      end if
      if mm .eq. 0
         res@gsnRightString = "WP base region"
      else if mm .eq. 1
	 delete(res@gsnRightString)
      end if
      end if
   
      plot(mm+3)  = gsn_csm_hov(wks,bigrw(:,{lonw:lone},mm),res)
   end do

   bigri_all = bigri12
   bigri_all = bigri_all + bigri34

   plotm(0) = gsn_csm_contour(wks,bigri_all(:,:,0),ores)      ; create the olr plot 
   plotm(1) = gsn_csm_contour(wks,bigri_all(:,:,1),ores)      ; create the olr plot 
   plotm(2) = gsn_csm_contour(wks,bigri_all(:,:,2),ores)      ; create the olr plot 

   bigrw_all = bigrw12
   bigrw_all = bigrw_all + bigrw34

   plotm(3) = gsn_csm_contour(wks,bigrw_all(:,:,0),ores)      ; create the olr plot
   plotm(4) = gsn_csm_contour(wks,bigrw_all(:,:,1),ores)      ; create the olr plot 
   plotm(5) = gsn_csm_contour(wks,bigrw_all(:,:,2),ores)      ; create the olr plot 

do jj=0,nyr*2-1
   overlay(plot(jj),plotm(jj))
end do

   lres = True
   lres@gsLineThicknessF  = 4.0
   lres@gsLineDashPattern = 0
   lres@gsLineColor = "grey"

   y1 = fspan(-nlag,nlag,nl)
   x1 = y1
   x1 = 110.
   poly  = new(nyr*2,graphic)                ; create a plot array
do jj=0,nyr*2-1
   poly(jj) = gsn_add_polyline(wks,plot(jj),x1,y1,lres)
end do

;************************************************
; create panel
;************************************************
  resP = True               ; modify the panel plot
  resP@gsnMaximize        = True                                              
  resP@gsnPanelLabelBar   = True                               ; add common colorbar
  resP@gsnPanelRowSpec    = True
  ;resP@gsnPanelYWhiteSpacePercent = 5.                         ; add some white space between panels

  ;resP@txString           = "Lag-regressed filtered OLR (Indian Ocean)" ; add a title

  resP@lbTitleString      = "W m~S~-2"
  resP@lbTitlePosition    = "Right"  
  resP@lbTitleDirection   = "Across"
  resP@lbTitleFontHeightF = .012 
  resP@lbLabelFontHeightF = .012
  resP@pmLabelBarOrthogonalPosF = -0.017
  ;resP@gsnPanelYF = (/0.975,0.665,0.355/) ;  to control placement of panels
  resP@gsnPanelXF = (/0.1,0.53,0.1,0.53,0.1,0.53/) ;  to control placement of panels

  idp = (/0,3,1,4,2,5/)
  gsn_panel(wks,plot(idp),(/2,2,2/),resP)      ; now draw as one plot
  delete([/wks,plot/])
;***********************************************************

print_clock("Finishing program")
;print(clock_st-clock_en)
end
