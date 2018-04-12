load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

; generic script to remove land
;
function removeLand(X)
local a, lsdata, lsm, data
begin
  data = X
  data@_FillValue = getVarFillValue(X)
  dims = getvardims(X)
  latDim = ind( dims.eq."lat")
  lonDim = ind( dims.eq."lon")
  ;; All values, expect for ocean points, will be set to _FillValue in the array data after creating a new land/sea mask from one-dimensional lat/lon variables:
  ; assume data is 3D (time,lat,lon)
  a    = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
  lsdata = a->LSMASK
  lsm = landsea_mask(lsdata,X&$dims(latDim)$,X&$dims(lonDim)$)

  ; lsm is a 2D array, in order to use it in mask, we must conform it
  ; to the size of the 3D array "data". 
  data = (/ where(conform(data,lsm,(/latDim,lonDim/)).ge.1, data@_FillValue, data) /)
  return(data)
end


function avg_hemis( X )
local rank, X2
begin
  if .not.isatt(X,"anti_symmetric") then
    asym = False
  else
    if X@anti_symmetric .eq. False then
      asym = False
    else
      asym = True
    end if
  end if
  if asym then
    fac = -1.
  else
    fac = 1.
  end if

  rank = dimsizes(dimsizes( X ) )
  if rank .eq. 1 then
    X2 = (X+(fac*X(::-1)))/2.
  end if

  if rank.eq.3 then
    ;; THIS IS THE  TIME LAT LON VERSION
    ;; SAME FOR THE  LEV LAT LON VERSION
    X2 = (X + (fac*X(:,::-1,:)))/2.
  end if
  copy_VarMeta( X , X2 )
  return(X2)
end

begin
;; use the old APE o3 in lower atm, but merge with modern upper levels from WACCM
wcm_o3_name ="/glade/p/cesm/chwg_dev/tilmes/climo/oxid_ozone_WACCM_CCMI_1850_2010_serial_3D_monthly.nc" 
FW = addfile(wcm_o3_name,"r")
ape_o3_name = "/glade/p/cesmdata/cseg/inputdata/atm/cam/ozone/apeozone_cam3_5_54.nc"
FA = addfile(ape_o3_name,"r")

;; inside FA
;; - P0  : reference pressure Pa
;; - lev : hybrid level at layer midpoints, positive down (59)
;; - lon (128)
;; - lat (64)
;; - date , seconds (12)
;; - datesec
;; - time
;; - OZONE_old
;; - OZONE(time,lev,lat,lon) "fraction"


;; inside FW
;; - lev  : hybrid level at midpoints (66)
;; - hyam
;; - hybm
;; - ilev
;; - hyai
;; - hybi
;; - P0
;; - time (216)
;; - date (YYYYMMDD) mid-month 1850 - 2010
;; - datesec
;; - lat (96)
;; - lon (144)
;; - slat
;; - slon
;; - time_bnds
;; - gw
;; - H2O2
;; - HO2
;; - NO3
;; - O3 : "O3 concentration" (mol/mol)
;; - OH
;; - PS

;; first thing to do is to look at the zonal mean climatology

Ozone_ape = FA->OZONE
Ozone_wcm = FW->O3

alat = FA->lat
alon = FA->lon
wlat = FW->lat
wlon = FW->lon

;; compare reference pressure:
P0a = FA->P0
P0w = FA->P0
print( "REFERENCE PRESSURE THE SAME: " + (P0a .eq. P0w) )

;; pressure levels
APE_p = FA->lev

;; need a surface pressure, use Jim's aquaplanet IC:
AIC = addfile("/glade/u/home/benedict/ys/datain/ape.IC.0003-01-01-00000.nc","r")
APE_ps = AIC->PS(0,:,:)

;; interpolate PS to WACCM grid:
APE_ps_wgrd = area_conserve_remap_Wrap( APE_ps&lon, APE_ps&lat, APE_ps, wlon, wlat, False)
print("PUT APE SURFACE PRESSURE ON WACCM GRID")

;; need the hybrid coefficients for WACCM grid:
hyao = FW->hyam
hybo = FW->hybm

;; time and zonal averages
Ozone_ape_climo = dim_avg_n_Wrap( Ozone_ape, 0 )
Ozone_wcm_climo = dim_avg_n_Wrap( Ozone_wcm, 0 ) ;; lev lat lon
Ozone_ape_avg = dim_avg_n_Wrap( Ozone_ape(lev|:, lat|:, lon|:, time|:), (/2,3/) ) ;; lev lat
Ozone_wcm_avg = dim_avg_n_Wrap( Ozone_wcm(lev|:, lat|:, lon|:, time|:), (/2,3/) )

;print(Ozone_ape_avg&lev)
;print(Ozone_wcm_avg&lev)

res = True
res@trYReverse = True
res@cnFillOn = True
res@cnFillPalette = "MPL_gnuplot"
res@gsnDraw = False
res@gsnFrame = False
res@tmYLMode = "Explicit"
res@tmYLValues = (/ 900., 700., 500., 200., 50., 25., 10., 5., 1., 0.5, 0.1, 0.05, 0.01, 0.005, 0.001/)
res@tmYLLabels = res@tmYLValues

wks = gsn_open_wks("png","~/plots/new_aqua_ozone_plot1")

plot = new(2,graphic)

res@gsnLeftString = "APE O3"
plot(0) = gsn_csm_contour( wks, Ozone_ape_avg, res )
res@gsnLeftString = "WACCM O3"
plot(1) = gsn_csm_contour( wks, Ozone_wcm_avg, res )

gsn_panel(wks, plot, (/1,2/), False)
delete( [/wks, plot/] )

;; interpolate the old data to the high-top grid
;; NOTE: APE ozone is actually interpreted as pressure levels.

;; First move to same horizontal grid
Ozone_ape_hi1 = area_conserve_remap_Wrap( alon, alat, Ozone_ape, wlon, wlat, False)
print("Put APE ozone on WACCM lat lon")
;; Then take climo avg
Ozone_ape_hi1_climo = dim_avg_n_Wrap( Ozone_ape_hi1, 0 )
;; Now interpolate vertical grid
;; --- note that P must be same units as PS and P0
Ozone_ape_hi = pres2hybrid_Wrap( APE_p, APE_ps_wgrd/100., P0a/100., Ozone_ape_hi1_climo, hyao, hybo, 0 )
Ozone_ape_hi&lev = Ozone_wcm&lev
print("Put APE ozone on WACCM lev")

Ozone_ape_hi_avg = dim_avg_n_Wrap( Ozone_ape_hi, 2 )

;; check if grids actually look the same:
printVarSummary( Ozone_ape_hi_avg )
printVarSummary( Ozone_wcm_avg )
;; same size grid

;; plot on same vert grid
wks = gsn_open_wks("png","~/plots/new_aqua_ozone_plot2")
plot = new(3,graphic)
res@gsnLeftString = "APE O3, WACCM grid"
plot(0) = gsn_csm_contour( wks, Ozone_ape_hi_avg, res )
res@gsnLeftString = "WACCM O3"
plot(1) = gsn_csm_contour( wks, Ozone_wcm_avg, res )

D = Ozone_ape_hi_avg ;; double
D = (/ Ozone_wcm_avg - Ozone_ape_hi_avg /)
res2 = res
res2@cnFillPalette = "NCV_jaisnd"
res2@gsnLeftString = "DIFFERENCE (WACCM - APE)"
plot(2) = gsn_csm_contour( wks, D, res2 )

gsn_panel(wks, plot, (/2,2/), False)
delete( [/wks, plot/] )
;; Make the symmetric version of the WACCM data
;; AND remove the land points ? -- AMBIGUOUS B/C ANTARCTICA
;; Ozone_wcm_climo = removeLand( Ozone_wcm_climo )
Ozone_wcm_climo_sym = avg_hemis( Ozone_wcm_climo )
Ozone_wcm_sym_avg = dim_avg_n_Wrap( Ozone_wcm_climo_sym, ind( getvardims(Ozone_wcm_climo_sym).eq. "lon") )
;; plot on same vert grid
wks = gsn_open_wks("png","~/plots/new_aqua_ozone_plot3")
plot = new(3,graphic)
res@gsnLeftString = "APE O3, WACCM grid"
plot(0) = gsn_csm_contour( wks, Ozone_ape_hi_avg, res )
res@gsnLeftString = "WACCM O3, symmetric"
plot(1) = gsn_csm_contour( wks, Ozone_wcm_sym_avg, res )
D = Ozone_ape_hi_avg ;; double
D = (/ Ozone_wcm_sym_avg - Ozone_ape_hi_avg /)
res2 = res
res2@cnFillPalette = "NCV_jaisnd"
res2@gsnLeftString = "DIFFERENCE (WACCM - APE)"
plot(2) = gsn_csm_contour( wks, D, res2 )
gsn_panel(wks, plot, (/2,2/), False)
delete( [/wks, plot/] )

;; NOW "BLEND" THEM
Ozone_new = Ozone_ape_hi ;; lev lat lon on WACCM, double
Ozone_new = (/ where( ismissing(Ozone_new), todouble( Ozone_wcm_climo_sym ), Ozone_new ) /)
;; make sure symmetric:
Ozone_new = (/avg_hemis( Ozone_new )/)
;; zonal average
Ozone_new_avg = dim_avg_n_Wrap( Ozone_new, 2 )

;; plot the blended product, comare to what?
wks = gsn_open_wks("png","~/plots/new_aqua_ozone_plot4")
plot = new(3,graphic)
res@gsnLeftString = "APE O3, WACCM grid"
plot(0) = gsn_csm_contour( wks, Ozone_ape_hi_avg, res )
res@gsnLeftString = "WACCM O3, no land, symmetric"
plot(1) = gsn_csm_contour( wks, Ozone_wcm_avg, res )
res@gsnLeftString = "Blended APE O3"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 1e-7
res@cnMaxLevelValF = 1e-6
res@cnLevelSpacingF = (res@cnMaxLevelValF - res@cnMinLevelValF)/25.
plot(2) = gsn_csm_contour( wks, Ozone_new_avg, res )
gsn_panel(wks, plot, (/2,2/), False)
delete( [/wks, plot/] )

;; look at equatorial profiles:
res3 = True
res3@xyLineColors = (/"black", "red", "orange"/)
res3@xyLineThicknessF = 3.
wks = gsn_open_wks("png","~/plots/new_aqua_ozone_plot5")
res3@tiMainString = "Equatorial Profiles"
res3@gsnLeftString = "WACCM:black, BLEND:red, Diff:orange"
res3@trYReverse = True
res3@trYAxisType = "LogAxis"
plot = gsn_csm_xy(wks, (/Ozone_wcm_sym_avg(:,{0}), Ozone_new_avg(:,{0}),Ozone_wcm_sym_avg(:,{0}) - Ozone_new_avg(:,{0})/), Ozone_new&lev, res3)

;; NOW WE WANT TO WRITE THE OUTPUT FILE
;; First attempt, just put on pressure levels.
pout =  (/ 5e-06, 1e-05, 1.5e-05, 2.5e-05, 4.5e-05, 7.e-05, \
            0.0001, 0.0002, 0.0003, 0.0005, \
            0.001,  0.002,  0.003,  0.005, \
            0.01,   0.02,   0.03,   0.05,  \
            0.1,    0.2,    0.3,    0.5,   \
            1.,     2.0,    3.0,    5.0,   \ 
           10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 90., \
          100., 110.,  120., 130., 140., 150., 160., 170., 180., 190., \
          200., 225., 250., 275., \
          300., 350., \
          400., 450., \
          500., 550., \
          600., 650., \
          700., 750., 775., \
          800., 850., 875., \
          900., 925., 950., 1000. /)
pout@long_name = "pressure"
pout@units = "hPa"
pout!0 = "lev"
pout&lev = pout

printVarSummary( Ozone_new )
printVarSummary( APE_ps_wgrd)

Ozone_plev = vinth2p( Ozone_new, hyao, hybo, pout, APE_ps_wgrd, 1, P0a/100., 1, True )
printVarSummary(Ozone_plev)

;; Ozone_new is (lev lat lon)
climo_size = dimsizes( Ozone_plev)
out_size = new( dimsizes(climo_size)+1, integer)
out_size(0) = 12
out_size(1::) = climo_size
print(out_size)
ozone_output = conform_dims( out_size, Ozone_plev, (/1,2,3/) )
ozone_output!0 = "time"
ozone_output!1 = "lev"
ozone_output!2 = "lat"
ozone_output!3 = "lon"
ozone_output&time = FA->time
ozone_output&lev  = pout
ozone_output&lat  = wlat
ozone_output&lon  = wlon
ozone_output@long_name = "ozone"
ozone_output@units     = "fraction (mol/mol)"
printVarSummary( ozone_output )

;; look at equatorial profiles:
res3 = True
res3@gsnDraw = False
res3@gsnFrame = False
res3@xyLineColors = (/"black", "red", "orange"/)
res3@xyLineThicknessF = 3.
wks = gsn_open_wks("png","~/plots/new_aqua_ozone_plot6")
res3@tiMainString = "Equatorial Profiles (2)"
res3@gsnLeftString = "WACCM:black, BLEND:red, PLEV:orange"
res3@trYReverse = True
res3@trYAxisType = "LogAxis"
plot = gsn_csm_xy(wks, (/Ozone_wcm_sym_avg(:,{0}), Ozone_new_avg(:,{0}) /), Ozone_new&lev, res3)
gres = True
gres@gsLineColor = "orange"
gres@gsLineThicknessF = 4.
gres@gsLineDashPattern = 2
lll = gsn_add_polyline( wks, plot, ozone_output( 0, :, {0},{0}), pout, gres )
draw(plot)
frame(wks)

;; FORM THE OUTPUT FILE:
foutname = "/glade/p/work/brianpm/model_data/aquaplanet_ozone_hightop"
foutname = foutname+"_c"+systemfunc("date +%y%m%d")+".nc"
Fout = addfile(foutname, "c")
Fout@title = "Aquaplanet ozone data set"
Fout@institution = "NCAR"
Fout@source = "Brian Medeiros <brianpm@ucar.edu>"
Fout@comment = "Blended data based on APE ozone (from AMIP II) and WACCM climatology (upper levels). \n APE ozone from "+ape_o3_name+" \n WACCM ozone from "+wcm_o3_name

Fout->OZONE   = ozone_output
Fout->P0      = FA->P0
Fout->date    = FA->date
Fout->datesec = FA->datesec

;; inside FA
;; - P0  : reference pressure Pa
;; - lev : hybrid level at layer midpoints, positive down (59)
;; - lon (128)
;; - lat (64)
;; - date , seconds (12)
;; - datesec
;; - time
;; - OZONE_old
;; - OZONE(time,lev,lat,lon) "fraction"


end
