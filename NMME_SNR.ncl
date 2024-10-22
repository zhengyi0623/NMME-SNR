begin
;读取数据
f1    = systemfunc("ls /public1/home/scb6369/rzy/NMME/NASA-GMAO/*nc")
f3    = systemfunc("ls /public1/home/scb6369/rzy/NMME/CCSM3/*nc")
f3sst=addfile(f3(0),"r")
sst3  =f3sst->sst
resst=sst3(0,0,0,:,:)
n1    =dimsizes(f1)
f_sst=addfile(f1(0),"r")
sst  =f_sst->sst
f_land=addfile("NMME_land.nc","r")
land=f_land->land
time  =f_sst->S
lat  =f_sst->Y
lon  =f_sst->X
adlat=dimsizes(lat)
adlon=dimsizes(lon)
n1size=dimsizes(sst)
tim1=new(12*30, typeof(time))
sst_prii=new((/12*30,n1size(1),n1size(2),n1size(3),n1size(4)/), typeof(sst))
delete([/sst,f_sst/])
number=0

do i=0,n1-1,1
    f_sst=addfile(f1(i),"r")
    sst  =f_sst->sst
    tim  =f_sst->S
    sst_p  = sst
    size = dimsizes(tim)
    sst_p2= sst_p(:,:,:,:,:)
    do j=0,size-1,1
        tim1(j+number)=tim(j)
        sst_prii(j+number,:,:,:,:) = sst_p2(j,:,:,:,:)
    end do
    number=number+size
    delete([/sst,f_sst,tim,sst_p,sst_p2,size/])     
end do
sst_pri1 = sst_prii(12:359,:,:,:,:)  ;去除1981年
time1 = tim1(12:359)
delete(sst_prii)
copy_VarAtts(time, time1) 
delete(time)
; resst1=reshape(resst,(/n1size(3)*n1size(4)/))
; loc_nfv   = ind(ismissing(resst1)) 
; resst2=sst_pri1
; resst3=reshape(resst2,(/12*29,n1size(1),n1size(2),n1size(3)*n1size(4)/))
; do i=0,12*29-1
;    do j=0,n1size(1)-1
;      do k=0,n1size(2)-1
;     resst3(i,j,k,loc_nfv)=999
;    end do
; end do
; end do
; sst_pri2=reshape(resst3,(/12*29,n1size(1),n1size(2),n1size(3),n1size(4)/))
; delete([/resst,resst1,resst2,resst3/]) 
if (any(isnan_ieee(sst_pri1))) then
      value = 1.e20
      replace_ieeenan(sst_pri1, value, 0)
      sst_pri1@_FillValue = value
end if
sst_pri = where(abs(sst_pri1).gt.100,sst_pri1@_FillValue,sst_pri1)
sst_pri@units = "degC"
sst_pri!0="time"
sst_pri!1="lead"
sst_pri!2="member"
sst_pri!3="lat"
sst_pri!4="lon"
sst_pri&lat=lat
sst_pri&lon=lon
sst_pri&time=time1
sst_ave = dim_avg_n(sst_pri,2)
copy_VarCoords(sst_pri(:,:,0,:,:), sst_ave)
;;;;计算VMI
latS   =  20
latN   =  66
lonL   =  124    
lonR   =  261 
vmlon      = lon({lonL:lonR})
vmlat      = lat({latS:latN})
sst_selc=sst_ave(:,0,{latS:latN}, {lonL:lonR})
;copy_VarCoords(sst_pri(:,0,0,:,:), sst_selc)
sstClm   = clmMonTLL(sst_selc)               ;monthly avg(monthly climatology)   
ssta = calcMonAnomTLL (sst_selc, sstClm)    ; replace with anonamlies;Compute the anomalies from the climatology just computed
copy_VarCoords(sst_selc, ssta)
rad    = 4.0*atan(1.0)/180.0
clat   = sqrt( cos(rad*vmlat) )              ; gw for gaussian gri
wsst  = ssta                               ; copy meta data
wsst  = ssta*conform(wsst, clat, 1)
sglb = sst_ave(:,0,:,:) 
copy_VarCoords(sst_pri(:,0,0,:,:), sglb)
sstClm_g   = clmMonTLL(sglb)               ;monthly avg(monthly climatology)   
sta_g = calcMonAnomTLL (sglb, sstClm_g) 
copy_VarCoords(sglb, sta_g)
nlat   = dimsizes( lat )
rad    = 4.0*atan(1.0)/180.0
re     = 6371220.0                          
rr     = re*rad    
dlon   = abs(lon(2)-lon(1))*rr
dx     = dlon*cos(lat*rad)
dy     = new ( nlat, typeof(dx))
dy(0)         = abs(lat(2)-lat(1))*rr
dy(1:nlat-2)  = abs(lat(2:nlat-1)-lat(0:nlat-3))*rr*0.5
dy(nlat-1)    = abs(lat(nlat-1)-lat(nlat-2))*rr
area   = dx*dy                              
wsta_g = wgt_areaave(sta_g, area, 1.0, 0)
copy_VarCoords(sst_selc, wsta_g)
wsta_g@long_name = "Wgt: "+"globle SSTA"
sst_eof =wsst - conform(sst_selc, wsta_g, 0) 
copy_VarCoords(sst_selc, sst_eof)
;————————————————————————————————————EOF分析——————————————————————————
neof   = 2        ; number of EOFs
optEOF = True       
optEOF@jopt = 0   ; 
;;optEOF@jopt = 1   ; 
optETS = False
X      = sst_eof(lat|:, lon|:, time|:)
eof    = eofunc_Wrap(X, neof, optEOF)      
eof_ts1 = eofunc_ts_Wrap(X, eof, optETS)
eof_ts=new((/2,29*12/),float) 
if(eof(1,0,0).lt.0) then
    eof_ts = eof_ts1
else 
    eof_ts = (/eof_ts1(:,:)/)*-1
end if
VMI1=eof_ts(1,:)
sst_YYYYMM =cd_calendar(time1, -1)
copy_VarCoords(sst_pri(:,0,0,0,0), VMI1)
VMI = dtrend_msg_n(sst_YYYYMM, VMI1, False, False, 0)
VMIstd=dim_standardize_n_Wrap(VMI,0,0)
;;;;sel_sst
latS   =  -20
latN   =  20
lonL   =  120   
lonR   =  281 
slon      = lon({lonL:lonR})
slat      = lat({latS:latN})
dimlat=dimsizes(slat)
dimlon=dimsizes(slon)
lsst_selc=sst_pri(:,:,:,{latS:latN}, {lonL:lonR})
lssta=new((/29*12,n1size(1),n1size(2),dimlat,dimlon/),float)
do j=0,n1size(1)-1
do i=0,n1size(2)-1
lsstClm   = clmMonTLL(lsst_selc(:,j,i,:,:))               ;monthly avg(monthly climatology)   
lssta(:,j,i,:,:)= calcMonAnomTLL (lsst_selc(:,j,i,:,:), lsstClm) 
end do
end do   ; replace with anonamlies;Compute the anomalies from the climatology just computed
copy_VarCoords(lsst_selc, lssta)
rad    = 4.0*atan(1.0)/180.0
lclat   = sqrt( cos(rad*slat) )              ; gw for gaussian gri
lwsst  = lssta                               ; copy meta data
lwsst  = lssta*conform(lwsst, lclat, 3)
nino_341 =  dim_avg_n_Wrap(lwsst(:,:,:,{-5:5},{190:240}), (/3,4/))
nino_342=new((/29*12,n1size(1),n1size(2)/),float)  ;nino34
do j=0,n1size(1)-1
do i=0,n1size(2)-1
 nino= dtrend_msg_n(sst_YYYYMM, nino_341(:,j,i), False, False, 0)
 nino_342(:,j,i)=nino
end do
end do
copy_VarCoords(nino_341, nino_342)
nino_34=nino_342(time|:, member|:, lead|:)
delete(n1size)
n1size=dimsizes(nino_34)
VMI_FMA=new(29,float)
do i=0,28
   VMI_FMA(i)=(VMIstd(1+12*i)+VMIstd(2+12*i)+VMIstd(3+12*i))/3.0
end do
f_av=new((/3,n1size(2)/),float)
f_av=0
do lead=0,n1size(2)-1
    do n=0,n1size(1)-1
    do year=0,28
    f_av(0,lead)=nino_34(1+12*year,n,lead)+f_av(0,lead)
    f_av(1,lead)=nino_34(2+12*year,n,lead)+f_av(1,lead)
    f_av(2,lead)=nino_34(3+12*year,n,lead)+f_av(2,lead)
      end do
  end do
  f_av(0,lead)=f_av(0,lead)/((n1size(1)-1)*28)
  f_av(1,lead)=f_av(1,lead)/((n1size(1)-1)*28)
  f_av(2,lead)=f_av(2,lead)/((n1size(1)-1)*28)
end do

fn_av=new((/3,n1size(2),29/),float)
fn_av=0
do lead=0,n1size(2)-1
    do year=0,28
      do n=0,n1size(1)-1
     fn_av(0,lead,year)=nino_34(1+12*year,n,lead)+fn_av(0,lead,year)
     fn_av(1,lead,year)=nino_34(2+12*year,n,lead)+fn_av(1,lead,year)
     fn_av(2,lead,year)=nino_34(3+12*year,n,lead)+fn_av(2,lead,year)
      end do
  fn_av(0,lead,year)=fn_av(0,lead,year)/(n1size(1)-1)
  fn_av(1,lead,year)=fn_av(1,lead,year)/(n1size(1)-1)
  fn_av(2,lead,year)=fn_av(2,lead,year)/(n1size(1)-1)
    end do
end do
sgmn=new((/3,n1size(2)/),float)
sgmn=0
do lead=0,n1size(2)-1
    do year=0,28
      do n=0,n1size(1)-1
    sgmn(0,lead)=(nino_34(1+12*year,n,lead)-fn_av(0,lead,year))^2+sgmn(0,lead)
    sgmn(1,lead)=(nino_34(2+12*year,n,lead)-fn_av(1,lead,year))^2+sgmn(1,lead)
    sgmn(2,lead)=(nino_34(3+12*year,n,lead)-fn_av(2,lead,year))^2+sgmn(2,lead)
      end do
  end do
   sgmn(0,lead)= (sgmn(0,lead)/((n1size(1)-1)*28))^0.5
   sgmn(1,lead)= (sgmn(1,lead)/((n1size(1)-1)*28))^0.5
   sgmn(2,lead)= (sgmn(2,lead)/((n1size(1)-1)*28))^0.5
end do
sgms=new((/3,n1size(2)/),float)
sgms=0
do lead=0,n1size(2)-1
    do year=0,28
        sgms(0,lead)=(fn_av(0,lead,year)-f_av(0,lead))^2+sgms(0,lead)
        sgms(1,lead)=(fn_av(1,lead,year)-f_av(1,lead))^2+sgms(1,lead)
        sgms(2,lead)=(fn_av(2,lead,year)-f_av(2,lead))^2+sgms(2,lead)
    end do
    sgms(0,lead)= (sgms(0,lead)/(28))^0.5
    sgms(1,lead)= (sgms(1,lead)/(28))^0.5
    sgms(2,lead)= (sgms(2,lead)/(28))^0.5
end do
SNR=sgms/sgmn

loc=ind(VMI_FMA.ge.0.5.or.VMI_FMA.le.-0.5)
locd=dimsizes(loc)

vf_av=new((/3,n1size(2)/),float)
vf_av=0
do lead=0,n1size(2)-1
    do n=0,n1size(1)-1
    do year=0,locd-1
    vf_av(0,lead)=nino_34(1+12*loc(year),n,lead)+vf_av(0,lead)
    vf_av(1,lead)=nino_34(2+12*loc(year),n,lead)+vf_av(1,lead)
    vf_av(2,lead)=nino_34(3+12*loc(year),n,lead)+vf_av(2,lead)
      end do
  end do
  vf_av(0,lead)=vf_av(0,lead)/((n1size(1)-1)*(locd-1))
  vf_av(1,lead)=vf_av(1,lead)/((n1size(1)-1)*(locd-1))
  vf_av(2,lead)=vf_av(2,lead)/((n1size(1)-1)*(locd-1))
end do
vfn_av=new((/3,n1size(2),locd/),float)
vfn_av=0
do lead=0,n1size(2)-1
    do year=0,locd-1
      do n=0,n1size(1)-1
     vfn_av(0,lead,year)=nino_34(1+12*loc(year),n,lead)+vfn_av(0,lead,year)
     vfn_av(1,lead,year)=nino_34(2+12*loc(year),n,lead)+vfn_av(1,lead,year)
     vfn_av(2,lead,year)=nino_34(3+12*loc(year),n,lead)+vfn_av(2,lead,year)
      end do
  vfn_av(0,lead,year)=vfn_av(0,lead,year)/(n1size(1)-1)
  vfn_av(1,lead,year)=vfn_av(1,lead,year)/(n1size(1)-1)
  vfn_av(2,lead,year)=vfn_av(2,lead,year)/(n1size(1)-1)
    end do
end do
vsgmn=new((/3,n1size(2)/),float)
vsgmn=0
do lead=0,n1size(2)-1
    do year=0,locd-1
      do n=0,n1size(1)-1
    vsgmn(0,lead)=(nino_34(1+12*loc(year),n,lead)-vfn_av(0,lead,year))^2+vsgmn(0,lead)
    vsgmn(1,lead)=(nino_34(2+12*loc(year),n,lead)-vfn_av(1,lead,year))^2+vsgmn(1,lead)
    vsgmn(2,lead)=(nino_34(3+12*loc(year),n,lead)-vfn_av(2,lead,year))^2+vsgmn(2,lead)
      end do
  end do
   vsgmn(0,lead)= (vsgmn(0,lead)/((n1size(1)-1)*(locd-1)))^0.5
   vsgmn(1,lead)= (vsgmn(1,lead)/((n1size(1)-1)*(locd-1)))^0.5
   vsgmn(2,lead)= (vsgmn(2,lead)/((n1size(1)-1)*(locd-1)))^0.5
end do
vsgms=new((/3,n1size(2)/),float)
vsgms=0
do lead=0,n1size(2)-1
    do year=0,locd-1
        vsgms(0,lead)=(vfn_av(0,lead,year)-vf_av(0,lead))^2+vsgms(0,lead)
        vsgms(1,lead)=(vfn_av(1,lead,year)-vf_av(1,lead))^2+vsgms(1,lead)
        vsgms(2,lead)=(vfn_av(2,lead,year)-vf_av(2,lead))^2+vsgms(2,lead)
    end do
    vsgms(0,lead)= (vsgms(0,lead)/(locd-1))^0.5
    vsgms(1,lead)= (vsgms(1,lead)/(locd-1))^0.5
    vsgms(2,lead)= (vsgms(2,lead)/(locd-1))^0.5
end do
vSNR=vsgms/vsgmn

;print(SNR)
;print(vSNR)
;print(SNR(0,:))
;print(vSNR(0,:))
 path_out = "NASA-GMAO_SNR.nc"
system("rm -f " + path_out)
ncdf = addfile(path_out, "c")
; eof_ts!1 = "year"
; eof_ts@year = ispan(1979, 2012, 1)
ncdf->SNR = SNR;储存pc1 命名为VMI
ncdf->vSNR =vSNR
ncdf->VMI_FMA = VMI_FMA ;储存pc1 命名为VMI
ncdf->nino_34=nino_34
;;;;plot
; x=new(9,integer)
; y=new(9,integer)
; do i=0,8
;     x(i)=i+1
;     y(i)=1
; end do
; ;绘图
; wks = gsn_open_wks("png","SNRbar_spread1") 

; res                   = True    
; res@gsnDraw           = False         
; res@gsnFrame          = False   
; res@gsnMaximize         =True     ; 页面最大化（large format）     
; res@gsnLeftString = " "
; res@gsnRightString = " "
; ;;视图参数
; res@vpHeightF            = 0.5        ; Changes the aspect ratio
; res@vpWidthF             = 0.7
; ;res@vpXF                 = 0.1         ; change start locations
; ;res@vpYF                 = 0.75         ; the plot


; res@tmBorderThicknessF = 4.0
; res@tmYLOn = True
; ;res@tmYLFormat            = "0@;*.1f"  
; res@tmYLMajorLengthF      = 0.012           
; res@tmYLMajorThicknessF   = 4
; res@tmYLMajorOutwardLengthF = 0.012
; res@tmYLMinorLengthF      = 0.005
; res@tmYLMinorThicknessF   = 2.2
; res@tmYLMinorOutwardLengthF = 0.005

; res@tmYLLabelFontHeightF  = 0.020
; res@tmYLLabelFontThicknessF = 0.20
; ;---------------------------------------------------------
; res@tmXBOn = True
; res@tmXBLabelFontHeightF  = 0.020
; res@tmXBLabelFontThicknessF = 0.20    
; res@tmXBTickSpacingF = 5
; res@tmXBMinorOn          = False      ;不设次刻度线
; ;res@tmXBFormat = "0@;*.1f" 
; res@tmXBMajorLengthF      = 0.02         
; res@tmXBMajorThicknessF   = 4
; res@tmXBMajorOutwardLengthF = 0.012

; res@tmXBMinorLengthF      = 0.005
; res@tmXBMinorThicknessF   = 2.2
; res@tmXBMinorOutwardLengthF = 0.005
; res@tiXAxisOn              = True       
; res@tiXAxisString          = "lead month"  ;F后表示选择字体
; res@tiXAxisFontHeightF     = 0.025 
; res@tiXAxisFontThicknessF  = 1
; res@tiYAxisOn              = True       
; res@tiYAxisString          = "SNR"  ;F后表示选择字体
; res@tiYAxisFontHeightF     = 0.025 
; res@tiYAxisFontThicknessF  = 1
; res@tmXTOn    = False
; res@tmYROn    = False
; ;--设置坐标标签
; res@trYMinF               = 1.6
; res@trYMaxF               = 4
; res@tmXBMode = "Explicit"


; ;res@tmXBLabelDeltaF = -0.6   
; ;res@tmXBLabelFont = 22


; res1=True
; res1=res
; res1@gsnLeftStringFontHeightF = 0.03
; res1@gsnLeftString = "(a) February initialized forecast"
; res1@tmXBMode = "Explicit"
; lon_value = ispan(1,9,1)
; res1@tmXBValues = lon_value
; res1@tmXBLabels = (/"May","Jun",\
;             "Jul","Aug","Sep","Oct","Nov","Dec","Jan"/)
; res1@xyMarkLineMode     = "Lines"
; res1@xyLineColors      = (/"blue"/)        ; line colors
; res1@xyLineThicknesses = (/3/)                      ; line thicknesses 
; res1@xyDashPatterns    = (/0/)                      ; keep all solid
; res2=True
; res2=res
; res2@xyMarkLineMode     = "Lines"
; res2@xyLineColors      = (/"red"/)        ; line colors
; res2@xyLineThicknesses = (/3/)                      ; line thicknesses 
; res2@xyDashPatterns    = (/0/)         

; plot1 = gsn_csm_xy(wks,x,SNR(0,3:11),res1)
; overplot1  = gsn_csm_xy(wks,x,vSNR(0,3:11),res2)
; overlay(plot1,overplot1)
; res3=True
; res3=res
; res3@gsnLeftStringFontHeightF = 0.03
; res3@gsnLeftString = "(b) March initialized forecast"
; res3@tmXBMode = "Explicit"
; lon_value = ispan(1,9,1)
; res3@tmXBValues = lon_value
; res3@tmXBLabels = (/"Jun",\
;             "Jul","Aug","Sep","Oct","Nov","Dec","Jan","Feb"/)
; res3@xyMarkLineMode     = "Lines"
; res3@xyLineColors      = (/"blue"/)        ; line colors
; res3@xyLineThicknesses = (/3/)                      ; line thicknesses 
; res3@xyDashPatterns    = (/0/)                      ; keep all solid
; plot2 = gsn_csm_xy(wks,x,SNR(1,3:11),res3)
; overplot2 = gsn_csm_xy(wks,x,vSNR(1,3:11),res2)
; overlay(plot2,overplot2)
; res4=True
; res4=res
; res4@gsnLeftStringFontHeightF = 0.03
; res4@gsnLeftString = "(c) April initialized forecast"
; res4@tmXBMode = "Explicit"
; lon_value = ispan(1,9,1)
; res4@tmXBValues = lon_value
; res4@tmXBLabels = (/"Jul",\
;            "Aug","Sep","Oct","Nov","Dec","Jan","Feb","March"/)
; res4@xyMarkLineMode     = "Lines"
; res4@xyLineColors      = (/"blue"/)        ; line colors
; res4@xyLineThicknesses = (/3/)                      ; line thicknesses 
; res4@xyDashPatterns    = (/0/)                      ; keep all solid
; plot3 = gsn_csm_xy(wks,x,SNR(2,3:11),res4)
; overplot3 = gsn_csm_xy(wks,x,vSNR(2,3:11),res2)
; overlay(plot3,overplot3)

; genres                         = True
;     genres@XPosPercent             = 10                     ; move to the right
;     genres@ItemSpacePercent        = 6
;     textres                        = True
;     textres@lgLabelFontHeights     = 0.02
;     textres@lgLabels               = (/"VM","all"/)
;     textres@lgPerimOn              = False                   ; no perimeter
;     textres@lgItemCount            = 5                       ; how many
;     lineres                        = True
;     lineres@lgLineThicknesses      = 3                     ; line thickness
;     lineres@LineLengthPercent      = 8                       ; expressed as %, 0->100, length of line
;     lineres@lgLineLabelFontHeights = 0.015                   ; font height
;     ;lineres@lgDashIndexes          = (/0,1,0,1/)             ; line patterns
;     lineres@lgLineColors           = (/"red","blue"/)
; dum = simple_legend(wks, plot1, genres, lineres, textres)
; dum = simple_legend(wks, plot2, genres, lineres, textres)
; dum = simple_legend(wks, plot3, genres, lineres, textres)

; resP                  = True                   ; modify the panel plot
; resP@gsnPanelCenter       = True         ; 居中对齐每行图（默认）
; resP@gsnMaximize         = True      ; 页面最大化（large format）
; resP@gsnPanelLabelBar    = True         ; 多附图共用一个色标（add common colorbar）

; resP@lbPerimOn       = False
; resP@gsnPanelYWhiteSpacePercent = 4   ; 每两个相邻子图之间Y轴方向的空白大小
; resP@gsnPanelXWhiteSpacePercent =  4    ; 每两个相邻子图之间X轴方向的空白大小
; ;color code
; ;resP@gsnPanelFigureStrings= (/"(a)MAM(0)","(b)MAM(0)",\
;    ;                          "(c)JJA(0)","(d)JJA(0)","(e)SON(0)","(f)SON(0)","(g)DJF(1)","(h)DJF(1)"/) ;"(a)DJF(0)","(b)DJF(0)", to panel;resP@gsnPanelFigureStringsFontHeightF = 0.008
; ;resP@amJust   = "TopLeft"
; ;resP@gsnPanelFigureStringsFontHeightF = 0.008
; resP@gsnPanelRowSpec = True 
; gsn_panel(wks,(/plot1,plot2,plot3/),(/1,1,1/),resP);,p
end