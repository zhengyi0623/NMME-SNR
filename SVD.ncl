begin
latS   =  0
latN   =  66
lonL   =  120
lonR   =  280
;读取数据
f    = systemfunc("ls /THL7/home/dingruiqiang/rzy/sst/*nc")
f1   = systemfunc("ls /THL7/home/dingruiqiang/rzy/u/*nc")
f2   = systemfunc("ls /THL7/home/dingruiqiang/rzy/v/*nc")
n    =dimsizes(f)
f_sst=addfile(f(0),"r")
lat  =f_sst->Y
lon  =f_sst->X
sst  =f_sst->sst
tim =f_sst->S
slon      = lon({lonL:lonR})
sizlon=dimsizes(sizlon)
slat      = lat({latS:latN})
sizlat=dimsizes(sizlat)
sst_sell=new((/5,10,12,sizlat,sizlon/), typeof(sst))
sst_prii=new((/5*n,10,12,sizlat,sizlon/), typeof(sst))
sst_hindd=new((/5*n*12,10,sizlat,sizlon/), typeof(sst))
u_sell=new((/5,10,12,sizlat,sizlon/),typeof(sst))
u_prii=new((/5*n,10,12,sizlat,sizlon/), typeof(sst))
u_hind=new((/5*n*12,10,sizlat,sizlon/), typeof(sst))
v_sell=new((/5,10,12,sizlat,sizlon/),typeof(sst))
v_prii=new((/5*n,10,12,sizlat,sizlon/), typeof(sst))
v_hind=new((/5*n*12,10,sizlat,sizlon/), typeof(sst))
time=new(60*n, typeof(tim))
do i=0,n-1,1
    f_sst=addfile(f(i),"r")
    f_u=addfile(f1(i),"r")
    f_v=addfile(f2(i),"r")
    sst  =f_sst->sst
    u  =f_u->u(:,0,:,:,:,:)
    v =f_v->v(:,0,:,:,:,:)
    tim =f_sst->S
    do t=0,4,1 
    sst_sell(t,:,:,:,:) =sst(1+12*t,:,:,{latS:latN}, {lonL:lonR}) ;提取2月开始预测数据
    u_sell(t,:,:,:,:) =u(1+12*t,:,:,{latS:latN}, {lonL:lonR})
    v_sell(t,:,:,:,:) =v(1+12*t,:,:,{latS:latN}, {lonL:lonR})
    end do
    sst_pri1  = sst_sell - 273.15       ;转为摄氏度
    do j=0,4,1
        sst_prii(j+i*5,:,:,:,:) =sst_pri1(j,:,:,:,:)  ;将多个文件数据变为一个
        u_prii(j+i*5,:,:,:,:) =u_sell(j,:,:,:,:)
        v_prii(j+i*5,:,:,:,:) =v_sell(j,:,:,:,:)
    end do
    do j=0,59,1
        time(j+i*60)=tim(j)
    end do     
end do
do i=0,5*n-1,1
    do j=0,11,1
    sst_hindd(12*i+j,:,:,:)=sst_prii(i,:,j,:,:)   ;按时间顺序排序
    u_hind(12*i+j,:,:,:)=u_prii(i,:,j,:,:)
    u_hind(12*i+j,:,:,:)=v_prii(i,:,j,:,:)
    end do
end do
time(:)=time(:)+1
sst_hind = where(abs(sst_hindd).eq.273.15,sst_hindd@_FillValue,sst_hindd)
sst_hind@units = "degC"
sst_hind!0="time"
sst_hind!1="member"
sst_hind!2="lat"
sst_hind!3="lon"
u_hind!0="time"
u_hind!1="member"
u_hind!2="lat"
u_hind!3="lon"
v_hind!0="time"
v_hind!1="member"
v_hind!2="lat"
v_hind!3="lon"
sst_hind&lat=slat
sst_hind&lon=slon
sst_hind&time=time
u_hind&lat=slat
u_hind&lon=slon
u_hind&time=time
v_hind&lat=slat
v_hind&lon=slon
v_hind&time=time
;去趋势化
sst_YYYYMM =cd_calendar(time, -1) 
sst_dpri = dtrend_msg_n(sst_YYYYMM, sst_hind, False, False, 0) 
u_dpri = dtrend_msg_n(sst_YYYYMM, u_hind, False, False, 0) 
v_dpri = dtrend_msg_n(sst_YYYYMM, v_hind, False, False, 0) 
;减去集合平均
sst_hind1 = dim_avg_n(sst_hind,1)
u_hind1 = dim_avg_n(u_hind,1)
v_hind1 = dim_avg_n(v_hind,1)
sst_hind2=new((/5*n*12,10,sizlat,sizlon/), typeof(sst))
u_hind2=new((/5*n*12,10,sizlat,sizlon/), typeof(sst))
v_hind2=new((/5*n*12,10,sizlat,sizlon/), typeof(sst))
do m=0,9,1
    sst_hind2(:,m,:,:)=sst_hind(:,m,:,:)-sst_hind1
    u_hind2(:,m,:,:)=u_hind(:,m,:,:)-u_hind1
    v_hind2(:,m,:,:)=v_hind(:,m,:,:)-v_hind1
end
;reshape以便进行SVD分析
copy_VarCoords(sst_hind, sst_hind2)
copy_VarCoords(u_hind, u_hind2)
copy_VarCoords(v_hind, v_hind2)
sst_hind3=reshape(sst_hind2,(/5*n,12,10,sizlat,sizlon/))
u_hind3=reshape(u_hind2,(/5*n,12,10,sizlat,sizlon/))
v_hind3=reshape(v_hind2,(/5*n,12,10,sizlat,sizlon/))
sst_hind3!0="year"
sst_hind3!1="month"
sst_hind3!2="member"
sst_hind3!3="lat"
sst_hind3!4="lon"
u_hind3!0="year"
u_hind3!1="month"
u_hind3!2="member"
u_hind3!3="lat"
u_hind3!4="lon"
v_hind3!0="year"
v_hind3!1="month"
v_hind3!2="member"
v_hind3!3="lat"
v_hind3!4="lon"
sst_hind4  = sst_hind3(lat|:, lon|:, month|:,member|:,year|:)
u_hind4  = u_hind3(lat|:, lon|:, month|:,member|:,year|:)
v_hind4  = v_hind3(lat|:, lon|:, month|:,member|:,year|:)
sst_pSVD=reshape(sst_hind4,(/sizlat*sizlon*12,10*5*n/))
u_pSVD=reshape(u_hind4,(/sizlat*sizlon*12,10*5*n/))
v_pSVD=reshape(v_hind4,(/sizlat*sizlon*12,10*5*n/))
;缺省值
sst_1d   = sst_pSVD(:,0)
loc_nfv   = ind(.not.ismissing(sst_1d))   ;不是缺省值的位置
num_nfv   = dimsizes(loc_nfv)
sst_SVD = new((/num_nfv,10*5*n/), typeof(sst)) ;non FillValue
sst_SVD= sst_pSVD(loc_nfv,:)     ;;已剔除缺省值,loc_nfv的位置都不是缺省值
u_SVD = new((/num_nfv,10*5*n/), typeof(sst))
v_SVD = new((/num_nfv,10*5*n/), typeof(sst))
u_SVD= u_pSVD(loc_nfv,:) 
v_SVD= v_pSVD(loc_nfv,:) 
uv_SVD = new((/num_nfv*2,10*5*n/), typeof(sst))         
do i=0,:num_nfv-1,1
    uv_SVD(i,:)=u_SVD(i,:)
    uv_SUD(i+num_nfv)=v_SVD(i,:)
end do
;SVD分析
nsvd = 2
homlft = new((/nsvd,num_nfv*2/),float)
hetlft = new((/nsvd,num_nfv*2/),float)
homrgt = new((/nsvd,num_nfv/),float)
hetrgt = new((/nsvd,num_nfv/),float)
X = svdcov(uv_SVD,sst_SVD,nsvd,homlft,hetlft,homrgt,hetrgt)
  ak   = onedtond(X@ak,(/nsvd,10*5*n/))
  bk   = onedtond(X@bk,(/nsvd,10*5*n/))
;获取异性场
sst_helt = new(/sizlat*sizlon*12/,float)
sst_helt = 999   ;记为缺省值
sst_helt(loc_nfv) = hetrgt(0,:)
sst_hmonth1 =reshape(sst_helt,(/sizlat,sizlon,12/))
sst_hmonth = where(abs(sst_hmonth1).eq.999,sst_hmonth1@_FillValue,sst_hmonth1)
sst_hmonth!0="lat"
sst_hmonth!1="lon"
sst_hmonth&lat=slat
sst_hmonth&lon=slon
u_hel = new(num_nfv,float)
v_hel = new(num_nfv,float)
do i=0,num_nfv-1,1
    u_hel(i)=hetlft(0,i)
    v_hel(i)=hetlft(0,i+num_nfv)
end do
u_helt = new(sizlat*sizlon*12,float)
v_helt = new(sizlat*sizlon*12,float)
u_helt = 999   ;记为缺省值
v_helt = 999   ;记为缺省值
u_helt(loc_nfv) = u_hel(:)
v_helt(loc_nfv) = v_hel(:)
u_hmonth1 =reshape(u_helt,(/sizlat,sizlon,12/))
u_hmonth = where(abs(u_hmonth1).eq.999,u_hmonth1@_FillValue,u_hmonth1)
u_hmonth!0="lat"
u_hmonth!1="lon"
u_hmonth&lat=slat
u_hmonth&lon=slon
v_hmonth1 =reshape(v_helt,(/sizlat,sizlon,12/))
v_hmonth = where(abs(v_hmonth1).eq.999,v_hmonth1@_FillValue,v_hmonth1)
v_hmonth!0="lat"
v_hmonth!1="lon"
v_hmonth&lat=slat
v_hmonth&lon=slon
;绘图
wks = gsn_open_wks("png","plot_SVD")
gsn_define_colormap(wks,"GMT_no_green")  
overplot = new(12,graphic)
plot = new(12,graphic)
res                   = True    
res@gsnAddCyclic      = False
res@gsnDraw           = False         
res@gsnFrame          = False        

;;视图参数
res@vpHeightF            = 0.6         ; Changes the aspect ratio
res@vpWidthF             = 0.8
res@vpXF                 = 0.05         ; change start locations
res@vpYF                 = 0.75         ; the plot

;;地图参数
res@mpFillOn        = False             ; 关闭地图填色
; res@mpLandFillColor      = 0            ;land white  
res@mpMinLatF       = latS              ; 地图框的最南纬度
res@mpMaxLatF       = latN              ; 地图框的最北纬度
res@mpCenterLonF    = 190               ; 地图框的中心经度
res@mpMinLonF       = lonL
res@mpMaxLonF       = lonR
res@mpGeophysicalLineThicknessF = 3         ; 地图边界的粗细
; res@pmTickMarkDisplayMode= "Always"         ; 坐标上标签上添加度符号
res@mpShapeMode              = "FreeAspect" ;!!!
res@tmXBLabelFontHeightF = 0.022;0.022
res@tmYLLabelFontHeightF = 0.022;0.022
res@tmXBLabelDeltaF = -0.6
res@tmYLLabelDeltaF = -0.6

; res@cnLevelSelectionMode  = "ExplicitLevels"
; res@cnLevels              = (/-20,20/) ; -20与20均对应着0.05置信度，但前者对应负值异常，后者对应正值异常
; res@cnFillColors          = (/"gray","white","gray"/)

;;等值线的参数
res@cnFillOn              = True        ;开启等值线填色
res@cnLinesOn             = False       ;关闭等值线
res@cnInfoLabelOn         = False       ;关闭等值线信息标签 
res@cnLineLabelsOn        = False       ;关闭等值线数值标签 
res@cnLevelSelectionMode = "ExplicitLevels"
res@cnLevels             = (/-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9/)

;;色标参数
res@lbLabelBarOn          = False        ;开启绘制色标            
; res@pmLabelBarOrthogonalPosF = .20      ;调整色标垂直方向的位置
res@lbBoxEndCapStyle = "TriangleBothEnds"
; res@tmXBLabelFontHeightF = 0.02         ;字体大小
res@tmXBLabelFontThicknessF = 50        ;粗细
; res@tmYLLabelFontHeightF = 0.02         ;字体大小
res@tmXBMajorLengthF     = 0.01         ;主刻度线的长度
res@tmYLMajorLengthF     = 0.01         ;

res@tmXBMinorOn          = False        ;不设次刻度线
res@tmYLMinorOn          = False        ;
res@tmXTOn               = False        ;不设XT轴
res@tmYROn               = False        ;不设YR轴
;--设置坐标标签
res@tmXBMode = "Explicit"
lon_value = ispan(120,280,20)
; lon_lable = lon_value+"E";"~S~o~N~E"
res@tmXBValues = lon_value
res@tmXBLabels = (/"140~S~o~N~E","160~S~o~N~E","180","160~S~o~N~W",\
                  "140~S~o~N~W","120~S~o~N~W","100~S~o~N~W","80~S~o~N~W"/)
lat_value = ispan(0,65,20)
; lat_lable = lat_value+"N";"~S~o~N~N"
res@tmYLValues = lat_value
res@tmYLLabels = (/"0","20~S~o~N~N","40~S~o~N~N","60~S~o~N~N","65~S~o~N~N"/)
;风场
res@vcVectorDrawOrder         = "PreDraw"
res@vcPositionMode            = "ArrowTail"  ;格点的位置对应着箭头的尾部
res@vcGlyphStyle              = "FillArrow"  ;LineArrow;FillArrow;WindBarb;CurlyVector
res@vcLineArrowThicknessF     = 3.0
; resv@vcMinDistanceF            = 0.033      ; 箭头之间的最小距离
; resv@vcMinMagnitudeF           = 1.0       ; 不绘制小于该数值的箭头
;---------------------------------
res@vcLevelSelectionMode    = "ManualLevels"
;res@vcMaxLevelValF      = 18
;res@vcMinLevelValF      = 2
res@vcLevelSpacingF     = 2
res@vcMonoLineArrowColor = False
res@vcLineArrowColor   = "green"
;resv@vcLineArrowEdgeThicknessF = 6.0         ; 箭头边界的粗细
;resv@vcLineArrowEdgeColor      = "white"   ; 箭头边界的颜色
;resv@vcLineArrowFillColor      = "black"   ; 箭头内部填充的颜色
;resv@vcLineArrowWidthF         = 1.0       ; 尾部宽度
;resv@vcLineArrowHeadXF         = 1.0       ; 见P211 施宁实习教程
; resv@vcLineArrowHeadYF         = 1.0      ;
; resv@vcLineArrowHeadInteriorXF = 1.0      ;
;--------------------------------
res@vcFillArrowEdgeThicknessF = 1         ; 箭头边界的粗细
res@vcFillArrowEdgeColor      = "black"   ; 箭头边界的颜色
res@vcFillArrowFillColor      = "black"   ; 箭头内部填充的颜色
res@vcFillArrowWidthF         = 0.05       ; 尾部宽度
; resv@vcFillArrowWidthF         = 0.1       ; 尾部宽度
; resv@vcFillArrowHeadXF         = 0.2       ; 见P211 施宁实习教程
res@vcFillArrowHeadYF         = 0.15     ;
res@vcFillArrowHeadInteriorXF = 0.05      ; 
         
; resv@vcFillArrowMinFracWidthF =1.0 
; resv@vcFillArrowHeadMinFracXF =1.0  
; resv@vcFillArrowHeadMinFracYF =1.0 

;-------------------------------------------------------- 
; *******设定参考箭头****************
;refmag                         = 0.2        ;标准长度的参考箭头所表示的风速大小
res@vcRefAnnoOn               = False
res@vcRefMagnitudeF           = 0.2    ;单位长度的箭头所表示的实际物理量的大小
res@vcRefLengthF              = 0.065
; resv@vcRefLengthF              = 0.07
res@vcRefAnnoBackgroundColor  = "white"
res@vcRefAnnoPerimOn          = False     ;关闭边框                                         
res@vcRefAnnoFontHeightF      = 0.008         

res@vcRefAnnoString1On     = False         ;不绘制参考箭头上方的字符        
res@vcRefAnnoString2On     = False         ;绘制参考箭头下方的字符
res@vcRefAnnoString2       = "0.2 m/s"  
res@vcRefAnnoFontThicknessF = 0.02
res@vcRefAnnoFontHeightF = 0.03       
res@vcRefAnnoJust          = "TopRight"  ;参考箭头放至图右侧
res@vcRefAnnoOrthogonalPosF  = -0.05  ;调整其垂直位置
res@vcRefAnnoParallelPosF    = 0.98   ;调整其水平位置

res@gsnLeftString = " "
res@gsnRightString = " "
resR@vcRefAnnoOn               = True

do i=0,11,1
   plot(i) = gsn_csm_contour_map_ce(wks, sst_hmonth(:,:,i), res) 
   overplot(i)=gsn_csm_vector_map_ce(wks,u_hmonth(:,:,i),v_hmonth(:,:,i),res)
   overlay(plot(i), overplot(i))
end do
resP                  = True                   ; modify the panel plot
resP@gsnPanelCenter       = True         ; 居中对齐每行图（默认）
resP@gsnMaximize         = True         ; 页面最大化（large format）
resP@gsnPanelLabelBar    = True         ; 多附图共用一个色标（add common colorbar）
resP@gsnPanelYWhiteSpacePercent = 13    ; 每两个相邻子图之间Y轴方向的空白大小
resP@gsnPanelXWhiteSpacePercent = 4     ; 每两个相邻子图之间X轴方向的空白大小
;color code
resP@pmLabelBarOrthogonalPosF =0.37     ;色bar位置
resP@pmLabelBarHeightF =0.05  
resP@lbLabelFontHeightF =0.013
resP@gsnPanelFigureStrings= (/"a(FEB)","b(MAR)","c(APR)","d(MAY)","e(JUN)",\
                             "f(JUL)","g(AUG)","h(SEP)","i(OCT)","j(NOV)","k(DEC)","l(JAN)"/) ; add strings to panel
resP@amJust   = "TopLeft"
esP@gsnPanelRowSpec = True 
gsn_panel(wks,plot,(/4,3/),resP)
end

