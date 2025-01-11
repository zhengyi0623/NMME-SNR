begin
latS   =  -10
latN   =  30

lonL   =  150
lonR   =  270

;读取数据
f    = systemfunc("ls /public1/home/scb6369/rzy/CanCM4i/*nc")
n    =dimsizes(f)
f_sst=addfile(f(0),"r")
lat  =f_sst->Y
lon  =f_sst->X
sst  =f_sst->sst
tim =f_sst->S
slon      = lon({lonL:lonR})
sizlon=dimsizes(slon)
slat      = lat({latS:latN})
sizlat=dimsizes(slat)
sst_sell=new((/38*12,10,12,sizlat,sizlon/), typeof(sst))
time=new(38*12,typeof(tim))
do i=0,n-1,1
    f_sst=addfile(f(i),"r")
    sst  =f_sst->sst
    tim =f_sst->S
    sst_pri1  = sst - 273.15       ;转为摄氏度
    copy_VarCoords(sst, sst_pri1)
    do j=0,59,1
        time(j+i*60)=tim(j)
        sst_sell(j+i*60,:,:,:,:) =sst_pri1(j,:,:,{latS:latN}, {lonL:lonR}) ;提取3月开始预测数据
    end do     
end do
delete([/sst,sst_pri1/])
delete(tim)
f_sst=addfile("sst2016_2018.nc","r")
tim =f_sst->S
sst  =f_sst->sst(:,:,:,:,:)
sst_pri1  = sst - 273.15 
copy_VarCoords(sst, sst_pri1)
do j=0,35,1
 time(j+12*35)=tim(j)
 sst_sell(j+12*35,:,:,:,:) =sst_pri1(j,:,:,{latS:latN}, {lonL:lonR}) 
end do
sst_pri = where(abs(sst_sell).eq.273.15,sst_sell@_FillValue,sst_sell)
sst_pri@units = "degC"
sst_pri!0="time"
sst_pri!1="member"
sst_pri!2="lead"
sst_pri!3="lat"
sst_pri!4="lon"
sst_pri&lat=slat
sst_pri&lon=slon
sst_pri&time=time
sst_ave = dim_avg_n(sst_pri,1)
copy_VarCoords(sst_pri(:,0,:,:,:), sst_ave)
nino_34_lead=new((/38*12,12/),float)
do i=0,11
    sst_selc=sst_ave(:,i,:,:)
    copy_VarCoords(sst_pri(:,0,0,:,:), sst_selc)
    sstClm   = clmMonTLL(sst_selc)               ;monthly avg(monthly climatology)   
    ssta = calcMonAnomTLL (sst_selc, sstClm)    ; replace with anonamlies;Compute the anomalies from the climatology just computed
    copy_VarCoords(sst_selc, ssta)
    rad    = 4.0*atan(1.0)/180.0
    clat   = sqrt( cos(rad*slat) )              ; gw for gaussian gri
    wsst  = ssta                               ; copy meta data
    wsst  = ssta*conform(wsst, clat, 1)
    nino_341 =  dim_avg_n_Wrap(wsst(:,{-5:5},{190:240}), (/1,2/))  ;nino34
    sst_YYYYMM =cd_calendar(time, -1)
    copy_VarCoords(sst_pri(:,0,0,0,0), nino_341)
    nino_34 = dtrend_msg_n(sst_YYYYMM, nino_341, False, False, 0) 
    nino_34_lead(:,i)=nino_34
end do
path_out = "NMME-nino34.nc"
system("rm -f " + path_out)
ncdf = addfile(path_out, "c")
ncdf->nino_34 = nino_34_lead ;储存pc1 命名为VMI
end