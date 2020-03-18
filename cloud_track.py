#!/usr/bin/env python
# coding: utf-8

import glob
import netCDF4
import numpy    as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import cv2
import os
import pandas as pd
import pickle
import itertools
import datetime
from tqdm import tqdm

def my_makedirs(path):
    """
    ディレクトリを作成する

    Parameters
    ----------
    path : str
        作成したいディレクトリのパス
    """
    if not os.path.isdir(path):
        os.makedirs(path)

def getObservationTime(path):
    """
    ファイル名から観測時刻を取得する
    ファイル名の例：lir_20180410_065012_pic_l3c_v20190401.nc

    Parameters
    ----------
    path : str
        対象となるnetCDF4ファイルのパス

    Returns
    -------
        '%Y%m%d_%H%M%S'形式の時刻
    """
    return datetime.datetime.strptime(path[-36:-21],'%Y%m%d_%H%M%S')

def mean_angle(angles):
    """
    複数の角度の平均を計算する

    Parameters
    ----------
    angles : list
        度数法で記された角度のリスト

    Returns
    -------
        度数法で記された角度の平均
    """
    vec = lambda a : np.cos(2*np.pi*a/360)+np.sin(2*np.pi*a/360)*1j
    m_angle = np.angle(np.array(list(map(vec,angles))).mean())*180/np.pi
    if m_angle>0:
        return m_angle
    else:
        return m_angle+360

def high_pass_filtering(img,size=81):
    """
    cv2.GaussianBlurを用いたハイパス処理

    Parameters
    ----------
    img : ndarray
        緯度経度展開された輝度温度マップ
    size : 2*int+1
        取り出したい構造の大きさを決めるパラメーター
        デフォルトではsize=81であり、1440x720の輝度温度マップでは81/720*180~20度以下の構造を強調する。

    Returns
    -------
        ハイパス処理された輝度温度マップ
    """
    kernel_size = (size,size)
    highPassFilter = lambda img : img - cv2.GaussianBlur(img, kernel_size,0)
    mask           = lambda img : ma.masked_greater(ma.masked_less(img,-10),10)
    return mask(highPassFilter(img))

def getSequencesGroup(ncFileList,mimGroupSize = 7,intervaltime = 3600,buffer_rate_min = 0.8,buffer_rate_max = 1.2):
    """
    時間的に連続した輝度温度マップを特定の枚数、時間間隔ごとにグルーピングする。
    1グループ内の画像は、画像を重ね合わせ、雲追跡を行って各地点で1つの風速を求めるための一連の画像群になる。

    Parameters
    ----------
    ncFileList : list
        対象の輝度温度マップがあるパスのリスト
    mimGroupSize : int
        最小のグループサイズ
    intervaltime : int (second)
        グループ内の時間間隔（の目安）デフォルトでは1時間となっており、1時間を大きく超える撮影間隔の連続した2枚の画像は、同じグループに属さない。
        1時間を大きく下回る撮影間隔の2枚の画像があったら、後に撮影された画像は用いない。
    buffer_rate_min : double (0~1)
        intervaltime*buffer_rate_minがグループ内連続した画像の最小許容時間間隔となる。
    buffer_rate_max : double (>1)
        intervaltime*buffer_rate_maxがグループ内連続した画像の最大許容時間間隔となる。
    Returns
    -------
    seqencesGroup : list of list
        ncFileListをグループ分けしたもの。内側のリストには画像のパスではなく、ncFileListにおける画像の番号が格納されている。 len(seqencesGroup)=雲追跡を行う画像グループ数
    """
    listObservationTime = list( map( getObservationTime, ncFileList))

    obsInterval = np.array(listObservationTime)-np.roll(listObservationTime,1) #N番目：(N-1)とNの間隔、0番目：最初と最後の間隔

    validIntervalPair = []
    for i,interval in  enumerate(obsInterval):
        if (interval > datetime.timedelta(seconds=intervaltime*buffer_rate_min) ) & (interval < datetime.timedelta(seconds=intervaltime*buffer_rate_max) ):
            validIntervalPair += [(i-1,i)]

    seqences = []
    seqencesGroup = []
    for i in range(1,len(validIntervalPair)):
        if validIntervalPair[i-1][1]==validIntervalPair[i][0]:
            seqences += [validIntervalPair[i][0]]
        #連続画像が途切れたら
        elif len(seqences) >= 2*mimGroupSize :
            seqencesGroup += [ [seqences[0]-1] + seqences[:mimGroupSize-1] ]
            seqencesGroup += [ seqences[-mimGroupSize:] ]
            seqences = []
        elif len(seqences) >= mimGroupSize :
            seqencesGroup += [ [seqences[0]-1] + seqences]
            seqences = []
        elif seqences != [] :
            seqences = []

    return seqencesGroup


# In[3]:

def getAveragedObservationData(ncFileList,sequences,plot=False):
    """
    複数枚の画像を背景東西風と共に動く座標系で重ね合わせる。ncFileListとsequencesで指定した全ての画像を用いるのではなく、
    ncFileListとsequencesで指定した用いる可能性のある画像のうち金星により近い連続した12枚の画像から連続した4枚の画像を重ね合わせ、
    計9枚の重ねた画像および撮影時刻や距離の平均を計算する。

    TODO：背景東西風の緯度分布を外部パラメーターにする。

    Parameters
    ----------
    ncFileList : list
        対象の輝度温度マップがあるパスのリスト
    sequences　: int
        ncFileListのうち、重ね合わせに用いる可能性のある画像の番号
    plot : bool
        Trueならば重ね合わせた画像を描画

    Returns
    -------
    averaged_observations : list of dict
        平均後の画像および撮影時刻や距離の平均などを格納したdictのリスト

    """
    # u_angle  = np.array([ 80.        ,  80.01218624,  80.04876354,  80.10978768,
    #         80.19535185,  80.305587  ,  80.44066237,  80.60078604,
    #         80.7862058 ,  80.99721006,  81.23412895,  81.4973356 ,
    #         81.78724759,  82.10432862,  82.44909035,  82.82209443,
    #         83.22395487,  83.65534052,  84.11697794,  84.60965449,
    #         85.1342218 ,  85.69159949,  86.28277941,  86.90883019,
    #         87.57090228,  88.27023352,  89.00815524,  89.78609901,
    #         90.60560406,  91.46832543,  92.37604307,  93.33067178,
    #         94.33427227,  95.38906343,  96.49743588,  97.6619671 ,
    #         98.8854382 , 100.17085265, 101.52145721, 102.94076527,
    #        104.43258315, 106.00103947, 107.65061837, 109.38619689,
    #        111.21308728, 113.13708499, 115.16452317, 117.30233485,
    #        119.55812399, 121.94024694, 124.45790615, 127.12125833,
    #        129.72131412, 132.04455686, 134.09497807, 135.87635643,
    #        137.39226342, 138.6460685 , 139.64094411, 140.37987034,
    #        140.86563936, 141.10085957, 141.08795951, 140.82919153,
    #        140.32663524, 139.58220073, 138.59763153, 137.37450746,
    #        135.91424717, 134.21811056, 132.28720094, 130.12246708,
    #        127.72470501, 125.09455969, 122.23252645, 119.13895232,
    #        115.81403716, 112.25783462, 108.47025295, 104.45105564,
    #        100.19986192,  95.71614708,  90.99924262,  86.04833631,
    #         80.86247202,  75.44054945,  69.78132368,  63.88340453,
    #         57.7452559 ,  51.36519476])

    u_angle  = np.array([ 95. for i in  range(90)])

    wind_speeds = np.array(list(reversed([ u_angle[int(i/4)] for i in range(360) ])) + [ u_angle[int(i/4)] for i in range(360) ])
    initial_time = getObservationTime(ncFileList[sequences[0]])
    #約1時間ごとに取得された12枚の観測を記録
    observations = []
    SR_subsolar_longtitudes = []

    beginDistance = netCDF4.Dataset(ncFileList[sequences[0]], 'r').variables['S_DISTAV'][0].data/1e4
    endDistance    = netCDF4.Dataset(ncFileList[sequences[-1]], 'r').variables['S_DISTAV'][0].data/1e4

    if beginDistance<endDistance :
        seq = sequences[:12]
        print('bgnSeq',beginDistance,endDistance)
    else :
        seq = sequences[-12:]
        print('endSeq',beginDistance,endDistance)

    for num in seq:

        ncData = netCDF4.Dataset(ncFileList[num], 'r')

        #番地が小さい方が南なので、y方向は反転して、北が上になるように読み込み
        btemp = ncData.variables['btemp'][0]
        if btemp.shape == (1440,2880) :
            img = high_pass_filtering(btemp[::-2,::2],size=65)
        elif btemp.shape == (720,1440) :
            img = high_pass_filtering(btemp[::-1,:],size=65)

        time = getObservationTime(ncFileList[num])
        rolled_img = img

        distance = ncData.variables['S_DISTAV'][0].data/1e4 #1e4km
        subsolar_longtitude = ncData.variables['S_SOLLON'][0].data
        akatsuki_longtitude = ncData.variables['S_SSCLON'][0].data
        akatsuki_localtime = ncData.variables['S_SSCLT'][0].data


        interval_second = (time-initial_time).days*24*3600+(time-initial_time).seconds
    #     print( displacement_pixcel(interval_second) )
#         rolled_img = np.roll(img,displacement_pixcel(interval_second), axis=1)

        displacement_pixcel = lambda U_background,interval_second : int(U_background / 100. * 360 / (2 * np.pi * 6050 / 1440) * interval_second/3600.)
        rolled_img = np.array([np.roll(strip,displacement_pixcel(U,interval_second)) for strip,U in zip (img,wind_speeds)])
        mask           = lambda img : ma.masked_greater(ma.masked_less(img,-10),10)

    #SRにおける個々のsubsolar_longtitude
        SR_subsolar_longtitudes += [(subsolar_longtitude+displacement_pixcel(95,interval_second)*360/1440.)%360]
        # print(SR_subsolar_longtitudes[-1],subsolar_longtitude)

        observation = {'img' : mask(rolled_img),
                               'distance' : distance,
                               'time' : time,
                               'subsolar_longtitude': subsolar_longtitude,
                               'akatsuki_longtitude': akatsuki_longtitude,
                               'akatsuki_localtime': akatsuki_localtime }

        observations += [observation]

    averaged = lambda key,j :sum([observations[i][key] for i in range(0+j,4+j)])/4
    averaged_observations = []

    #SR座標系におけるsubsolar_longtitudeも計算するが、個々の重ね合わせ画像の平均ではなく、12枚のsequences全部の平均を計算する


    subsolar_SRlongtitude = mean_angle(SR_subsolar_longtitudes)

    for j in range(9):
        averaged_observation = {'img' : averaged('img',j),
                               'distance' : averaged('distance',j),
                               'time' : observations[j+2]['time'],
                               'subsolar_SRlongtitude':subsolar_SRlongtitude,
                               #地理座標における3枚目の画像の経度
                               'subsolar_longtitude':  observations[j+2]['subsolar_longtitude'],
                               'akatsuki_longtitude':  observations[j+2]['akatsuki_longtitude'],
                               'akatsuki_localtime':  observations[j+2]['akatsuki_localtime'] }

        if plot==True:
            plotObservstion(averaged_observation)

        averaged_observations += [averaged_observation]
    return averaged_observations

# def plotObservstion(observation):
#
#     img  = observation['img']
#     time = observation['time']
#     distance = observation['distance']
#     subsolar_longtitude = observation['subsolar_longtitude']
#     akatsuki_localtime = observation['akatsuki_localtime']
#     akatsuki_longtitude = observation['akatsuki_longtitude']
#
#     fig, ax = plt.subplots(1, 1,
#                            figsize=(8,6))
#     mappable = ax.imshow(img,cmap='inferno',vmin=-1,vmax=1)
#     ax.set_title(time.strftime('%Y-%m-%d-%H:%M:%S'),fontsize=18)
#     ax.set_ylabel('latitude',fontsize=18)
#
#     plt.gcf().text(0,0,
#                    'distance {:.3g} 1e4km '.format(distance) +
#                    'subsolar_longtitude {:.3g} '.format(subsolar_longtitude) +
#                    'akatsuki_longtitude {:.3g} '.format(akatsuki_longtitude) +
#                    'akatsuki_localtime {:.3g} '.format(akatsuki_localtime) )
#
#     cbar=fig.colorbar(mappable,shrink=0.5)
#     cbar.set_label('Brightness Temperature (K)',size=18)
#     fig.tight_layout()
#     return


# In[4]:


#画像が端で切れているとright<leftと逆転して返す
def imageLocation(img):
    """
    画像の緯度、経度方向の端を見つける（テンプレートエリア、サーチエリアを設定する際に用いる）
    緯度経度マップに展開された時、写っていない緯度経度領域はmaskされていることを用いている。

    Parameters
    ----------
    img : ndarray
        緯度経度展開された輝度温度マップ

    Returns
    -------
    left,right,top,bottom : int,int,int,int
        画像の左(西)、右(東)、上(北)、下(南)端のピクセル番地
    """

#     print(np.arange(start=0,stop=img[360].size)[img[360].mask==False])
    if img[360].mask[0]==True: #画像が端で切れていない
        try:
            left = np.arange(start=0,stop=img[360].size)[img[360].mask==False][0]
            right = np.arange(start=0,stop=img[360].size)[img[360].mask==False][-1]
            top = np.arange(start=0,stop=img.T[int((left+right)/2)].size)[img.T[int((left+right)/2)].mask==False][0]
            bottom = np.arange(start=0,stop=img.T[int((left+right)/2)].size)[img.T[int((left+right)/2)].mask==False][-1]
        except:
            left = 0#np.arange(start=0,stop=img[360].size)[img[360].mask==False][0]
            right = 1#np.arange(start=0,stop=img[360].size)[img[360].mask==False][-1]
            top = 0
            bottom = 1
        return left,right,top,bottom

    else: #画像が端で切れている
        _img = img
        img =  np.roll(_img,720, axis=1)
        try:
            left = np.arange(start=0,stop=img[360].size)[img[360].mask==False][0]
            right = np.arange(start=0,stop=img[360].size)[img[360].mask==False][-1]
            top = np.arange(start=0,stop=img.T[int((left+right)/2)].size)[img.T[int((left+right)/2)].mask==False][0]
            bottom = np.arange(start=0,stop=img.T[int((left+right)/2)].size)[img.T[int((left+right)/2)].mask==False][-1]
        except:
            left = 0#np.arange(start=0,stop=img[360].size)[img[360].mask==False][0]
            right = 1#np.arange(start=0,stop=img[360].size)[img[360].mask==False][-1]
            top = 0
            bottom = 1
        return left+720,right-720,top,bottom

def sizeOfSearchArea_for_SRcoordinate(delta_hour=None,size_templateX=None,size_templateY=None,
                                                                max_displacementXperHour=7,max_displacementYperHour=7):
    max_displacementX = max_displacementXperHour*delta_hour
    max_displacementY = max_displacementYperHour*delta_hour
    size_serchX = 2*max_displacementX+size_templateX
    size_serchY = 2*max_displacementY+size_templateY
    return (size_serchX,size_serchY)

def listTemplateCenter(img,max_delta_hour=None,size_templateX=None,size_templateY=None,
                       max_displacementXperHour=7,max_displacementYperHour=7):

    """
    - 東西、南北風速場測定範囲：-50 ~ +50 m/s
        - pixcel換算で `-7 *delta_hour`~`+7 *delta_hour`

    - テンプレートはもっと極域に伸ばせるかもしれない
    """

    left,right,top,bottom = imageLocation(img)
    size_serchX,size_serchY = sizeOfSearchArea_for_SRcoordinate(max_delta_hour,size_templateX,size_templateY)
#     print( size_serchX,size_serchY)
#     まずは緯度方向
    c00y= top + size_serchY/2
    cY = []
    N = 0
    while c00y + size_templateY/2 * N + size_serchY/2 < bottom :
        cY += [c00y + size_templateY/2 * N]
        N += 1

    #画像が端にまたがっていない場合
    if left<right:
        c00x= left + size_serchX/2
        cX = []
        N = 0
        while c00x + size_templateX/2 * N + size_serchX/2  < right :
            cX += [c00x + size_templateX/2 * N]
            N += 1

    #画像が端にまたがっている場合
    else:
        c00x= 0 + size_serchX/2
        cX = []
        N = 0
        while c00x + size_templateX/2 * N + size_serchX/2 < right:
            cX += [c00x + size_templateX/2 * N]
            N += 1

        c00x= 1440 - size_serchX/2
        N = 0
        while c00x - size_templateX/2 * N - size_serchX/2 > left:
            cX += [c00x - size_templateX/2 * N]
            N += 1
    return list(itertools.product(cX, cY))

def plotCenters(img,centers):
    plt.scatter(np.array(centers).T[0],np.array(centers).T[1])
    plt.imshow(img)
    plt.show()
    return


# In[5]:


#相関係数を求める
def subCCmapMatrix(observations,templateCenter,num_observations,size_templateX=None,size_templateY=None):
    cX,cY = templateCenter
    #subMapsMatrixを作成
    subMapList =  [ [ [] for i in range(5)] for j in range(5)]

    for delta in range(4,10):
        for num1 in range(9-delta):
            num2 = num1 + delta
    #         print(num2)
            #早い時刻に撮影されたもの
            templateObservation = observations[num1]['img']
            #遅い時刻に撮影されたもの
            searchObservation = observations[num2]['img']

            size_searchX,size_searchY = sizeOfSearchArea_for_SRcoordinate(delta_hour=delta,
                                                                         size_templateX=size_templateX,
                                                                         size_templateY=size_templateY)

            templateImg = observations[num1]['img'][int(cY-size_templateY/2):int(cY+size_templateY/2),
                                                                          int(cX-size_templateX/2):int(cX+size_templateX/2)]

            #右端はテンプレートの右端に合わせる：書き方が汚いので変えたい
            searchImg = observations[num2]['img'][int(cY-size_searchY/2):int(cY+size_searchY/2),
                                                                       int(cX-size_searchX/2):int(cX+size_searchX/2)]

            """
            何もしなければcc_submap 00 には displacements_x[0]=もっとも右方向にシフト（u最小）、
            displacements_y[0]=もっとも下方向にシフトしたもの（v最小）が入る
            → u軸は左から右(正)、v軸は上から下（逆）
            →最初からu,vの軸を正しいものに置き換えておく
            """
    #         cc_submap = cv2.matchTemplate(searchImg, templateImg, cv2.TM_CCOEFF_NORMED)

            cc_submap = cv2.resize(cv2.matchTemplate(searchImg, templateImg, cv2.TM_CCOEFF_NORMED),
                                           (2*7*12+1,2*7*12+1))[::-1,:]

            subMapList[num1][num2-4] += [cc_submap]



    for delta in range(4,10):
        arr = np.zeros_like(subMapList[0][0][0])
        for num1 in range(9-delta):
            num2 = num1 + delta
            arr += subMapList[num1][num2-4][0]
        if delta == 4:
            subMapList[2][0] += [arr/(9-delta)]
        if delta == 5:
            subMapList[3][1] += [arr/(9-delta)]
        if delta == 6:
            subMapList[4][2] += [arr/(9-delta)]
        if delta == 7:
            subMapList[3][0] += [arr/(9-delta)]
        if delta == 8:
            subMapList[4][1] += [arr/(9-delta)]

    arr = np.zeros_like(subMapList[0][0][0])
    for (i,j) in [(2,0),(3,1),(4,2),(3,0),(4,1)]:
        arr += subMapList[i][j][0]

    ccmap = arr/5
    subMapList[4][0] += [ccmap]

    # with open(savePath('subCCmapMatrix')+'.pickle', 'wb') as f:
    #         pickle.dump(subMapList, f, protocol=2)
    with open(savePath('ccmap')+'.pickle', 'wb') as f:
            pickle.dump(ccmap, f, protocol=2)

    return subMapList,ccmap


# In[ ]:





# In[ ]:





# In[6]:


#描画
def plotSubCCmapMatrix(subMapList,savePath=None):
    fs = 12
    fs2=14

    num = len(subMapList)

    fig, axs = plt.subplots(num, num, figsize=(8,8))

    for num1 in range(num):
        for num2 in range(num):
            ax = axs[num1,num2]
            ax.set_xticks([])
            ax.set_yticks([])
            if subMapList[num1][num2] != []:
                    mappable = ax.pcolormesh(subMapList[num1][num2][0],cmap='jet',vmin=-0.3,vmax=0.3)
    plt.tight_layout()
#     plt.show()
    plt.savefig(savePath)

    plt.clf()
    plt.cla()
    plt.close()

    return

def plotCCmap(ccmap,savePath=None):
    fs = 22
    fs2= 24

    fig, ax = plt.subplots(1, 1, figsize=(10,8))
    mappable = ax.pcolormesh(ccmap,cmap='jet')#jet')#,vmin=-0.5,vmax=0.5)

    upper = lambda parcent : np.sort(ccmap.reshape((1,np.array(ccmap).size)))[0][int(parcent/100.*np.array(ccmap).size)]
    try :
        ax.contour(ccmap, levels = [upper(90),upper(95),upper(99)],colors=['black','green','blue'])
    except :
        print('contour error')


    #tick = np.linspace(-1,1,41)

    cbar = fig.colorbar(mappable)#,ticks=tick)
    # cbar.set_label('correlation coefficient',fontsize=fs2)
    # cbar.ax.set_yticklabels(map('{:.1g}'.format,tick) ,fontsize=fs)

    ax.set_xticks(np.arange(7*12+1.637*45,7*12-1.637*46,-1.637*15))
    ax.set_xticklabels(np.arange(45,-46,-15),fontsize=fs)
    ax.set_xlabel(r'Eastward velocity / $\mathrm{\cos  \theta}$  (m/s)',fontsize=fs2)

    ax.set_yticks(np.arange(7*12+1.637*45,7*12-1.637*46,-1.637*15))
    ax.set_yticklabels(np.arange(45,-46,-15),fontsize=fs)
    ax.set_ylabel('Northward velocity (m/s)',fontsize=fs2)

    #最大値をマーキング
    maxPiY,maxPiX = np.argmax(ccmap)//ccmap.shape[1],np.argmax(ccmap)%ccmap.shape[1]
    ax.scatter(maxPiX,maxPiY,marker='x',c='black',s=300)

    plt.tight_layout()
    plt.savefig(savePath)
#     plt.show()

    plt.clf()
    plt.cla()
    plt.close()
    return


# In[7]:


def saveCCmapInfp(orbit,num_observations,groupNum,row,tempdir):
    header = ['savePath','center','SorcePaths']
    name = 'ccmapInfo_orbit'+'{0:04d}'.format(orbit)+'_mimGroupSize'+str(num_observations)+'_groupNum'+str(groupNum)
    df = pd.DataFrame(row,columns=header)
    df.to_csv(tempdir+'{0:02d}'.format(groupNum)+'/'+name+'.csv')

    with open(tempdir+'{0:02d}'.format(groupNum)+'/'+name+'.pickle', 'wb') as f:
            pickle.dump(df , f, protocol=2)
    return

#重ね合わせに使った画像の撮影状況のサマリー
def summerizeReferencedImages(orbit,observations,num_observations):
    distances,times,ss_SRlons,ss_lons,ak_lons,ak_lt = np.array(list(map(lambda di :list(di.values())[1:], observations[:num_observations]))).T

    row = [ [ orbit,
              groupNum,
              np.min(times),
              np.max(times),
              np.mean(distances),
              np.max(distances),
              np.min(distances),
              mean_angle(ss_SRlons),
              mean_angle(ss_lons),
              mean_angle(ak_lons),
              np.mean(ak_lt)] ]

    header = ['orbit','groupNum','time_begin','time_end','distance_mean','distance_max','distance_min','subsolar_SRlongtitude','subsolar_longtitude','akatsuki_longtitude','akatsuki_localtime']

    name = 'refarencesImagesSummary_orbit'+str(orbit)+'_mimGroupSize'+str(num_observations)+'_groupNum'+str(groupNum)
    df = pd.DataFrame(row,columns=header)
    df.to_csv(tempdir+'{0:02d}'.format(groupNum)+'/'+name+'.csv')

    with open(tempdir+'{0:02d}'.format(groupNum)+'/'+name+'.pickle', 'wb') as f:
            pickle.dump(df , f, protocol=2)
    return


# In[8]:


orbits = [i for i in range(79,80,1)]
for orbit in orbits :
    print('orbit='+str(orbit))
    # path = '/Users/kiichi/Google Drive/akatsuki/191030toMac/data/r'+"{0:04d}".format(orbit)+'/*.nc'
    path = 'input/r'+"{0:04d}".format(orbit)+'/*.nc'
    # tempdir =  '/Users/kiichi/Desktop/orbit'+"{0:04d}".format(orbit)+'/'
    tempdir = './output/orbit'+"{0:04d}".format(orbit)+'/'
    my_makedirs(tempdir)
    ncFileList = sorted(glob.glob(path))

    size_templateX=4*20
    size_templateY=4*20
    #一回の相関マップ重ね合わせに使う画像の枚数
    num_observations = 12

    #特定の時間間隔で連続して撮影されている画像番号を出力
    sequencesGroup=getSequencesGroup(ncFileList,mimGroupSize = num_observations,intervaltime = 3600,buffer_rate_min = 0.8,buffer_rate_max = 1.2)
#     print(sequencesGroup)

    for groupNum,selectedGroup in enumerate([sequencesGroup[0]]):
    # for groupNum,selectedGroup in enumerate(sequencesGroup):
        observations = getAveragedObservationData(ncFileList,selectedGroup,plot=False)

        templateCenters = listTemplateCenter(observations[0]['img'],
                                             max_delta_hour=8, #重ね合わせているので最大時間幅は12-2*2時間
                                             size_templateX=size_templateX,
                                             size_templateY=size_templateY)

        my_makedirs(tempdir+'{0:02d}'.format(groupNum))

        #付属情報保存用
        row = []
        for templateCenter in tqdm(templateCenters):
            savePath = lambda flag : tempdir+'{0:02d}'.format(groupNum)+'/'+flag+str(templateCenter)

            subMapList,ccmap = subCCmapMatrix(observations,templateCenter,num_observations,
                                                         size_templateX=size_templateX,
                                                         size_templateY=size_templateY)


            plotSubCCmapMatrix(subMapList,savePath=savePath('sub')+'.jpg' )
            plotCCmap(ccmap,savePath = savePath('cc')+'.jpg')

            #CSVで付属情報の保存
            paths = [ncFileList[index] for index in [selectedGroup[i] for i in range(num_observations) ] ]
            row += [ [savePath('cc')+'.jpg',
                      templateCenter,
                      paths] ]

            # break

        saveCCmapInfp(orbit,num_observations,groupNum,row,tempdir)
        summerizeReferencedImages(orbit,observations,num_observations)

        # break


# In[ ]:





# In[ ]:
