import cv2
import numpy    as np
import matplotlib.pyplot as plt
import pickle
import pandas as pd
import glob
from tqdm import tqdm

def plotFitedEllipse(ccmap,detected_ellipses):

    #変数名シノニムを作るのではなく新たにメモリを確保
    ccmap_add_fitted_ellipse =  np.copy(ccmap)

    ellipseText = ' '
    for ellipse in detected_ellipses:
        cv2.ellipse(ccmap_add_fitted_ellipse, (ellipse['center'], ellipse['axes'], ellipse['angle']),-1)
        ellipseText += str(ellipse['threshold']) + 'long axis {:.3g} '.format(ellipse['axes'][1]/1.637)

    fs =22
    fs2=24

    fig, ax = plt.subplots(1, 1, figsize=(10,8))
    mappable = ax.pcolormesh(ccmap_add_fitted_ellipse,cmap='jet',vmin=np.min(ccmap),vmax=np.max(ccmap))

    upper = lambda parcent : np.sort(ccmap.reshape((1,np.array(ccmap).size)))[0][int(parcent/100.*np.array(ccmap).size)]
    ax.contour(ccmap, levels = [upper(90),upper(95),upper(99)],colors=['black','green','blue'])

    tick = [-0.2,0.0,0.2,0.4,0.6]
    # tick = [-0.1,0.0,0.1,0.2,0.3]


    cbar = fig.colorbar(mappable,ticks=tick)
    cbar.set_label('correlation coefficient',fontsize=fs2)
    cbar.ax.set_yticklabels(tick ,fontsize=fs)

    # ax.set_xticks(np.arange(72+1.637*40,72-1.637*40,-1.637*10))
    # ax.set_xticklabels(np.arange(40,-41,-10),fontsize=fs)
    # ax.set_xlabel(r'Eastward velocity / $\cos(\theta)$  (m/s)',fontsize=fs2)
    #
    # ax.set_yticks(np.arange(72+1.637*40,72-1.637*40,-1.637*10))
    # ax.set_yticklabels(np.arange(40,-41,-10),fontsize=fs)
    # ax.set_ylabel('Northward velocity (m/s)',fontsize=fs2)

    #最大値をマーキング
    maxPiY,maxPiX = np.argmax(ccmap)//ccmap.shape[1],np.argmax(ccmap)%ccmap.shape[1]
    ax.scatter(maxPiX,maxPiY,marker='x',c='black',s=300)

    U,V = (maxPiX-72)/1.637,(maxPiY-72)/1.637

    info_text = 'u {:.3g} '.format( U ) + 'v {:.3g} '.format( V ) + ellipseText


    plt.gcf().text(0,0,info_text)

    plt.tight_layout()
    plt.show()
    return info_text


def recoedFitedEllipse(ccmap,detected_ellipses=[],threshold = 99):
    upper_ = lambda ccmap,parcent : np.sort(ccmap.reshape((1,np.array(ccmap).size)))[0][int(parcent/100.*np.array(ccmap).size)]
     #輪郭検出のための前処理：2値化
    _, binaryMap = cv2.threshold(ccmap,upper_(ccmap,threshold),1,cv2.THRESH_BINARY)
    _, contours, _ = cv2.findContours(binaryMap.astype('u1'),cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)

    #楕円近似
    for conter in contours:
        #There should be at least 5 points to fit the ellipse in function 'fitEllipse'
        if conter.shape[0]>5:
            center,axes,angle = cv2.fitEllipse(conter)

        #     中心がマップ内にあるもののみを取り出す
            # if min(center[0],center[1])>0 and center[0]<ccmap.shape[0] and center[1]<ccmap.shape[1]:

            ellipseInfo = {'center':center,
                                'axes':axes,
                                'angle':angle,
                                'threshold':threshold}
            detected_ellipses += [ellipseInfo]

    return detected_ellipses

#最大となるUV
def maxUV_SR(ccmap):
    maxPiY,maxPiX = np.argmax(ccmap)//ccmap.shape[1],np.argmax(ccmap)%ccmap.shape[1]
    #2*7*12+1 == 84+1+84
    U,V = (maxPiX-85)/1.637,-(maxPiY-85)/1.637
    return U,V

if __name__ == '__main__':

    # orbit番号で指定。グループ番号はfor文で回す
    #orbits = [101]
    # orbits = [99,100,101,102]
    orbits = [i for i in range(79,80,1)]
    row = []
    for orbit in orbits:
        # target_dir = '/Users/kiichi/Desktop/orbit'+'{0:04d}'.format(orbit)
        target_dir = './output/orbit'+"{0:04d}".format(orbit)+'/'

        for subDir in sorted(glob.glob(target_dir+'/??')):
            if len(glob.glob(subDir+'/refarencesImagesSummary_*.pickle'))>0:
                pathSummary = glob.glob(subDir+'/refarencesImagesSummary_*.pickle')[0]

                with open(pathSummary,'rb') as f:
                    dfSummary = pickle.load(f)
                '''
                subsolar_longtitudeが地理座標系になっているので変更する
                '''
                _, groupNum, time_begin, _, distance_mean, _, _, subsolar_SRlongtitude, _, _, _ = np.array(dfSummary)[0]
                print(groupNum, time_begin, distance_mean, subsolar_SRlongtitude)

                for pathCCmap in sorted(glob.glob(subDir+'/ccmap(*).pickle')):
                    a,b,c = pathCCmap.rfind('('),pathCCmap.rfind(','),pathCCmap.rfind(')')
                    cX,cY = float(pathCCmap[a+1:b]),float(pathCCmap[b+1:c])

                    cLON,cLAT = cX/1440 * 360, (360-cY)/360 * 90
                    cLT = (cLON - subsolar_SRlongtitude)%360 / 15

                    with open(pathCCmap,'rb') as f:
                        ccmap = pickle.load(f)

                    d99 = recoedFitedEllipse(ccmap,detected_ellipses=[],threshold = 99)
                    d95 = recoedFitedEllipse(ccmap,detected_ellipses=[],threshold = 95)
                    d90 = recoedFitedEllipse(ccmap,detected_ellipses=[],threshold = 90)

                    u,v = maxUV_SR(ccmap)

                    if len(d99)>0:
                        longestAxes99 = max(list(map(lambda di :di['axes'][1], d99)))
                    else:
                        longestAxes99 = 0
            #             print(pathCCmap)

                    if len(d95)>0:
                        longestAxes95 = max(list(map(lambda di :di['axes'][1], d95)))
                    else :
                        longestAxes95 = 0

                    if len(d90)>0:
                        longestAxes90 = max(list(map(lambda di :di['axes'][1], d90)))
                    else :
                        longestAxes90 = 0


                    row += [ [orbit ,groupNum, u, v,
                                   cLT,cLAT,cLON,
                                   time_begin, distance_mean,
                                   np.max(ccmap),len(d99),len(d95),len(d90),
                                   longestAxes99,longestAxes95,longestAxes90,
                                   pathCCmap] ]

        #         d99 = recoedFitedEllipse(ccmap[:,:-n],detected_ellipses=[],threshold = 99)
        #         d95 = recoedFitedEllipse(ccmap[:,:-n],detected_ellipses=d99,threshold = 95)
        #         detected_ellipses = recoedFitedEllipse(ccmap[:,:-n],detected_ellipses=d95,threshold = 90)
        #         plotFitedEllipse(ccmap[:,:-n],d99)

    header = ['orbit','groupNum','u','v',
            'lt','lat','lon','time_begin','distance_mean',
            'maxValue','region99','region95','region90',
            'longestAxes99','longestAxes95','longestAxes90',
            'path']

    df = pd.DataFrame(row,columns=header)
    df.to_csv(target_dir+'/../characteristicsCCmap_'+'{0:04d}'.format(orbit)+'.csv')
    with open(target_dir+'/../characteristicsCCmap_'+'{0:04d}'.format(orbit)+'.pickle','wb') as f:
        pickle.dump(df,f,protocol=2)

        #         print(cY,cLAT,cLT,u,v)
        #         print(d99,d95,d90)#を使って特徴量抽出を行う
        # /Users/kiichi/Desktop/orbit0101/00/refarencesImagesSummary_orbit0101_mimGroupSize6_groupNum0.csv
