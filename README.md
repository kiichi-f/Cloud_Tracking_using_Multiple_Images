## Cloud tracking using LIR images with averaging of multiple images

*Japanese follows English*

### Overview

The codes perform cloud tracking using thermal infrared images taken by LIR onboard the Venus orbiter Akatsuki. In this method, noises in the images are reduced by averaging multiple images in the coordinate system rotating with the typical background flow. There are two classes of codes that play the following roles.

- `cloud_track.py` : Compute the cross-correlation surface from the latitude-longitude mapped LIR images
- `summarize.py` : Create csv files for statistical analysis that contain the estimated wind speed and local time of the observation point

This code has been verified on Windows 10, Python 3.6.6, and Anaconda 4.5.11. The versions of each library installed in the virtual environment used are listed in a separate spec-file.txt file. (On Windows, you can create a similar environment by using "conda create --name [virtual environment name] --file spec-file.txt", but it will also install libraries that are not necessary to run this code. For more information about anaconda virtual environments, see [Managing environments](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html))

### Procedure to run cloud_track.py

1. Create a directory `input/` in the directory where the executable code is located.

2. Copy the LIR source file used for cloud tracking into `input/r????`. (???? contains a four-digit orbit number)

   - Supported files are L3c netCDF4 files with latitude/longitude projection. (assuming an image pixel size of 1440 in the longitude direction and 720 in the latitude direction)

   - The file names must remain the same as they were downloaded from the archive. Each file must be placed in a folder with all corresponding orbit numbers. 

   - - Example: `input/r0079/lir_20180410_065612_pic_l3c_v20190401.nc`

3. In `clud_track.py`, replace `orbits = [i for i in range(79,85)]` with a list of orbit numbers to be analyzed. In this example, we use images of orbits 79-84 for cloud tracking. These orbit numbers must be numbers that exist in `input/`.

4. In the directory where the executable code is located, type the command  `python cloud_track.py` to run it. 

5. Once executed, a directory `output/` is created in the directory where the executable code is located, and for each orbit number and each group within the orbit number, csv files containing the cross-correlation surface data (a pickle file containing jpeg images and a two-dimensional array representing the correlation surface), the path to the original data used to create the cross-correlation surface, and the exposure time of the original data are generated. The pickle files are needed when running `summerize.py`. These data are not directly referenced during statistical analysis. 
   - Example: `output/orbit0079/01/ccmap(806.0, 257.0).pickle`


### Procedure to run summarize.py

The code creates csv files for statistical analysis that contain the wind speeds estimated from cross-correlation surfaces and the local times of the observation points. It needs to be run after `cloud_track.py` has been excuted to generate correlation surface data.

1. In `summarize.py`, replace `orbits = [i for i in range(79,85)]` with a list of orbit numbers to be     analyzed. In this example, we use images of orbits 79-84 for cloud tracking. These orbit numbers must be numbers that exist in `output/`.
2. In the directory where the executable code is located, type the command `python summarize.py` to run it. 
3. Once executed, two files for statistical analysis, `summaryCCmap__????.csv` and `.pickle`, are generated. These files are different in format and contain the same information. These files contain values characterizing cross-correlation surfaces and the associated metadata. 

#### Contents of the file for statistical analysis

| **Name**                    | **Explanation**                                              | **Unit**                                                     |
| --------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| (#)                         | Line  number                                                 |                                                              |
| orbit                       | Orbit  number                                                |                                                              |
| groupNum                    | Group  number in the orbit                                   |                                                              |
| u                           | Zonal  velocity that maximizes the cross-correlation coefficient | m/s  : A value obtained by subtracting the assumed background zonal velocity from  the zonal velocity and dividing it by cosθ (θ is latitude). This quantity is  proportional to the angular velocity. |
| v                           | Meridional  velocity velocity that maximizes the cross-correlation coefficient | m/s                                                          |
| lt                          | Average  of the local times of the images used               | hour                                                         |
| lat                         | Average  of the latitudes of the images used                 | degree                                                       |
| lon                         | Average  of the longitudes of the images used                | degree                                                       |
| time_begin                  | The  time of the first image used                            |                                                              |
| distance_mean               | Average  distance between Venus and the spacecraft           | 1e4  km                                                      |
| maxValue                    | Maximum  cross-correlation coefficient                       |                                                              |
| region  x (x=99,95,90)      | The  number of pixels in the cross-correlation surface that have top 100-x %  correlation coefficients |                                                              |
| longestAxes  x (x=99,95,90) | The  major axis of the region in the cross-correlation surface that has top 100-x  % correlation coefficients | pixel/h  (Number of pixels in the cross-correlation surface) |
| path                        | Path  to the cross-correlation data                          |                                                              |




---

## 画像の重ね合わせによるLIR雲追跡の実行方法


### 概要

本研究で用いた、重ね合わせによるLIR雲追跡のコードは2種類あり、それぞれ以下の役割を果たす。
- `cloud_track.py`  緯度経度マッピングされたLIR画像から、相互相関曲面を計算する。(本研究オリジナルの解析手法)
- `summarize.py`  得られた相互相関曲面から推定風速や観測地点のローカルタイムなどをまとめた統計解析用のcsvファイルを作成する(補助)

本コードは `Windows 10`, `Python 3.6.6`, `Anaconda 4.5.11` で実行できることを確認している。用いた仮想環境にインストールした各ライブラリのバージョンは別途`spec-file.txt`にまとめた。(Windowsであれば、`conda create --name (任意の仮想環境名) --file spec-file.txt`とすれば同様の環境を作ることができるが、このコードを動かすために必要でないライブラリも同時にインストールされる。anacondaの仮想環境については[Managing environments](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)を参照)

### cloud_track.py 実行手順

1. 実行コードのあるディレクトリに、ディレクトリ`input/`を作成
2. `input/r????内`に雲追跡に用いるLIRソースファイルをコピーする。ただし、????には4桁の軌道番号が入る。
    - 対応ファイルはLIRの緯度経度展開されたL3cのnetCDF4ファイル（画像のピクセルサイズは経度方向1440、緯度方向720を想定）
    - ファイル名はアーカイブからダウンロードしたままでなければならない。また、各ファイルは全て対応する軌道番号のフォルダ内に置く必要がある。
	    - 例）`input/r0079/lir_20180410_065612_pic_l3c_v20190401.nc`
3. `clud_track.py`内の`orbits = [i for i in range(79,85)]`を解析したい軌道番号が記されたリストに書き換える。この例では軌道番号79〜84の画像を用いて雲追跡を行う。この軌道番号は`input/`内に存在する番号でなければならない。
4. 実行コードのあるディレクトリで`python cloud_track.py`とコマンドを打ち、実行 
5. 実行すると実行コードのあるディレクトリに`output/`が作成され、その中に各軌道番号、軌道番号内グループごとに、相互相関曲面(jpeg形式の画像と曲面を表す2次元配列を格納されたpickleファイル)データと相互相関曲面を作成するときに使った元データのパスや元データの撮影時刻などの情報が記されたcsv形式のファイルが出力される。出力されたpickleファイルは`summerize.py`実行時に必要になる。統計解析時にこれらのデータを直接参照することはない。
    - 例）`output/orbit0079/01/ccmap(806.0, 257.0).pickle`

### summarize.py 実行手順

得られた相互相関曲面群から推定風速や観測地点のローカルタイムなどをまとめた統計解析用のcsvファイルを作成するためのコードであるため、`cloud_track.py`を実行し、必要な相関曲面データが出力された後に実行する必要がある。
1. `summarize.py`内の`orbits = [i for i in range(79,85)]`を解析したい軌道番号が記されたリストに書き換える。この例では軌道番号79〜84の画像を用いて雲追跡を行う。この軌道番号は`output/`内に存在する番号でなければならない。
2. 実行コードのあるディレクトリで`python summarize.py`とコマンドを打ち、実行 
3. 実行すると`output/`以下に、 統計解析用の`summaryCCmap_????.csv`,`.pickle`の2つのファイルが出力される。これらのファイルはどちらも同じ情報が格納された異なるファイル形式である。このファイルに読み込んだ全ての相互相関曲面の特徴量とそれに付随するメタデータが書き込まれている。

#### 統計解析用ファイルの見方
| 名前                       | 説明                                                | 単位など                                    |
|----------------------------|-----------------------------------------------------|---------------------------------------------|
| (#)                        | 行番号                                              |                                             |
| orbit                      | 軌道番号                                            |                                             |
| groupNum                   | 同じ軌道内のグループ番号                            |                                             |
| u                          | 相互相関係数が最大となる東西風速                    | m/s, 対地東西風速から仮定した背景東西風速度を引いた値をcosθ (θは緯度)で割ったもので、角速度に比例する。               |
| v                          | 相互相関係数が最大となる南北風速                    | m/s                                         |
| lt                         | LIR画像における中心のローカルタイムの平均           | hour                                        |
| lat                        | LIR画像における中心の経度の平均                     | degree                                      |
| lon                        | LIR画像における中心の緯度                           | degree                                      |
| time_begin                 | 用いたLIR画像のうち、一番最初に撮影された画像の時刻 |                                             |
| distance_mean              | 撮影時の金星-探査機間の距離の平均                   | 1e4 km                                      |
| maxValue                   | 相互相関係数の最大値                                |                                             |
| region x (x=99,95,90)      | 相関曲面内の上位 100-x %以上の領域の個数            |                                             |
| longestAxes x (x=99,95,90) | 相関曲面内の上位 100-x %以上の領域の長軸の最大値    | pixel/h (相互相関曲面でのピクセル数に換算) |
| path                       | 相互相関曲面2次元データのパス                       |                                             |
