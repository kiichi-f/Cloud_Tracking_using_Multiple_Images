### 1. `git clone` と 仮想環境作成
- github https://github.com/kiichi-f/Cloud_Tracking_using_Multiple_Images
    - 1軌道分のサンプルの入力ファイルも含まれています: `Cloud_Tracking_using_Multiple_Images/input/r0079`
```sh
git clone git@github.com:kiichi-f/Cloud_Tracking_using_Multiple_Images.git
```
- anacondaでライブラリインストール (コマンド例...`spec-file.txt`はバージョンが古いです。)
```sh
conda create -n demo_cloud_track python=3.8 numpy jupyter matplotlib pandas tqdm
conda install -c conda-forge opencv
conda install -c anaconda netcdf4
```

### 2. 雲追跡実行
```sh
python cloud_track.py
```
- githubの入力サンプルのみでも10分くらいかかるかもしれません。
```python
# コンソール出力例
orbit=79
bgnSeq 9.1209796875 16.4234
100%|████████████████████████████████████████████████████████████████████████████| 121/121 [08:22<00:00,  4.15s/it]
```
- 出力ファイルは `Cloud_Tracking_using_Multiple_Images/output/orbit0079/00`
    - 軌道番号`orbit0079`の`00`番目のグループの出力
    - `cmap(806.0, 257.0).pickle`や`cmap(806.0, 257.0).pickle`は個々の雲追跡(重ね合わせ後)の結果
        - pickleは生の情報(2次元配列)でjpgはそれを簡易的に可視化したもの
        - `(806.0, 257.0)`は画像内の座標
    -  `sub(806.0, 257.0).pickle/jpg`は個々の雲追跡(重ね合わせ前)の結果
        - 今後使用はしないが、デバッグ用

### 3. 実行結果をまとめ上げる
```sh
python summarize.py
```
- githubのサンプルだと数秒
```sh
# コンソール出力例
0 2018-04-10 17:01:14 12.78787443576389 276.8680661955423
```
- 出力ファイルは `Cloud_Tracking_using_Multiple_Images/output/characteristicsCCmap_0079.{csv,pickle}`
    - 個々の相互相関マップの特徴量を取り出している。
- 表の見方については[readme](https://github.com/kiichi-f/Cloud_Tracking_using_Multiple_Images#%E7%B5%B1%E8%A8%88%E8%A7%A3%E6%9E%90%E7%94%A8%E3%83%95%E3%82%A1%E3%82%A4%E3%83%AB%E3%81%AE%E8%A6%8B%E6%96%B9)に記載
- orbitの番号を指定して実行(112行目)
```python
orbits = [i for i in range(79,80,1)]
```

### 4. 最終的に用いたデータ
- `vector_for_statistical_analysis_10x10_all.csv`
- 他