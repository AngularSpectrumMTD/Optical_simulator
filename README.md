# Optical_simulator
本プログラムはヘルムホルツ方程式にのっとった
等方的かつ均一な媒質中の波動の回折伝搬をシミュレートするものです．

伝搬計算手法としては角スペクトル法を用いています．
この手法と，開口を設定する関数/レンズの位相を生成する関数を用いることで
種々のシミュレーションを行えます．
こちらは非平行平面間の回折伝搬にも対応しています．
更に通常の角スペクトル法で生じるエイリアシング軽減策・消滅策として
帯域制限角スペクトル法および4倍拡張各スペクトル法を採用しています．

以上のような基本的な光学的シミュレーションから，
ポリゴン法を用いた計算機合成ホログラムを作成するための光波を生成できます．
物体光波の結像再生シミュレーションにも対応しています．

また，マルチスレッド対応をしています．
＊画像のロード・出力に外部ヘッダオンリライブラリを使用します．
詳しくはreferenceをご覧ください．

サンプルコードである
samplecode propagation.cpp
は平行平面間の伝搬計算をシミュレーションする用のコードであり，
samplecode propagation between non parallel planes.cpp
は非平行な平面間の伝搬計算をシミュレーションする用のコードです．
また，
samplecode model calculation.cpp
はmqo(メタセコイア)のモデルをロードし，物体光波を計算するためのものです．
24bitBMPのテクスチャには対応しています．現在はモデル・テクスチャ共にslnと同一階層にあることが前提です．

計算例

・平行平面間の回折伝搬計算

![prop](https://user-images.githubusercontent.com/65929274/91724081-b9726780-ebd7-11ea-8cdf-a6ec5f6bd4ef.gif)


・非平行平面間の回折伝搬計算

![prop2](https://user-images.githubusercontent.com/65929274/91724222-f2124100-ebd7-11ea-8b38-07849d783e93.gif)


・物体光波の回折伝搬計算

![物体光波](https://user-images.githubusercontent.com/65929274/91642074-56f85a80-ea63-11ea-8d8d-5ae9f5232313.png)


・隠面消去

![図6](https://user-images.githubusercontent.com/65929274/91662330-be261580-eb1c-11ea-89c9-521075bcdaeb.png)

