#include"..\include\WaveFront.h"

int main()
{
	//openmpを用いる前提です
	//長さの単位は m です

	WaveFront input(256, 256, 2e-6, 2e-6, 630e-9);
	//↑引数: 横ピクセル数 縦ピクセル数 横標本間隔 縦標本間隔 波長

	double radius = 256 * 2e-6 / 8;//円形開口の半径 or ガウス関数の半値幅
	double rectw = 256 * 2e-6 / 8;//矩形開口の横
	double recth = 256 * 2e-6 / 8;//矩形開口の縦
	input.GenerateCirc(radius);//光波として円形開口を通過した平行光を使用
	//input.GenerateGaussian(radius, 2);//光波として2次ガウシアンビームを使用
	//input.GenerateRect(rectw, recth);//光波として矩形開口を通過した平行光を使用
	
	input.Normalize();//画像出力用に正規化(用意している開口には不要)
	input.SaveBmp("input test field AMPLITUDE.bmp", AMPLITUDE);//振幅像
	input.SaveBmp("input test field PHASE.bmp", PHASE);//位相像
	//↑引数: ファイル名 出力形式(REAL/IMAGE/PHASE/AMPLITUDE/INTENSITY)//

	//input.ModRandomphase();位相の乱数化 これにより光波が散乱する
	
	double distance = 1e-3;//伝搬距離 1 mm
	input.AsmProp(distance);//角スペクトル法を用いて指定した距離自由空間で回折伝搬
	//input.ExactAsmProp(distance);//より厳密な解を求める場合(メモリ4倍消費)

	input.SaveAsCsv("test field.csv",X_AXIS,input.GetNy()/2);//csv形式で光波を保存
	//↑引数: ファイル名 出力に用いる軸(X_AXIS/Y_AXIS) どの位置の一次元分布を用いるか

	input.Normalize();//画像出力用に正規化
	input.SaveBmp("test field AMPLITUDE.bmp",AMPLITUDE);
	input.SaveBmp("test field PHASE.bmp", PHASE);
}