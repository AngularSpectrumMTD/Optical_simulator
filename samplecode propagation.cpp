#include"WaveFront.h"
#include"Model.h"

double lambda = 630e-9;//�g��[m]
int WIDTH = 1024;//�摜����
int HEIGHT = 1024;//�摜�c��
double Dx = 1e-6;//�W�{�Ԋu
double Dy = 1e-6;//�W�{�Ԋu

int main()
{
	//openmp��p����O��ł�
	//�����̒P�ʂ� m �ł�

	WaveFront input(256, 256, 2e-6, 2e-6, lambda);
	//������: ���s�N�Z���� �c�s�N�Z���� ���W�{�Ԋu �c�W�{�Ԋu �g��

	double radius = 256 * 2e-6 / 8;//�~�`�J���̔��a or �K�E�X�֐��̔��l��
	double rectw = 256 * 2e-6 / 8;//��`�J���̉�
	double recth = 256 * 2e-6 / 8;//��`�J���̏c
	input.SetCirc(radius);//���g�Ƃ��ĉ~�`�J����ʉ߂������s�����g�p
	//input.SetGaussian(radius, 2);//���g�Ƃ���2���K�E�V�A���r�[�����g�p
	//input.SetRect(rectw, recth);//���g�Ƃ��ċ�`�J����ʉ߂������s�����g�p
	
	input.Normalize();//�摜�o�͗p�ɐ��K��(�p�ӂ��Ă���J���ɂ͕s�v)
	input.SaveBmp("input test field AMPLITUDE.bmp", AMPLITUDE);//�U����
	input.SaveBmp("input test field PHASE.bmp", PHASE);//�ʑ���
	//������: �t�@�C���� �o�͌`��(REAL/IMAGE/PHASE/AMPLITUDE/INTENSITY)

	//input.ModRandomphase();�ʑ��̗����� ����ɂ����g���U������
	
	double distance = 1e-3;//�`������ 1 mm
	input.AsmProp(distance);//�p�X�y�N�g���@��p���Ďw�肵���������R��Ԃŉ�ܓ`��
	//input.ExactAsmProp(distance);//��茵���ȉ������߂�ꍇ(������4�{����)

	input.SaveAsCsv("test field.csv",X_AXIS,input.GetNy()/2);//csv�`���Ō��g��ۑ�
	//������: �t�@�C���� �o�͂ɗp���鎲(X_AXIS/Y_AXIS) �ǂ̈ʒu�̈ꎟ�����z��p���邩

	input.Normalize();//�摜�o�͗p�ɐ��K��
	input.SaveBmp("test field AMPLITUDE.bmp",AMPLITUDE);
	input.SaveBmp("test field PHASE.bmp", PHASE);
}