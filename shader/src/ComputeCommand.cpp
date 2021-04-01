#include "../include/ComputeFilterApp.h"

int ComputeFilterApp::ExecuteFFTCommand()
{
	ExecuteCopyCommand(m_RfullsizeTex.at(0), m_RWfullsizeTex.at(0));
	ExecuteCopyCommand(m_RfullsizeTex.at(1), m_RWfullsizeTex.at(1));
	ExecuteClearCommand(m_RWfullsizeTex.at(1));//このパスではいらない
	ExecuteFFTCommand(m_RWfullsizeTex.at(0), m_RWfullsizeTex.at(1));
	ExecuteCalcurateAmplitudeCommand(
		m_RWfullsizeTex.at(0), m_RWfullsizeTex.at(1)
		, m_RWfullsizeTex.at(2), m_RWfullsizeTex.at(3));//あえて振幅に

	ExecuteOrdinalNormalizeCommand(m_RWfullsizeTex.at(2));

	m_commandList->SetComputeRootDescriptorTable(1, m_RWfullsizeTex.at(2).ReadState());
	m_commandList->SetComputeRootDescriptorTable(3, m_RWfullsizeTex.at(3).WriteState());
	m_commandList->SetPipelineState(m_pipelines["raiseCS"].Get());
	m_commandList->Dispatch(1, m_texheight, 1);

	return 3;
}

int ComputeFilterApp::ExecuteIFFTCommand()
{
	ExecuteCopyCommand(m_RfullsizeTex.at(0), m_RWfullsizeTex.at(0));
	ExecuteCopyCommand(m_RfullsizeTex.at(1), m_RWfullsizeTex.at(1));
	ExecuteFFTCommand(m_RWfullsizeTex.at(0), m_RWfullsizeTex.at(1));
	ExecuteIFFTCommand(m_RWfullsizeTex.at(0), m_RWfullsizeTex.at(1));

	ExecuteOrdinalNormalizeCommand(m_RWfullsizeTex.at(0));

	return 0;
}

int ComputeFilterApp::ExecuteConvolutionCommand()
{
	ExecuteCopyCommand(m_RfullsizeTex.at(0), m_RWfullsizeTex.at(0));
	ExecuteCopyCommand(m_RfullsizeTex.at(1), m_RWfullsizeTex.at(4));

	ExecuteClearCommand(m_RWfullsizeTex.at(1));
	ExecuteClearCommand(m_RWfullsizeTex.at(3));

	ExecuteCalcMaxMinCommand(m_RWfullsizeTex.at(4), m_RWmaxminTex.at(0), m_RWmaxminTex.at(1));//片方のカーネルについて情報収集 今は4
	ExecuteDivideMaxAmpCommand(
		m_RWmaxminTex.at(0), m_RWmaxminTex.at(1)
		, m_RWfullsizeTex.at(4), m_RWfullsizeTex.at(1)
		, m_RWfullsizeTex.at(5), m_RWfullsizeTex.at(6));//そのカーネルの振幅を1以下にする

	ExecuteColorScaleCommand(m_RWfullsizeTex.at(5), m_RWfullsizeTex.at(2));//スケーリングで振幅をより適切な値に変換

	ExecuteConvolutionCommand(
		m_RWfullsizeTex.at(0), m_RWfullsizeTex.at(1)
		, m_RWfullsizeTex.at(2), m_RWfullsizeTex.at(3)
		, m_RWfullsizeTex.at(4), m_RWfullsizeTex.at(5));

	ExecuteCalcurateAmplitudeCommand(
		m_RWfullsizeTex.at(4), m_RWfullsizeTex.at(5)
		, m_RWfullsizeTex.at(6), m_RWfullsizeTex.at(7));

	ExecuteCalcMaxMinCommand(m_RWfullsizeTex.at(6), m_RWmaxminTex.at(0), m_RWmaxminTex.at(1));
	ExecuteDivideMaxAmpCommand(
		m_RWmaxminTex.at(0), m_RWmaxminTex.at(1)
		, m_RWfullsizeTex.at(6), m_RWfullsizeTex.at(7)
		, m_RWfullsizeTex.at(4), m_RWfullsizeTex.at(5));

	return 4;
}

int ComputeFilterApp::ExecuteASMCommand()
{
	ComputeParameters cp;
	cp.intervalx = m_interval_x * 1e-6;
	cp.intervaly = cp.intervalx * (m_texwidth * 1.0f) / (m_texheight * 1.0f);
	cp.distance0 = m_propdistance[0] / 1000.0f;
	cp.gausssigma = m_gausssigma;
	cp.gausssigma2 = m_gausssigma2;
	cp.glareintensity = m_glareintensity;
	cp.threshold = m_threshold;
	cp.focus0 = m_focus[0] / 1000.0f;
	cp.N = m_N;
	cp.r = m_r;

	WriteToUploadHeapMemory(m_mainComputeCB[0].Get(), sizeof(cp), &cp);

	//ExecuteCopyCommand(m_RfullsizeTex.at(0), m_RWfullsizeTex.at(0));
	//ExecuteCopyCommand(m_RfullsizeTex.at(0), m_RWfullsizeTex.at(1));
	ExecuteClearCommand(m_RWfullsizeTex.at(1));
	ExecuteCopyCommand(m_RfullsizeTex.at(1), m_RWfullsizeTex.at(2));

	switch (m_inputfrag)
	{
	case InputFrag_Image:
		ExecuteCopyCommand(m_RfullsizeTex.at(0), m_RWfullsizeTex.at(0));
		break;
	case InputFrag_Aperture:
		ExecuteDrawPolygonCommand(m_RWfullsizeTex.at(0));
		break;
	}

	if (m_focusfrag == FocusFrag_ON && m_focus[0] != 0 && m_propdistance[0] != 0)
	{
		ExecuteDrawQuadraticPhaseCommand(m_RWfullsizeTex.at(6), m_RWfullsizeTex.at(7));

		ExecuteMultiplyCommand(
			m_RWfullsizeTex.at(0), m_RWfullsizeTex.at(1)
			, m_RWfullsizeTex.at(6), m_RWfullsizeTex.at(7)
			, m_RWfullsizeTex.at(4), m_RWfullsizeTex.at(5));

		ExecuteClearCommand(m_RWfullsizeTex.at(0));
		ExecuteClearCommand(m_RWfullsizeTex.at(1));
		ExecuteClearCommand(m_RWfullsizeTex.at(6));
		ExecuteClearCommand(m_RWfullsizeTex.at(7));

		ExecuteCopyCommand(m_RWfullsizeTex.at(4), m_RWfullsizeTex.at(0));
		ExecuteCopyCommand(m_RWfullsizeTex.at(5), m_RWfullsizeTex.at(1));
	}

	ExecuteFFTCommand(m_RWfullsizeTex.at(0), m_RWfullsizeTex.at(1));

	ExecuteDrawFRFCommand(m_RWfullsizeTex.at(6), m_RWfullsizeTex.at(7));
	ExecuteMultiplyCommand(
		m_RWfullsizeTex.at(0), m_RWfullsizeTex.at(1)
		, m_RWfullsizeTex.at(6), m_RWfullsizeTex.at(7)
		, m_RWfullsizeTex.at(4), m_RWfullsizeTex.at(5));
	ExecuteIFFTCommand(m_RWfullsizeTex.at(4), m_RWfullsizeTex.at(5));

	ExecuteCalcurateAmplitudeCommand(
		m_RWfullsizeTex.at(4), m_RWfullsizeTex.at(5)
		, m_RWfullsizeTex.at(6), m_RWfullsizeTex.at(7));
	ExecuteCalcMaxMinCommand(m_RWfullsizeTex.at(6), m_RWmaxminTex.at(0), m_RWmaxminTex.at(1));
	ExecuteDivideMaxAmpCommand(
		m_RWmaxminTex.at(0), m_RWmaxminTex.at(1)
		, m_RWfullsizeTex.at(4), m_RWfullsizeTex.at(5)
		, m_RWfullsizeTex.at(2), m_RWfullsizeTex.at(3));

	if (m_focusfrag == FocusFrag_ON && m_focus[0] != 0)
	{
		ExecuteRaiseRICommand(
			m_RWfullsizeTex.at(2), m_RWfullsizeTex.at(3)
			, m_RWfullsizeTex.at(6), m_RWfullsizeTex.at(7));//グレアの輝度を底上げ
	}
	else
	{
		ExecuteCopyCommand(m_RWfullsizeTex.at(2), m_RWfullsizeTex.at(6));
		ExecuteCopyCommand(m_RWfullsizeTex.at(3), m_RWfullsizeTex.at(7));
	}

	ExecuteConvertToFormatCommand(m_RWfullsizeTex.at(6), m_RWfullsizeTex.at(7), m_RWfullsizeTex.at(4), m_RWfullsizeTex.at(5));
	ExecuteNormalizeWaveFrontCommand(m_RWfullsizeTex.at(4), m_RWfullsizeTex.at(5), m_RWfullsizeTex.at(0), m_RWfullsizeTex.at(2));

	return 0;
}

int ComputeFilterApp::ExecuteKamehamehaCommand()
{
	ComputeParameters cp;
	cp.intervalx = m_interval_x * 1e-6;
	cp.intervaly = cp.intervalx * (m_texwidth * 1.0f) / (m_texheight * 1.0f);
	cp.distance0 = m_propdistance[0] / 1000;
	cp.gausssigma = m_gausssigma;
	cp.glareintensity = m_glareintensity;
	cp.threshold = m_threshold;
	cp.focus0 = m_focus[0] / 1000.0f;
	cp.N = m_N;
	cp.r = m_r;
	cp.N1 = m_N1;
	cp.r1 = m_r1;
	cp.posx = m_posx;
	cp.posy = m_posy;

	WriteToUploadHeapMemory(m_mainComputeCB[0].Get(), sizeof(cp), &cp);

	ExecuteCopyCommand(m_RfullsizeTex.at(0), m_RWfullsizeTex.at(0));
	ExecuteCopyCommand(m_RfullsizeTex.at(0), m_RWfullsizeTex.at(1));

	ExecuteClearCommand(m_RWfullsizeTex.at(1));

	ExecuteBinaryThresholdCommand(m_RWfullsizeTex.at(0), m_RWfullsizeTex.at(2));

	//絞り開口
	ExecuteDrawMovingPolygonCommand(m_RWfullsizeTex.at(3));
	ExecuteCopyCommand(m_RWfullsizeTex.at(3), m_RWdisplayTex.at(1));
	ExecuteClearCommand(m_RWfullsizeTex.at(4));

	ExecuteMultiplyCommand(
		m_RWfullsizeTex.at(2), m_RWfullsizeTex.at(1)
		, m_RWfullsizeTex.at(3), m_RWfullsizeTex.at(4)
		, m_RWfullsizeTex.at(5), m_RWfullsizeTex.at(6));

	ExecuteCopyCommand(m_RWfullsizeTex.at(5), m_RWfullsizeTex.at(2));
	ExecuteCopyCommand(m_RWfullsizeTex.at(6), m_RWfullsizeTex.at(1));

	if (m_focusfrag == FocusFrag_ON && m_focus[0] != 0 && m_propdistance[0] != 0)
	{
		ExecuteDrawQuadraticPhaseCommand(m_RWfullsizeTex.at(6), m_RWfullsizeTex.at(7));

		ExecuteMultiplyCommand(
			m_RWfullsizeTex.at(2), m_RWfullsizeTex.at(1)
			, m_RWfullsizeTex.at(6), m_RWfullsizeTex.at(7)
			, m_RWfullsizeTex.at(4), m_RWfullsizeTex.at(5));

		ExecuteClearCommand(m_RWfullsizeTex.at(0));
		ExecuteClearCommand(m_RWfullsizeTex.at(1));
		ExecuteClearCommand(m_RWfullsizeTex.at(6));
		ExecuteClearCommand(m_RWfullsizeTex.at(7));

		ExecuteCopyCommand(m_RWfullsizeTex.at(4), m_RWfullsizeTex.at(0));
		ExecuteCopyCommand(m_RWfullsizeTex.at(5), m_RWfullsizeTex.at(2));
	}

	ExecuteClearCommand(m_RWfullsizeTex.at(3));

	ExecuteFFTCommand(m_RWfullsizeTex.at(2), m_RWfullsizeTex.at(3));

	ExecuteDrawFRFCommand(m_RWfullsizeTex.at(6), m_RWfullsizeTex.at(7));
	ExecuteMultiplyCommand(
		m_RWfullsizeTex.at(2), m_RWfullsizeTex.at(3)
		, m_RWfullsizeTex.at(6), m_RWfullsizeTex.at(7)
		, m_RWfullsizeTex.at(4), m_RWfullsizeTex.at(5));
	ExecuteIFFTCommand(m_RWfullsizeTex.at(4), m_RWfullsizeTex.at(5));

	ExecuteCalcurateAmplitudeCommand(
		m_RWfullsizeTex.at(4), m_RWfullsizeTex.at(5)
		, m_RWfullsizeTex.at(6), m_RWfullsizeTex.at(7));
	ExecuteCalcMaxMinCommand(m_RWfullsizeTex.at(6), m_RWmaxminTex.at(0), m_RWmaxminTex.at(1));
	ExecuteDivideMaxAmpCommand(
		m_RWmaxminTex.at(0), m_RWmaxminTex.at(1)
		, m_RWfullsizeTex.at(6), m_RWfullsizeTex.at(7)
		, m_RWfullsizeTex.at(4), m_RWfullsizeTex.at(5));

	ExecuteCopyCommand(m_RfullsizeTex.at(0), m_RWfullsizeTex.at(0));
	ExecuteAddCommand(m_RWfullsizeTex.at(0), m_RWfullsizeTex.at(4), m_RWfullsizeTex.at(2));

	//退避
	ExecuteCopyCommand(m_RWfullsizeTex.at(4), m_RWdisplayTex.at(0));

	return 2;
}

int ComputeFilterApp::ExecuteOtherCommand()
{
	switch (m_mode)
	{
	case Mode_Sepia:
		m_commandList->SetPipelineState(m_pipelines["sepiaCS"].Get());
		break;
	case Mode_Sobel:
		m_commandList->SetPipelineState(m_pipelines["sobelCS"].Get());
		break;
	case Mode_swap:
		m_commandList->SetPipelineState(m_pipelines["swapCS"].Get());
		break;
	default:
		m_commandList->SetPipelineState(m_pipelines["sepiaCS"].Get());
		break;
	}

	//WriteToUploadHeapMemory(m_mainComputeCB[0].Get(), sizeof(cp), &cp);

	m_commandList->SetComputeRootDescriptorTable(
		1, m_RfullsizeTex.at(0).handleRead
	);
	m_commandList->SetComputeRootDescriptorTable(
		2, m_RfullsizeTex.at(1).handleRead
	);
	m_commandList->SetComputeRootDescriptorTable(
		3, m_RWfullsizeTex.at(0).WriteState()
	);
	m_commandList->SetComputeRootDescriptorTable(
		4, m_RWfullsizeTex.at(1).WriteState()
	);

	m_commandList->Dispatch(m_texwidth, m_texheight, 1);

	return 0;
}

int ComputeFilterApp::ExecuteGaussianBlurCommand()
{
	ComputeParameters cp;
	cp.intervalx = m_interval_x * 1e-6;
	cp.intervaly = cp.intervalx * (m_texwidth * 1.0f) / (m_texheight * 1.0f);
	cp.distance0 = m_propdistance[0];
	cp.gausssigma = m_gausssigma;
	cp.glareintensity = m_glareintensity;
	cp.threshold = m_threshold;

	WriteToUploadHeapMemory(m_mainComputeCB[0].Get(), sizeof(cp), &cp);

	ExecuteCopyCommand(m_RfullsizeTex.at(0), m_RWfullsizeTex.at(0));
	ExecuteCopyCommand(m_RfullsizeTex.at(0), m_RWfullsizeTex.at(1));

	ExecuteBinaryThresholdCommand(m_RWfullsizeTex.at(0), m_RWfullsizeTex.at(2));
	ExecuteClearCommand(m_RWfullsizeTex.at(3));

	ExecuteDrawGaussianCommand(m_RWfullsizeTex.at(6), m_RWfullsizeTex.at(7));

	ExecuteConvolutionCommand(
		m_RWfullsizeTex.at(6), m_RWfullsizeTex.at(7)
		, m_RWfullsizeTex.at(2), m_RWfullsizeTex.at(3)
		, m_RWfullsizeTex.at(0), m_RWfullsizeTex.at(1));
	ExecuteCalcurateAmplitudeCommand(
		m_RWfullsizeTex.at(0), m_RWfullsizeTex.at(1)
		, m_RWfullsizeTex.at(2), m_RWfullsizeTex.at(3));

	ExecuteCalcMaxMinCommand(m_RWfullsizeTex.at(2), m_RWmaxminTex.at(0), m_RWmaxminTex.at(1));
	ExecuteDivideMaxAmpCommand(
		m_RWmaxminTex.at(0), m_RWmaxminTex.at(1)
		, m_RWfullsizeTex.at(2), m_RWfullsizeTex.at(3)
		, m_RWfullsizeTex.at(4), m_RWfullsizeTex.at(5));

	ExecuteCopyCommand(m_RfullsizeTex.at(0), m_RWfullsizeTex.at(0));
	ExecuteAddCommand(m_RWfullsizeTex.at(0), m_RWfullsizeTex.at(4), m_RWfullsizeTex.at(1));

	//退避
	ExecuteCopyCommand(m_RWfullsizeTex.at(4), m_RWdisplayTex.at(0));

	return 1;//もともと
	//return 4;
}

int ComputeFilterApp::ExecuteGlareCommand()
{
	ComputeParameters cp;
	cp.intervalx = 1e-6;
	cp.intervaly = cp.intervalx * (m_texwidth * 1.0f) / (m_texheight * 1.0f);
	cp.distance0 = m_propdistance[0];
	cp.gausssigma = m_gausssigma;
	cp.glareintensity = m_glareintensity;
	cp.threshold = m_threshold;
	cp.N = m_N;
	cp.r = m_r;
	cp.glarelambdasamplenum = m_glarelambdasamplenum;

	WriteToUploadHeapMemory(m_mainComputeCB[0].Get(), sizeof(cp), &cp);

	ExecuteCopyCommand(m_RfullsizeTex.at(0), m_RWfullsizeTex.at(0));//グレアがかかる画像
	ExecuteBinaryThresholdCommand(m_RWfullsizeTex.at(0), m_RWfullsizeTex.at(1));
	ExecuteClearCommand(m_RWfullsizeTex.at(2));

	//グレア用画像
	switch (m_glaremode)
	{
	case GlareMode_EyeRush:
		ExecuteCopyCommand(m_RfullsizeTex.at(1), m_RWfullsizeTex.at(3));
		break;
	case GlareMode_Cross:
		ExecuteDrawPolygonCommand(m_RWfullsizeTex.at(3));

		//ごみ画像
		ExecuteCopyCommand(m_RfullsizeTex.at(3), m_RWfullsizeTex.at(4));

		//ExecuteClearCommand(m_RWfullsizeTex.at(2));
		ExecuteClearCommand(m_RWfullsizeTex.at(5));

		ExecuteMultiplyCommand(
			m_RWfullsizeTex.at(4), m_RWfullsizeTex.at(5)
			, m_RWfullsizeTex.at(3), m_RWfullsizeTex.at(2)
			, m_RWfullsizeTex.at(0), m_RWfullsizeTex.at(6));

		if (m_aperturedust == ApertureDust_ON)
		{
			ExecuteCopyCommand(m_RWfullsizeTex.at(0), m_RWdisplayTex.at(2));
			ExecuteCopyCommand(m_RWfullsizeTex.at(0), m_RWfullsizeTex.at(3));
		}
		else
		{
			ExecuteCopyCommand(m_RWfullsizeTex.at(3), m_RWdisplayTex.at(2));
		}

		break;
	case GlareMode_Halo:
		ExecuteCopyCommand(m_RfullsizeTex.at(2), m_RWfullsizeTex.at(3));
		break;
	default:
		ExecuteCopyCommand(m_RfullsizeTex.at(1), m_RWfullsizeTex.at(3));
		break;
	}

	ExecuteClearCommand(m_RWfullsizeTex.at(4));

	ExecuteFFTCommand(m_RWfullsizeTex.at(3), m_RWfullsizeTex.at(4));//グレア生成

	ExecuteCalcurateAmplitudeCommand(
		m_RWfullsizeTex.at(3), m_RWfullsizeTex.at(4)
		, m_RWfullsizeTex.at(5), m_RWfullsizeTex.at(6));

	//この画像をずらしまくって加算したい(もはやシェーダでずらす)

	ExecuteSpectrumScalingCommand(
		m_RWfullsizeTex.at(5), m_RWfullsizeTex.at(6)
		, m_RWfullsizeTex.at(3), m_RWfullsizeTex.at(4));//波長スケーリング

	ExecuteCopyCommand(m_RWfullsizeTex.at(3), m_RWfullsizeTex.at(5));
	ExecuteCopyCommand(m_RWfullsizeTex.at(4), m_RWfullsizeTex.at(6));

	//正規化部分
	//ExecuteCalcMaxMinCommand(m_RWfullsizeTex.at(5), m_RWmaxminTex.at(0), m_RWmaxminTex.at(1));
	//ExecuteDivideMaxAmpCommand(
	//    m_RWmaxminTex.at(0), m_RWmaxminTex.at(1)
	//    , m_RWfullsizeTex.at(5), m_RWfullsizeTex.at(6)
	//    , m_RWfullsizeTex.at(3), m_RWfullsizeTex.at(4));
	//正規化部分

	 //トーンマッピング正規化部分
	ExecuteToneMappingCommand(m_RWfullsizeTex.at(5), m_RWfullsizeTex.at(3));
	//トーンマッピング正規化部分

	ExecuteRaiseRICommand(
		m_RWfullsizeTex.at(3), m_RWfullsizeTex.at(4)
		, m_RWfullsizeTex.at(5), m_RWfullsizeTex.at(6));//グレアの輝度を底上げ

	ExecuteCopyCommand(m_RWfullsizeTex.at(5), m_RWfullsizeTex.at(3));
	ExecuteCopyCommand(m_RWfullsizeTex.at(6), m_RWfullsizeTex.at(4));

	//退避
	ExecuteCopyCommand(m_RWfullsizeTex.at(3), m_RWdisplayTex.at(1));
	//return 3;//ここで返すとグレア単体
	//二値化画像とグレア画像を畳み込む
	ExecuteConvolutionCommand(
		m_RWfullsizeTex.at(1), m_RWfullsizeTex.at(2)
		, m_RWfullsizeTex.at(3), m_RWfullsizeTex.at(4)
		, m_RWfullsizeTex.at(5), m_RWfullsizeTex.at(6));

	ExecuteCalcurateAmplitudeCommand(
		m_RWfullsizeTex.at(5), m_RWfullsizeTex.at(6)
		, m_RWfullsizeTex.at(0), m_RWfullsizeTex.at(1));

	ExecuteCalcMaxMinCommand(m_RWfullsizeTex.at(0), m_RWmaxminTex.at(0), m_RWmaxminTex.at(1));
	ExecuteDivideMaxAmpCommand(
		m_RWmaxminTex.at(0), m_RWmaxminTex.at(1)
		, m_RWfullsizeTex.at(0), m_RWfullsizeTex.at(1)
		, m_RWfullsizeTex.at(3), m_RWfullsizeTex.at(4));

	ExecuteCopyCommand(m_RfullsizeTex.at(0), m_RWfullsizeTex.at(0));
	ExecuteAddCommand(m_RWfullsizeTex.at(0), m_RWfullsizeTex.at(3), m_RWfullsizeTex.at(1));

	//退避
	ExecuteCopyCommand(m_RWfullsizeTex.at(3), m_RWdisplayTex.at(0));

	return 1;
}

int ComputeFilterApp::ExecuteGodRayCommand()
{
	ComputeParameters cp;
	cp.threshold = m_threshold;
	cp.raylength = m_raylength;
	cp.gausssigma = m_gausssigma;
	cp.posx = m_posx;
	cp.posy = m_posy;

	WriteToUploadHeapMemory(m_mainComputeCB[0].Get(), sizeof(cp), &cp);

	ExecuteCopyCommand(m_RfullsizeTex.at(0), m_RWfullsizeTex.at(0));
	ExecuteCopyCommand(m_RfullsizeTex.at(0), m_RWfullsizeTex.at(1));

	ExecuteClearCommand(m_RWfullsizeTex.at(1));

	ExecuteBinaryThresholdCommand(m_RWfullsizeTex.at(0), m_RWfullsizeTex.at(2));

	for (int i = 0; i < m_godrayquarity; ++i)
	{
		ExecuteRadialBlurCommand(m_RWfullsizeTex.at(2), m_RWfullsizeTex.at(3));
		ExecuteCopyCommand(m_RWfullsizeTex.at(3), m_RWfullsizeTex.at(2));
	}

	ExecuteDrawGaussianCommand(m_RWfullsizeTex.at(6), m_RWfullsizeTex.at(7));

	ExecuteClearCommand(m_RWfullsizeTex.at(2));

	ExecuteConvolutionCommand(
		m_RWfullsizeTex.at(6), m_RWfullsizeTex.at(7)
		, m_RWfullsizeTex.at(3), m_RWfullsizeTex.at(2)
		, m_RWfullsizeTex.at(0), m_RWfullsizeTex.at(1));
	ExecuteCalcurateAmplitudeCommand(
		m_RWfullsizeTex.at(0), m_RWfullsizeTex.at(1)
		, m_RWfullsizeTex.at(2), m_RWfullsizeTex.at(3));

	ExecuteCalcMaxMinCommand(m_RWfullsizeTex.at(2), m_RWmaxminTex.at(0), m_RWmaxminTex.at(1));
	ExecuteDivideMaxAmpCommand(
		m_RWmaxminTex.at(0), m_RWmaxminTex.at(1)
		, m_RWfullsizeTex.at(2), m_RWfullsizeTex.at(3)
		, m_RWfullsizeTex.at(4), m_RWfullsizeTex.at(5));

	ExecuteCopyCommand(m_RfullsizeTex.at(0), m_RWfullsizeTex.at(0));

	ExecuteAddCommand(m_RWfullsizeTex.at(0), m_RWfullsizeTex.at(4), m_RWfullsizeTex.at(1));

	//退避
	ExecuteCopyCommand(m_RWfullsizeTex.at(4), m_RWdisplayTex.at(0));

	return 1;
}

void ComputeFilterApp::generateBurst()
{
	PIXBeginEvent(m_commandList.Get(), 0, "GenerateBurst");

	ExecuteDrawFixPolygonCommand(m_RWfullsizeTex.at(1));

	//ごみ画像
	ExecuteCopyCommand(m_RfullsizeTex.at(3), m_RWfullsizeTex.at(3));

	ExecuteClearCommand(m_RWfullsizeTex.at(2));
	ExecuteClearCommand(m_RWfullsizeTex.at(4));

	ExecuteMultiplyCommand(
		m_RWfullsizeTex.at(1), m_RWfullsizeTex.at(2)
		, m_RWfullsizeTex.at(3), m_RWfullsizeTex.at(4)
		, m_RWfullsizeTex.at(5), m_RWfullsizeTex.at(6));

	ExecuteCopyCommand(m_RWfullsizeTex.at(5), m_RWdisplayTex.at(2));

	ExecuteClearCommand(m_RWfullsizeTex.at(6));

	ExecuteFFTCommand(m_RWfullsizeTex.at(5), m_RWfullsizeTex.at(6));//グレア生成

	ExecuteCalcurateAmplitudeCommand(
		m_RWfullsizeTex.at(5), m_RWfullsizeTex.at(6)
		, m_RWfullsizeTex.at(0), m_RWfullsizeTex.at(1));

	ExecuteSpectrumScalingCommand(
		m_RWfullsizeTex.at(0), m_RWfullsizeTex.at(1)
		, m_RWfullsizeTex.at(2), m_RWfullsizeTex.at(3));//波長スケーリング

	ExecuteToneMappingCommand(m_RWfullsizeTex.at(2), m_RWfullsizeTex.at(3));

	ExecuteRaiseRICommand(
		m_RWfullsizeTex.at(3), m_RWfullsizeTex.at(4)
		, m_RWfullsizeTex.at(5), m_RWfullsizeTex.at(6));

	ExecuteClearCommand(m_RWfullsizeTex.at(1));

	ExecuteClearCommand(m_RWfullsizeTex.at(6));

	/*ExecuteCalcurateAmplitudeCommand(
		m_RWfullsizeTex.at(5), m_RWfullsizeTex.at(6)
		, m_BurstCachedTex, m_RWfullsizeTex.at(1));*/

	//ExecutefadeByGaussCommand(m_RWfullsizeTex.at(0), m_BurstCachedTex);//キャッシュする

	ExecuteCalcurateAmplitudeCommand(
		m_RWfullsizeTex.at(5), m_RWfullsizeTex.at(6)
		, m_RWfullsizeTex.at(0), m_RWfullsizeTex.at(1));

	ExecutefadeByGaussCommand(m_RWfullsizeTex.at(0), m_BurstCachedTex);//キャッシュする

	m_burstKernelRegenerate = false;

	PIXEndEvent(m_commandList.Get());
}
void ComputeFilterApp::generateGhost()
{
	PIXBeginEvent(m_commandList.Get(), 0, "GenerateGhost");

	ExecuteDrawFixPolygonCommand(m_RWfullsizeTex.at(1));
	ExecuteClearCommand(m_RWfullsizeTex.at(2));

	//ごみ画像
	ExecuteCopyCommand(m_RfullsizeTex.at(3), m_RWfullsizeTex.at(3));
	ExecuteClearCommand(m_RWfullsizeTex.at(4));
	ExecuteMultiplyCommand(
		m_RWfullsizeTex.at(1), m_RWfullsizeTex.at(2)
		, m_RWfullsizeTex.at(3), m_RWfullsizeTex.at(4)
		, m_RWfullsizeTex.at(5), m_RWfullsizeTex.at(6));
	ExecuteCopyCommand(m_RWfullsizeTex.at(5), m_RWfullsizeTex.at(1));
	ExecuteCopyCommand(m_RWfullsizeTex.at(6), m_RWfullsizeTex.at(2));

	ExecuteFFTCommand(m_RWfullsizeTex.at(1), m_RWfullsizeTex.at(2));

	ExecuteDrawFRFCommand(m_RWfullsizeTex.at(3), m_RWfullsizeTex.at(4));
	ExecuteMultiplyCommand(
		m_RWfullsizeTex.at(1), m_RWfullsizeTex.at(2)
		, m_RWfullsizeTex.at(3), m_RWfullsizeTex.at(4)
		, m_RWfullsizeTex.at(5), m_RWfullsizeTex.at(6));
	ExecuteIFFTCommand(m_RWfullsizeTex.at(5), m_RWfullsizeTex.at(6));

	ExecuteCalcurateIntensityCommand(
		m_RWfullsizeTex.at(5), m_RWfullsizeTex.at(6)
		, m_RWfullsizeTex.at(3), m_RWfullsizeTex.at(4));

	ExecuteToneMappingCommand(m_RWfullsizeTex.at(3), m_RWfullsizeTex.at(4));

	ExecuteApplyCausticCommand(m_RWfullsizeTex.at(4), m_RWfullsizeTex.at(3));

	ExecuteBlurCommand(m_RWfullsizeTex.at(3), m_RWfullsizeTex.at(0));
	ExecuteBlurCommand(m_RWfullsizeTex.at(0), m_RWfullsizeTex.at(1));
	ExecuteBlurCommand(m_RWfullsizeTex.at(1), m_GhostCachedTex);//キャッシュする

	m_ghostKernelRegenerate = false;

	PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::constructLensFlare()
{
	PIXBeginEvent(m_commandList.Get(), 0, "One-Loop");
	{
		{
			ExecuteScalingSizeByRandomTblCommand(m_GhostCachedTex, m_RWfullsizeTex.at(5));
		}
		{
			/*ExecuteInverseRotateImageByRandomTblCommand(m_GhostCachedTex, m_RWfullsizeTex.at(5));
			ExecuteScalingSizeByRandomTblCommand(m_RWfullsizeTex.at(5), m_RWfullsizeTex.at(1));
			ExecuteRotateImageByRandomTblCommand(m_RWfullsizeTex.at(1), m_RWfullsizeTex.at(5));*/
		}
		ExecuteShiftImageByRandomTblansAddCommand(m_RWfullsizeTex.at(5), m_RWfullsizeTex.at(7));

		ExecuteIncrementRandomTblIndexCommand();
	}
	PIXEndEvent(m_commandList.Get());
}

int ComputeFilterApp::ExecuteLensFlareCommand()
{
	//定数バッファの設定
	ComputeParameters cp;
	cp.intervalx = m_interval_x * 1e-6;
	cp.intervaly = cp.intervalx * (m_texwidth * 1.0f) / (m_texheight * 1.0f);
	cp.gausssigma = m_width/2/2/2/2;
	cp.gausssigma2 = m_gausssigma2;
	cp.glareintensity = m_glareintensity;
	cp.threshold = m_threshold;
	cp.distance0 = m_r * m_propdistance[0] / 1000.0f;
	cp.focus0 = m_focus[0] / 1000.0f;
	//固定
	cp.N = m_N;
	cp.r = m_r;
	//動く
	cp.N1 = m_N1;
	cp.r1 = m_r1;
	cp.posx = m_posx;
	cp.posy = m_posy;
	//cp.ghostScale = m_ghostScale;
	cp.rotAngle = m_rotAngle;
	cp.baseColor = m_baseColor;

	WriteToUploadHeapMemory(m_mainComputeCB[0].Get(), sizeof(cp), &cp);
	WriteToUploadHeapMemory(m_ghostPosBuff[0].Get(), sizeof(float) * m_ghostPosTbl.size(), m_ghostPosTbl.data());

	std::vector<float> rotMat;
	float x = m_posx - 0.5, y = m_posy - 0.5;
	float theta = atan2(y , x);
	float Sin = sin(theta);
	float Cos = cos(theta);
	rotMat.push_back(Cos);
	rotMat.push_back(-Sin);
	rotMat.push_back(Sin);
	rotMat.push_back(Cos);
	WriteToUploadHeapMemory(m_ghostRotateMatBuff[0].Get(), sizeof(float) * 4, rotMat.data());
	m_commandList->SetComputeRootConstantBufferView(7, m_ghostPosBuff[0]->GetGPUVirtualAddress());//b1
	m_commandList->SetComputeRootDescriptorTable(8, m_RandomIndexCounter[0].WriteState());//u4
	m_commandList->SetComputeRootConstantBufferView(9, m_ghostRotateMatBuff[0]->GetGPUVirtualAddress());//u4

	if (m_burstKernelRegenerate)
	{
		generateBurst();
	}

	ExecuteScaleingByPosCommand(m_BurstCachedTex, m_RWfullsizeTex.at(2));
	ExecuteShiftImageByTargetPosCommand(m_RWfullsizeTex.at(2), m_RWfullsizeTex.at(1));

	ExecuteCopyCommand(m_RWfullsizeTex.at(1), m_RWfullsizeTex.at(9));
	ExecuteCopyCommand(m_RWfullsizeTex.at(1), m_RWdisplayTex.at(0));

	if (m_ghostKernelRegenerate)
	{
		generateGhost();
	}

	ExecuteCopyCommand(m_GhostCachedTex, m_RWdisplayTex.at(1));

	ExecuteResetRandomTblIndexCommand();

	ExecuteDrawGaussianNoNormalizeCommand2(m_RWfullsizeTex.at(8), m_RWfullsizeTex.at(7));

	PIXBeginEvent(m_commandList.Get(), 0, "Ghost Generate Loop");
	for (int i = 0; i < m_ghostNum; i++)
	{
		constructLensFlare();
	}
	PIXEndEvent(m_commandList.Get());

	ExecuteToneMappingCommand(m_RWfullsizeTex.at(7), m_RWfullsizeTex.at(2));//ないほうが味が出る説 白飛び感が出る
	ExecuteBlurCommand(m_RWfullsizeTex.at(2), m_RWfullsizeTex.at(0));//ないほうが味が出る説 白飛び感が出る
	ExecuteAddCommand(m_RWfullsizeTex.at(0), m_RWfullsizeTex.at(9), m_RWfullsizeTex.at(1));

	int ret = 0;
	switch (m_inputfrag)
	{
	case InputFrag_Image:
		ExecuteCopyCommand(m_RfullsizeTex.at(0), m_RWfullsizeTex.at(0));
		ExecuteAddCommand(m_RWfullsizeTex.at(0), m_RWfullsizeTex.at(1), m_RWfullsizeTex.at(2));
		ret =  2;
		break;
	case InputFrag_Aperture:
		ret =  1;
		break;
	default:
		ret =  1;
	}

	//for (int i = 0; i < m_RWfullsizeTex.size(); i++)
	//{
	//	ExecuteClearCommand(m_RWfullsizeTex.at(i));
	//}

	return ret;
}

