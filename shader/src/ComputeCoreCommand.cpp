#include "../include/ComputeFilterApp.h"

const int ComputeFilterApp::m_threadNum = 16;//‚±‚ê‚ª‘‚¢

void ComputeFilterApp::ExecuteCopyCommand(TextureData& In, AutoTransitionTexture& Out)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(1, In.handleRead);
    m_commandList->SetComputeRootDescriptorTable(3, Out.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "copyCS");
    m_commandList->SetPipelineState(m_pipelines["copyCS"].Get());
    //m_commandList->Dispatch(1, m_texheight, 1);
    m_commandList->Dispatch(m_texwidth / m_threadNum, m_texheight / m_threadNum, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteCopyCommand(AutoTransitionTexture& In, AutoTransitionTexture& Out)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(1, In.ReadState());
    m_commandList->SetComputeRootDescriptorTable(3, Out.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "copyCS");
    m_commandList->SetPipelineState(m_pipelines["copyCS"].Get());
    //m_commandList->Dispatch(1, m_texheight, 1);
    m_commandList->Dispatch(m_texwidth / m_threadNum, m_texheight / m_threadNum, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteClearCommand(AutoTransitionTexture& Tex)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(3, Tex.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "clearCS");
    m_commandList->SetPipelineState(m_pipelines["clearCS"].Get());
    m_commandList->Dispatch(1, m_texheight, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteSwapCommand(AutoTransitionTexture& InReal, AutoTransitionTexture& InImage, AutoTransitionTexture& OutReal, AutoTransitionTexture& OutImage)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(1, InReal.ReadState());
    m_commandList->SetComputeRootDescriptorTable(2, InImage.ReadState());
    m_commandList->SetComputeRootDescriptorTable(3, OutReal.WriteState());
    m_commandList->SetComputeRootDescriptorTable(4, OutImage.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "swapCS");
    m_commandList->SetPipelineState(m_pipelines["swapCS"].Get());
    m_commandList->Dispatch(1, m_texheight, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteAddCommand(AutoTransitionTexture& Tex1, AutoTransitionTexture& Tex2, AutoTransitionTexture& Out)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(1, Tex1.ReadState());
    m_commandList->SetComputeRootDescriptorTable(2, Tex2.ReadState());
    m_commandList->SetComputeRootDescriptorTable(3, Out.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "AddCS");
    m_commandList->SetPipelineState(m_pipelines["AddCS"].Get());
    //m_commandList->Dispatch(1, m_texheight, 1);
    m_commandList->Dispatch(m_texwidth / m_threadNum, m_texheight / m_threadNum, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteColorScaleCommand(AutoTransitionTexture& Tex, AutoTransitionTexture& Out)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(1, Tex.ReadState());
    m_commandList->SetComputeRootDescriptorTable(3, Out.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "colorScaleCS");
    m_commandList->SetPipelineState(m_pipelines["colorScaleCS"].Get());
    m_commandList->Dispatch(1, m_texheight, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteBinaryThresholdCommand(AutoTransitionTexture& In, AutoTransitionTexture& Out)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(1, In.ReadState());
    m_commandList->SetComputeRootDescriptorTable(3, Out.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "BTCS");
    m_commandList->SetPipelineState(m_pipelines["BTCS"].Get());
    m_commandList->Dispatch(1, m_texheight, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteCalcurateAmplitudeCommand(AutoTransitionTexture& InReal, AutoTransitionTexture& InImage, AutoTransitionTexture& OutReal, AutoTransitionTexture& OutImage)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(1, InReal.ReadState());
    m_commandList->SetComputeRootDescriptorTable(2, InImage.ReadState());
    m_commandList->SetComputeRootDescriptorTable(3, OutReal.WriteState());
    m_commandList->SetComputeRootDescriptorTable(4, OutImage.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "ampCS");
    m_commandList->SetPipelineState(m_pipelines["ampCS"].Get());
    m_commandList->Dispatch(1, m_texheight, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteCalcurateIntensityCommand(AutoTransitionTexture& InReal, AutoTransitionTexture& InImage, AutoTransitionTexture& OutReal, AutoTransitionTexture& OutImage)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(1, InReal.ReadState());
    m_commandList->SetComputeRootDescriptorTable(2, InImage.ReadState());
    m_commandList->SetComputeRootDescriptorTable(3, OutReal.WriteState());
    m_commandList->SetComputeRootDescriptorTable(4, OutImage.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "intenCS");
    m_commandList->SetPipelineState(m_pipelines["intenCS"].Get());
    m_commandList->Dispatch(1, m_texheight, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteConvertToFormatCommand(AutoTransitionTexture& InReal, AutoTransitionTexture& InImage, AutoTransitionTexture& OutReal, AutoTransitionTexture& OutImage)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(1, InReal.ReadState());
    m_commandList->SetComputeRootDescriptorTable(2, InImage.ReadState());
    m_commandList->SetComputeRootDescriptorTable(3, OutReal.WriteState());
    m_commandList->SetComputeRootDescriptorTable(4, OutImage.WriteState());
    switch (m_mode)
    {
    case Mode_PROP_REAL:
        PIXBeginEvent(m_commandList.Get(), 0, "realCS");
        m_commandList->SetPipelineState(m_pipelines["realCS"].Get());
        break;
    case Mode_PROP_IMAGE:
        PIXBeginEvent(m_commandList.Get(), 0, "imageCS");
        m_commandList->SetPipelineState(m_pipelines["imageCS"].Get());
        break;
    case Mode_PROP_AMPLITUDE:
        PIXBeginEvent(m_commandList.Get(), 0, "ampCS");
        m_commandList->SetPipelineState(m_pipelines["ampCS"].Get());
        break;
    case Mode_PROP_INTENSITY:
        PIXBeginEvent(m_commandList.Get(), 0, "intenCS");
        m_commandList->SetPipelineState(m_pipelines["intenCS"].Get());
        break;
    case Mode_PROP_PHASE:
        PIXBeginEvent(m_commandList.Get(), 0, "phaseCS");
        m_commandList->SetPipelineState(m_pipelines["phaseCS"].Get());
        break;
    default:
        PIXBeginEvent(m_commandList.Get(), 0, "intenCS");
        m_commandList->SetPipelineState(m_pipelines["intenCS"].Get());
        break;
    }
    m_commandList->Dispatch(1, m_texheight, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteCalcMaxMinCommand(AutoTransitionTexture& Tex, AutoTransitionTexture& OutOnePixReal_MAX, AutoTransitionTexture& OutOnePixImage_MIN)
{
    //Å‘å’lÅ¬’l‚ÌŒvŽZ
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(1, Tex.ReadState());
    m_commandList->SetComputeRootDescriptorTable(3, m_RWlineInnerTex.at(0).WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "maxminfirstCS");
    m_commandList->SetPipelineState(m_pipelines["maxminfirstCS"].Get());
    m_commandList->Dispatch(m_texwidth, 1, 1);
    PIXEndEvent(m_commandList.Get());

    m_commandList->SetComputeRootDescriptorTable(1, m_RWlineInnerTex.at(0).ReadState());
    m_commandList->SetComputeRootDescriptorTable(3, OutOnePixReal_MAX.WriteState());
    m_commandList->SetComputeRootDescriptorTable(4, OutOnePixImage_MIN.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "maxminsecondCS");
    m_commandList->SetPipelineState(m_pipelines["maxminsecondCS"].Get());
    m_commandList->Dispatch(1, 1, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteMaxMinNormalizeCommand(AutoTransitionTexture& OutOnePixReal_MAX, AutoTransitionTexture& OutOnePixImage_MIN, AutoTransitionTexture& Tex)
{
    ExecuteCopyCommand(Tex, m_RWfullsizeInnerTex.at(6));
    //³‹K‰»
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(1, OutOnePixReal_MAX.ReadState());
    m_commandList->SetComputeRootDescriptorTable(2, OutOnePixImage_MIN.ReadState());
    m_commandList->SetComputeRootDescriptorTable(3, Tex.WriteState());
    m_commandList->SetComputeRootDescriptorTable(4, m_RWfullsizeInnerTex.at(6).ReadState());
    PIXBeginEvent(m_commandList.Get(), 0, "normalizeCS");
    m_commandList->SetPipelineState(m_pipelines["normalizeCS"].Get());
    m_commandList->Dispatch(1, m_texheight, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteNormalizeWaveFrontCommand(AutoTransitionTexture& InReal, AutoTransitionTexture& InImage, AutoTransitionTexture& OutReal, AutoTransitionTexture& OutImage)
{
    //Šeo—Íƒpƒ^[ƒ“‚É‘Î‰ž‚µ‚½³‹K‰»
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(1, InReal.ReadState());
    m_commandList->SetComputeRootDescriptorTable(2, InImage.ReadState());
    m_commandList->SetComputeRootDescriptorTable(3, OutReal.WriteState());
    m_commandList->SetComputeRootDescriptorTable(4, OutImage.WriteState());
    switch (m_mode)
    {
    case Mode_PROP_REAL:
        PIXBeginEvent(m_commandList.Get(), 0, "REALIMAGE4DISPCS");
        m_commandList->SetPipelineState(m_pipelines["REALIMAGE4DISPCS"].Get());
        break;
    case Mode_PROP_IMAGE:
        PIXBeginEvent(m_commandList.Get(), 0, "REALIMAGE4DISPCS");
        m_commandList->SetPipelineState(m_pipelines["REALIMAGE4DISPCS"].Get());
        break;
    case Mode_PROP_AMPLITUDE:
        PIXBeginEvent(m_commandList.Get(), 0, "AMPINT4DISPCS");
        m_commandList->SetPipelineState(m_pipelines["AMPINT4DISPCS"].Get());
        break;
    case Mode_PROP_INTENSITY:
        PIXBeginEvent(m_commandList.Get(), 0, "AMPINT4DISPCS");
        m_commandList->SetPipelineState(m_pipelines["AMPINT4DISPCS"].Get());
        break;
    case Mode_PROP_PHASE:
        PIXBeginEvent(m_commandList.Get(), 0, "PHASE4DISPCS");
        m_commandList->SetPipelineState(m_pipelines["PHASE4DISPCS"].Get());
        break;
    default:
        PIXBeginEvent(m_commandList.Get(), 0, "AMPINT4DISPCS");
        m_commandList->SetPipelineState(m_pipelines["AMPINT4DISPCS"].Get());
        break;
    }
    m_commandList->Dispatch(1, m_texheight, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteOrdinalNormalizeCommand(AutoTransitionTexture& Tex)
{
    ExecuteCalcMaxMinCommand(Tex, m_RWmaxminInnerTex.at(0), m_RWmaxminInnerTex.at(1));

    ExecuteMaxMinNormalizeCommand(m_RWmaxminInnerTex.at(0), m_RWmaxminInnerTex.at(1), Tex);
}

void ComputeFilterApp::ExecuteFFTCommand(AutoTransitionTexture& Real, AutoTransitionTexture& Image)
{
    ExecuteSwapCommand(Real, Image, m_RWfullsizeInnerTex.at(0), m_RWfullsizeInnerTex.at(1));
    //c•ûŒü
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(1, m_RWfullsizeInnerTex.at(0).ReadState());
    m_commandList->SetComputeRootDescriptorTable(2, m_RWfullsizeInnerTex.at(1).ReadState());
    m_commandList->SetComputeRootDescriptorTable(3, m_RWfullsizeInnerTex.at(2).WriteState());
    m_commandList->SetComputeRootDescriptorTable(4, m_RWfullsizeInnerTex.at(3).WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "fftCS_ROW");
    m_commandList->SetPipelineState(m_pipelines["fftCS_ROW"].Get());
    m_commandList->Dispatch(1, m_texheight, 1);
    PIXEndEvent(m_commandList.Get());

    //‰¡•ûŒü
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(1, m_RWfullsizeInnerTex.at(2).ReadState());
    m_commandList->SetComputeRootDescriptorTable(2, m_RWfullsizeInnerTex.at(3).ReadState());
    m_commandList->SetComputeRootDescriptorTable(3, m_RWfullsizeInnerTex.at(0).WriteState());
    m_commandList->SetComputeRootDescriptorTable(4, m_RWfullsizeInnerTex.at(1).WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "fftCS_COL");
    m_commandList->SetPipelineState(m_pipelines["fftCS_COL"].Get());
    m_commandList->Dispatch(1, m_texwidth, 1);
    PIXEndEvent(m_commandList.Get());

    ExecuteSwapCommand(m_RWfullsizeInnerTex.at(0), m_RWfullsizeInnerTex.at(1), Real, Image);
}

void ComputeFilterApp::ExecuteIFFTCommand(AutoTransitionTexture& Real, AutoTransitionTexture& Image)
{
    ExecuteSwapCommand(Real, Image, m_RWfullsizeInnerTex.at(0), m_RWfullsizeInnerTex.at(1));
    //c•ûŒü
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(1, m_RWfullsizeInnerTex.at(0).ReadState());
    m_commandList->SetComputeRootDescriptorTable(2, m_RWfullsizeInnerTex.at(1).ReadState());
    m_commandList->SetComputeRootDescriptorTable(3, m_RWfullsizeInnerTex.at(2).WriteState());
    m_commandList->SetComputeRootDescriptorTable(4, m_RWfullsizeInnerTex.at(3).WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "ifftCS_ROW");
    m_commandList->SetPipelineState(m_pipelines["ifftCS_ROW"].Get());
    m_commandList->Dispatch(1, m_texheight, 1);
    PIXEndEvent(m_commandList.Get());
    //‰¡•ûŒü
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(1, m_RWfullsizeInnerTex.at(2).ReadState());
    m_commandList->SetComputeRootDescriptorTable(2, m_RWfullsizeInnerTex.at(3).ReadState());
    m_commandList->SetComputeRootDescriptorTable(3, m_RWfullsizeInnerTex.at(0).WriteState());
    m_commandList->SetComputeRootDescriptorTable(4, m_RWfullsizeInnerTex.at(1).WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "ifftCS_COL");
    m_commandList->SetPipelineState(m_pipelines["ifftCS_COL"].Get());
    m_commandList->Dispatch(1, m_texwidth, 1);
    PIXEndEvent(m_commandList.Get());
    ExecuteSwapCommand(m_RWfullsizeInnerTex.at(0), m_RWfullsizeInnerTex.at(1), Real, Image);
}

void ComputeFilterApp::ExecuteDrawFRFCommand(AutoTransitionTexture& Real, AutoTransitionTexture& Image, int num)
{
    //’†SŠî€‚É•`‰æ
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(3, Real.WriteState());
    m_commandList->SetComputeRootDescriptorTable(4, Image.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "frfCS0");
    m_commandList->SetPipelineState(m_pipelines["frfCS0"].Get());
    m_commandList->Dispatch(1, m_texheight, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteDrawGaussianCommand(AutoTransitionTexture& Real, AutoTransitionTexture& Image)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(3, Real.WriteState());
    m_commandList->SetComputeRootDescriptorTable(4, Image.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "gaussianCS");
    m_commandList->SetPipelineState(m_pipelines["gaussianCS"].Get());
    m_commandList->Dispatch(1, m_texheight, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteDrawGaussianNoNormalizeCommand(AutoTransitionTexture& Real, AutoTransitionTexture& Image)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(3, Real.WriteState());
    m_commandList->SetComputeRootDescriptorTable(4, Image.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "gaussianNoNormalizeCS");
    m_commandList->SetPipelineState(m_pipelines["gaussianNoNormalizeCS"].Get());
    m_commandList->Dispatch(1, m_texheight, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteDrawGaussianNoNormalizeCommand2(AutoTransitionTexture& Real, AutoTransitionTexture& Image)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(3, Real.WriteState());
    m_commandList->SetComputeRootDescriptorTable(4, Image.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "gaussianNoNormalize2CS");
    m_commandList->SetPipelineState(m_pipelines["gaussianNoNormalize2CS"].Get());
    m_commandList->Dispatch(1, m_texheight, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteConvolutionCommand(AutoTransitionTexture& InReal0, AutoTransitionTexture& InImage0, AutoTransitionTexture& InReal1, AutoTransitionTexture& InImage1, AutoTransitionTexture& OutReal, AutoTransitionTexture& OutImage)
{
    PIXBeginEvent(m_commandList.Get(), 0, "CONVOLUTION");
    //ˆê•û‚ÌFFT
    ExecuteFFTCommand(InReal0, InImage0);
    //‚à‚¤ˆê•û‚ÌFFT
    ExecuteFFTCommand(InReal1, InImage1);
    //æŽZ
    ExecuteMultiplyCommand(InReal0, InImage0, InReal1, InImage1, OutReal, OutImage);
    //‹tFFT
    ExecuteIFFTCommand(OutReal, OutImage);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteMultiplyCommand(AutoTransitionTexture& InReal0, AutoTransitionTexture& InImage0, AutoTransitionTexture& InReal1, AutoTransitionTexture& InImage1, AutoTransitionTexture& OutReal, AutoTransitionTexture& OutImage)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(1, InReal0.ReadState());
    m_commandList->SetComputeRootDescriptorTable(2, InImage0.ReadState());
    m_commandList->SetComputeRootDescriptorTable(3, InReal1.ReadState());
    m_commandList->SetComputeRootDescriptorTable(4, InImage1.ReadState());
    m_commandList->SetComputeRootDescriptorTable(5, OutReal.WriteState());
    m_commandList->SetComputeRootDescriptorTable(6, OutImage.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "mulCS");
    m_commandList->SetPipelineState(m_pipelines["mulCS"].Get());
    m_commandList->Dispatch(1, m_texheight, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteDivideMaxAmpCommand(AutoTransitionTexture& OutOnePixReal_MAX, AutoTransitionTexture& OutOnePixImage_MIN, AutoTransitionTexture& InReal, AutoTransitionTexture& InImage, AutoTransitionTexture& OutReal, AutoTransitionTexture& OutImage)
{
    //Å‘åU•‚É‚æ‚éœŽZ
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(1, OutOnePixReal_MAX.ReadState());
    m_commandList->SetComputeRootDescriptorTable(2, OutOnePixImage_MIN.ReadState());
    m_commandList->SetComputeRootDescriptorTable(3, InReal.ReadState());
    m_commandList->SetComputeRootDescriptorTable(4, InImage.ReadState());
    m_commandList->SetComputeRootDescriptorTable(5, OutReal.WriteState());
    m_commandList->SetComputeRootDescriptorTable(6, OutImage.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "divByMaxAMPCS");
    m_commandList->SetPipelineState(m_pipelines["divByMaxAMPCS"].Get());
    m_commandList->Dispatch(1, m_texheight, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteSpectrumScalingCommand(AutoTransitionTexture& InReal, AutoTransitionTexture& InImage, AutoTransitionTexture& OutReal, AutoTransitionTexture& OutImage)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(1, InReal.ReadState());
    m_commandList->SetComputeRootDescriptorTable(2, InImage.ReadState());
    m_commandList->SetComputeRootDescriptorTable(3, OutReal.WriteState());
    m_commandList->SetComputeRootDescriptorTable(4, OutImage.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "spectrumScalingCS");
    m_commandList->SetPipelineState(m_pipelines["spectrumScalingCS"].Get());
    m_commandList->Dispatch(1, m_texheight, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteSpectrumScalingCommand2(AutoTransitionTexture& InReal, AutoTransitionTexture& InImage, AutoTransitionTexture& Out)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(1, InReal.ReadState());
    m_commandList->SetComputeRootDescriptorTable(2, InImage.ReadState());
    m_commandList->SetComputeRootDescriptorTable(3, Out.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "spectrumScalingCS");
    m_commandList->SetPipelineState(m_pipelines["spectrumScalingCS"].Get());
    m_commandList->Dispatch(1, m_texheight, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteRaiseRICommand(AutoTransitionTexture& InReal, AutoTransitionTexture& InImage, AutoTransitionTexture& OutReal, AutoTransitionTexture& OutImage)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(1, InReal.ReadState());
    m_commandList->SetComputeRootDescriptorTable(2, InImage.ReadState());
    m_commandList->SetComputeRootDescriptorTable(3, OutReal.WriteState());
    m_commandList->SetComputeRootDescriptorTable(4, OutImage.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "raiseRICS");
    m_commandList->SetPipelineState(m_pipelines["raiseRICS"].Get());
    m_commandList->Dispatch(1, m_texheight, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteRadialBlurCommand(AutoTransitionTexture& In, AutoTransitionTexture& Out)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(1, In.ReadState());
    m_commandList->SetComputeRootDescriptorTable(3, Out.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "radialblurCS");
    m_commandList->SetPipelineState(m_pipelines["radialblurCS"].Get());
    m_commandList->Dispatch(1, m_texheight, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteDrawDecayLineCommand(AutoTransitionTexture& Tex)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(3, Tex.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "decaylineCS");
    m_commandList->SetPipelineState(m_pipelines["decaylineCS"].Get());
    m_commandList->Dispatch(1, m_texheight, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteDrawQuadraticPhaseCommand(AutoTransitionTexture& Real, AutoTransitionTexture& Image, int num)
{
    //’†SŠî€‚É•`‰æ
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(3, Real.WriteState());
    m_commandList->SetComputeRootDescriptorTable(4, Image.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "quadraticCS0");
    m_commandList->SetPipelineState(m_pipelines["quadraticCS0"].Get());

    m_commandList->Dispatch(1, m_texheight, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteDrawPolygonCommand(AutoTransitionTexture& Tex)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(3, Tex.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "drawpolygonCS");
    m_commandList->SetPipelineState(m_pipelines["drawpolygonCS"].Get());
    m_commandList->Dispatch(1, m_texheight, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteDrawFixPolygonCommand(AutoTransitionTexture& Tex)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(3, Tex.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "drawfixpolygonCS");
    m_commandList->SetPipelineState(m_pipelines["drawfixpolygonCS"].Get());
    m_commandList->Dispatch(1, m_texheight, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteDrawMovingPolygonCommand(AutoTransitionTexture& Tex)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(3, Tex.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "drawmovingpolygonCS");
    m_commandList->SetPipelineState(m_pipelines["drawmovingpolygonCS"].Get());
    m_commandList->Dispatch(1, m_texheight, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteAmplitudeAdjustmentCommand(AutoTransitionTexture& In, AutoTransitionTexture& Out)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(1, In.ReadState());
    m_commandList->SetComputeRootDescriptorTable(3, Out.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "amplitudeadjustmentCS");
    m_commandList->SetPipelineState(m_pipelines["amplitudeadjustmentCS"].Get());
    m_commandList->Dispatch(1, m_texheight, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteToneMappingCommand(AutoTransitionTexture& In, AutoTransitionTexture& Out)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(1, In.ReadState());
    m_commandList->SetComputeRootDescriptorTable(3, Out.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "tonemappingCS");
    m_commandList->SetPipelineState(m_pipelines["tonemappingCS"].Get());
    m_commandList->Dispatch(1, m_texheight, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteDrawMovingDotCommand(AutoTransitionTexture& Tex)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(3, Tex.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "drawMovingDotCS");
    m_commandList->SetPipelineState(m_pipelines["drawMovingDotCS"].Get());
    m_commandList->Dispatch(1, m_texheight, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteDrawSmallCircCommand(AutoTransitionTexture& Tex)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(3, Tex.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "drawSmallCircCS");
    m_commandList->SetPipelineState(m_pipelines["drawSmallCircCS"].Get());
    m_commandList->Dispatch(1, m_texheight, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteDrawDotByRandomTblCommand(AutoTransitionTexture& Tex)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(3, Tex.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "drawDotByRandomTblCS");
    m_commandList->SetPipelineState(m_pipelines["drawDotByRandomTblCS"].Get());
    m_commandList->Dispatch(1, m_texheight, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteIncrementRandomTblIndexCommand()
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    PIXBeginEvent(m_commandList.Get(), 0, "incrementRandomIndexCS");
    m_commandList->SetPipelineState(m_pipelines["incrementRandomIndexCS"].Get());
    m_commandList->Dispatch(1, 1, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteResetRandomTblIndexCommand()
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    PIXBeginEvent(m_commandList.Get(), 0, "resetRandomIndexCS");
    m_commandList->SetPipelineState(m_pipelines["resetRandomIndexCS"].Get());
    m_commandList->Dispatch(1, 1, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteScalingSizeByRandomTblCommand(AutoTransitionTexture& In, AutoTransitionTexture& Out)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(1, In.ReadState());
    m_commandList->SetComputeRootDescriptorTable(3, Out.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "scalingSizeByRandomTblCS");
    m_commandList->SetPipelineState(m_pipelines["scalingSizeByRandomTblCS"].Get());
    //m_commandList->Dispatch(1, m_texheight, 1);
    m_commandList->Dispatch(m_texwidth / m_threadNum, m_texheight / m_threadNum, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteShiftImageByRandomTblCommand(AutoTransitionTexture& In, AutoTransitionTexture& Out)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(1, In.ReadState());
    m_commandList->SetComputeRootDescriptorTable(3, Out.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "shiftImageByRandomTblCS");
    m_commandList->SetPipelineState(m_pipelines["shiftImageByRandomTblCS"].Get());
    //m_commandList->Dispatch(1, m_texheight, 1);
    m_commandList->Dispatch(m_texwidth / m_threadNum, m_texheight / m_threadNum, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteShiftImageByRandomTblansAddCommand(AutoTransitionTexture& In, AutoTransitionTexture& Out)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(1, In.ReadState());
    m_commandList->SetComputeRootDescriptorTable(3, Out.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "shiftImageByRandomTblansAddCS");
    m_commandList->SetPipelineState(m_pipelines["shiftImageByRandomTblansAddCS"].Get());
    //m_commandList->Dispatch(1, m_texheight, 1);
    m_commandList->Dispatch(m_texwidth / m_threadNum, m_texheight / m_threadNum, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteShiftImageByTargetPosCommand(AutoTransitionTexture& In, AutoTransitionTexture& Out)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(1, In.ReadState());
    m_commandList->SetComputeRootDescriptorTable(3, Out.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "shiftImageByTargetPosCS");
    m_commandList->SetPipelineState(m_pipelines["shiftImageByTargetPosCS"].Get());
    m_commandList->Dispatch(1, m_texheight, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteRotateImageByRandomTblCommand(AutoTransitionTexture& In, AutoTransitionTexture& Out)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(1, In.ReadState());
    m_commandList->SetComputeRootDescriptorTable(3, Out.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "rotateImageByRandomTblCS");
    m_commandList->SetPipelineState(m_pipelines["rotateImageByRandomTblCS"].Get());
    //m_commandList->Dispatch(1, m_texheight, 1);
   m_commandList->Dispatch(m_texwidth / m_threadNum, m_texheight / m_threadNum, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteInverseRotateImageByRandomTblCommand(AutoTransitionTexture& In, AutoTransitionTexture& Out)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(1, In.ReadState());
    m_commandList->SetComputeRootDescriptorTable(3, Out.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "invRotateImageByRandomTblCS");
    m_commandList->SetPipelineState(m_pipelines["invRotateImageByRandomTblCS"].Get());
    //m_commandList->Dispatch(1, m_texheight, 1);
    m_commandList->Dispatch(m_texwidth / m_threadNum, m_texheight / m_threadNum, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteApplyCausticCommand(AutoTransitionTexture& In, AutoTransitionTexture& Out)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(1, In.ReadState());
    m_commandList->SetComputeRootDescriptorTable(3, Out.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "causticCS");
    m_commandList->SetPipelineState(m_pipelines["causticCS"].Get());
    //m_commandList->Dispatch(1, m_texheight, 1);
    m_commandList->Dispatch(m_texwidth / m_threadNum, m_texheight / m_threadNum, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteApplyVignettingByRandomTblCommand(AutoTransitionTexture& KerareMask, AutoTransitionTexture& In, AutoTransitionTexture& Out)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(1, KerareMask.ReadState());
    m_commandList->SetComputeRootDescriptorTable(2, In.ReadState());
    m_commandList->SetComputeRootDescriptorTable(3, Out.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "applyVignettingByRandomTblCS");
    m_commandList->SetPipelineState(m_pipelines["applyVignettingByRandomTblCS"].Get());
    //m_commandList->Dispatch(1, m_texheight, 1);
    m_commandList->Dispatch(m_texwidth / m_threadNum, m_texheight / m_threadNum, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteBlurCommand(AutoTransitionTexture& In, AutoTransitionTexture& Out)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(1, In.ReadState());
    m_commandList->SetComputeRootDescriptorTable(3, Out.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "blurCS");
    m_commandList->SetPipelineState(m_pipelines["blurCS"].Get());
    //m_commandList->Dispatch(1, m_texheight, 1);
    m_commandList->Dispatch(m_texwidth / m_threadNum, m_texheight / m_threadNum, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecutefadeByGaussCommand(AutoTransitionTexture& In, AutoTransitionTexture& Out)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(1, In.ReadState());
    m_commandList->SetComputeRootDescriptorTable(3, Out.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "fadeByGaussCS");
    m_commandList->SetPipelineState(m_pipelines["fadeByGaussCS"].Get());
    //m_commandList->Dispatch(1, m_texheight, 1);
    m_commandList->Dispatch(m_texwidth / m_threadNum, m_texheight / m_threadNum, 1);
    PIXEndEvent(m_commandList.Get());
}

void ComputeFilterApp::ExecuteScaleingByPosCommand(AutoTransitionTexture& In, AutoTransitionTexture& Out)
{
    m_commandList->SetComputeRootConstantBufferView(0, m_mainComputeCB[0]->GetGPUVirtualAddress());
    m_commandList->SetComputeRootDescriptorTable(1, In.ReadState());
    m_commandList->SetComputeRootDescriptorTable(3, Out.WriteState());
    PIXBeginEvent(m_commandList.Get(), 0, "scaleByPosCS");
    m_commandList->SetPipelineState(m_pipelines["scaleByPosCS"].Get());
    //m_commandList->Dispatch(1, m_texheight, 1);
    m_commandList->Dispatch(m_texwidth / m_threadNum, m_texheight / m_threadNum, 1);
    PIXEndEvent(m_commandList.Get());
}