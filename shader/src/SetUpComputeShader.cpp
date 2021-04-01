#define _CRT_SECURE_NO_WARNINGS
#define _SILENCE_ALL_CXX17_DEPRECATION_WARNINGS

#include "../include/ComputeFilterApp.h"

#include <stdlib.h>
#include <DirectXTex.h>
#include <fstream>
#include <codecvt>

#include <utility>

#include <random>

using namespace std;
using namespace DirectX;

void ComputeFilterApp::randomGenerate()
{
    //乱数再生成
    random_device device;

    //0.5で光源位置 -0.5で対象な位置

    //normal_distribution<> normalDist(-0.5, 0.05);

    //normal_distribution<> normalDist1(0.5, 0.1);

    //normal_distribution<> normalDist2(-0.2, 0.1);

    normal_distribution<> normalDist(-1.0, m_randSigmaOppositeFar);

    normal_distribution<> normalDist1(0.0, m_randSigma);

    normal_distribution<> normalDist2(-0.5, m_randSigmaOppositeNear);

    mt19937 mt(device());

    m_ghostPosTbl.clear();
    m_ghostPosTbl.resize(ComputeFilterApp::m_ghostNum);

    int count = 0;
    for (auto& v : m_ghostPosTbl)
    {
        if (count < m_ghostPosTbl.size() / 2.0f)
        {
            v = normalDist2(mt);
            count++;
            continue;
        }
            
        if (count >= 3.0 * m_ghostPosTbl.size() / 4.0f)
        {
            v = normalDist1(mt);
            count++;
            continue;
        }

        else
        {
            v = normalDist(mt);
            count++;
            continue;
        }
    }

    m_randomgenerate = false;
}

void ComputeFilterApp::Prepare()//親のInitialize内で呼ばれる
{
    SetTitle("ComputeFilter");
    CreateRootSignatures();

    m_commandList->Reset(m_commandAllocators[0].Get(), nullptr);
    ID3D12DescriptorHeap* heaps[] = { m_heap->GetHeap().Get() };
    m_commandList->SetDescriptorHeaps(1, heaps);
    m_commandList->Close();

    //PrepareSimpleModel();
    PrepareComputeFilter();
    PrepareSimpleModel();
    PreparePipeline();
    PrepareConstantBuffers();
}

void ComputeFilterApp::PrepareConstantBuffers()
{
    auto cbDesc = CD3DX12_RESOURCE_DESC::Buffer(
        sizeof(SceneParameters)
    );
    m_mainSceneCB = CreateConstantBuffers(cbDesc);

    auto cbDescC = CD3DX12_RESOURCE_DESC::Buffer(
        sizeof(ComputeParameters)
    );

    m_mainComputeCB = CreateConstantBuffers(cbDescC);

    ComputeParameters cp;
    cp.intervalx = 1e-6;
    cp.intervaly = cp.intervalx * (m_texwidth * 1.0f) / (m_texheight * 1.0f);
    cp.distance0 = m_propdistance[0];
    cp.gausssigma = m_gausssigma;
    cp.glareintensity = m_glareintensity;
    cp.threshold = m_threshold;

    WriteToUploadHeapMemory(m_mainComputeCB[0].Get(), sizeof(cp), &cp);

    auto size = sizeof(float) * m_ghostNum;
    auto alignedsize = (size + (D3D12_CONSTANT_BUFFER_DATA_PLACEMENT_ALIGNMENT - 1))&~(D3D12_CONSTANT_BUFFER_DATA_PLACEMENT_ALIGNMENT - 1);

    auto cbDescR = CD3DX12_RESOURCE_DESC::Buffer(
        alignedsize
    );

    m_ghostPosBuff = CreateConstantBuffers(cbDescR, 1);

    auto cbDescRM = CD3DX12_RESOURCE_DESC::Buffer(
        sizeof(float) * 4
    );

    m_ghostRotateMatBuff = CreateConstantBuffers(cbDescRM, 1);
}

void ComputeFilterApp::CreateRootSignatures()//ルートシグネチャの作成 ディスクリプタテーブルをまとめるもの 　ここではPS/VSの定数などを定義
{
    //もともと
  // RootSignature
    D3D12_DESCRIPTOR_RANGE descSrv{};
    descSrv.RangeType = D3D12_DESCRIPTOR_RANGE_TYPE_SRV;
    descSrv.BaseShaderRegister = 0;
    descSrv.NumDescriptors = 1;

    array<CD3DX12_ROOT_PARAMETER, 2> rootParams;//ルートパラメタはディスクリプタテーブルの実体 
    ////////////////////////////////////////////////////////////////////////////////////////////////ディスクリプタヒープ(定数やテクスチャのアドレスのまとめ役)とシェーダのレジスタを関連付ける
    rootParams[0].InitAsConstantBufferView(0);//ConstantBuffer<SceneParameters> sceneConstants : register(b0);
    rootParams[1].InitAsDescriptorTable(1, &descSrv);//Texture2D imageTex : register(t0);

    array<CD3DX12_STATIC_SAMPLER_DESC, 1> samplerDesc;
    samplerDesc[0].Init(0);//SamplerState imageSampler : register(s0);

    CD3DX12_ROOT_SIGNATURE_DESC rootSignatureDesc{};
    rootSignatureDesc.Init(
        UINT(rootParams.size()), rootParams.data(),
        UINT(samplerDesc.size()), samplerDesc.data(),
        D3D12_ROOT_SIGNATURE_FLAG_ALLOW_INPUT_ASSEMBLER_INPUT_LAYOUT);

    ComPtr<ID3DBlob> signature, errBlob;
    D3D12SerializeRootSignature(&rootSignatureDesc,
        D3D_ROOT_SIGNATURE_VERSION_1_0, &signature, &errBlob);
    m_device->CreateRootSignature(0, signature->GetBufferPointer(), signature->GetBufferSize(), IID_PPV_ARGS(&m_rootSignature));
}

void ComputeFilterApp::PreparePipeline()//パイプラインステートの準備 パイプラインステートは頂点レイアウトや頂点シェーダのような表示までの流れに必要なパラメタをまとめるもの
{
    D3D12_INPUT_ELEMENT_DESC inputElementDesc[] = {
      { "POSITION",   0, DXGI_FORMAT_R32G32B32_FLOAT, 0, D3D12_APPEND_ALIGNED_ELEMENT, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0 },
      { "TEXCOORD",   0, DXGI_FORMAT_R32G32_FLOAT, 0, D3D12_APPEND_ALIGNED_ELEMENT, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0 },
    };

    std::vector<wstring> flags;
    std::vector<Shader::DefineMacro> defines;

    Shader shaderVS, shaderPS;
    shaderVS.load(L"ComputeFilter.hlsl", Shader::Vertex, L"mainVS", flags, defines);//mainVSのエントリを探す(エントリポイントは最初に実行される関数)
    shaderPS.load(L"ComputeFilter.hlsl", Shader::Pixel, L"mainPS", flags, defines);//mainPSのエントリを探す


    auto rasterizerState = CD3DX12_RASTERIZER_DESC(D3D12_DEFAULT);
    rasterizerState.CullMode = D3D12_CULL_MODE_NONE;

    auto psoDesc = book_util::CreateDefaultPsoDesc(
        DXGI_FORMAT_R8G8B8A8_UNORM,
        rasterizerState,
        inputElementDesc, _countof(inputElementDesc),
        m_rootSignature,
        shaderVS.getCode(), shaderPS.getCode()
    );

    HRESULT hr;
    PipelineState pipeline;
    hr = m_device->CreateGraphicsPipelineState(&psoDesc, IID_PPV_ARGS(&pipeline));
    ThrowIfFailed(hr, "CreateGraphicsPipelineState Failed.");
    m_pipelines["default"] = pipeline;
}

void ComputeFilterApp::createHeaderHLSL(UINT width, UINT height, bool raw, bool inv)
{
    ofstream file("param.hlsl");

    if (!file) {
        cerr << "Error: file not opened." << endl;
        exit(0);
    }

    file << "#define" << " " << "WIDTH" << " " << width << "\n";
    file << "#define" << " " << "HEIGHT" << " " << height << "\n";
    file << "#define" << " " << "PI" << " " << 3.14159265 << "\n";
    file << "#define" << " " << "RAD" << " " << "PI / 180.0" << "\n";
    file << "#define" << " " << "RANDNUM" << " " << m_ghostNum << "\n";
    file << "#define" << " " << "THREADNUM" << " " << m_threadNum << "\n";

    if (raw)
    {
        file << "#define" << " " << "ROW" << " " << "1" << "\n";
        file << "#define" << " " << "LENGTH" << " " << "WIDTH" << "\n";
        file << "#define" << " " << "BUTTERFLY_COUNT" << " " << inv2pow(width) << "\n";
    }
    else
    {
        file << "#define" << " " << "ROW" << " " << "0" << "\n";
        file << "#define" << " " << "LENGTH" << " " << "HEIGHT" << "\n";
        file << "#define" << " " << "BUTTERFLY_COUNT" << " " << inv2pow(height) << "\n";
    }

    if (inv)
    {
        file << "#define" << " " << "INVERSE" << " " << "1" << "\n";
    }
    else
    {
        file << "#define" << " " << "INVERSE" << " " << "0" << "\n";
    }

    file.close();
}

void ComputeFilterApp::PrepareComputeFilter()//コンピュートシェーダを使用するための準備 関数のコンパイル 画像はここで決める
{
    ImageSize size1, size2, size3, size4;
    vector<ImageSize> size;
   
    vector<wstring> names;

    std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> converter;

    m_RfullsizeTex.resize(0);

    const string imageslash = "image/";

    //const char* から wstringに変換
    for (auto& name : m_textureNames)
    {
        string s(imageslash + name);
        std::wstring wide = converter.from_bytes(s);
        ImageSize ss;
        TextureData tex = LoadTextureFromFile(wide, ss);
        m_RfullsizeTex.push_back(tex);
        size.push_back(ss);
    }

    if (!((size[0] == size[1]) && (size[1] == size[2]) && (size[2] == size[3]) &&(size[3] == size[0])))
    {
        throw std::runtime_error("image size missmatch");
    }

    int RWtexNum = 10;//読み書きテクスチャ枚数 FFTは一つに4枚か
    int RWmaxmintexNum = 2;//最大値最小値格納用テクスチャ枚数
    HRESULT hr;

    m_RWfullsizeTex.resize(0);
    m_RWfullsizeInnerTex.resize(0);
    m_RWmaxminTex.resize(0);
    m_RWmaxminInnerTex.resize(0);
    m_RWlineInnerTex.resize(0);
    m_RWdisplayTex.resize(0);

    m_RWfullsizeTex.shrink_to_fit();
    m_RWfullsizeInnerTex.shrink_to_fit();
    m_RWmaxminTex.shrink_to_fit();
    m_RWmaxminInnerTex.shrink_to_fit();
    m_RWlineInnerTex.shrink_to_fit();
    m_RWdisplayTex.shrink_to_fit();

    for (int i = 0; i < RWtexNum; i++)//フルサイズ
    {
        AutoTransitionTexture tex(m_heap, m_device, m_commandList, m_texwidth, m_texheight, DXGI_FORMAT_R16G16B16A16_FLOAT);
        wchar_t name[30];
        swprintf(name, 30, L"m_RWfullsizeTex[%d]", i);
        tex.GetTextureData().texture->SetName(name);
        m_RWfullsizeTex.push_back(tex);
    }

    for (int i = 0; i < RWtexNum; i++)//フルサイズ(内部計算用)
    {
        AutoTransitionTexture tex(m_heap, m_device, m_commandList, m_texwidth, m_texheight, DXGI_FORMAT_R16G16B16A16_FLOAT);
        wchar_t name[30];
        swprintf(name, 30, L"m_RWfullsizeInnerTex[%d]", i);
        tex.GetTextureData().texture->SetName(name);
        m_RWfullsizeInnerTex.push_back(tex);
    }

    for (int i = 0; i < 4; i++)//フルサイズ(表示計算用)
    {
        AutoTransitionTexture tex(m_heap, m_device, m_commandList, m_texwidth, m_texheight, DXGI_FORMAT_R16G16B16A16_FLOAT);
        wchar_t name[30];
        swprintf(name, 30, L"m_RWdisplayTex[%d]", i);
        tex.GetTextureData().texture->SetName(name);
        m_RWdisplayTex.push_back(tex);
    }

    for (int i = 0; i < RWmaxmintexNum; i++)//最大最小用
    {
        AutoTransitionTexture tex(m_heap, m_device, m_commandList, 1, 1, DXGI_FORMAT_R16G16B16A16_FLOAT);
        wchar_t name[30];
        swprintf(name, 30, L"m_RWmaxminTex[%d]", i);
        tex.GetTextureData().texture->SetName(name);
        m_RWmaxminTex.push_back(tex);
    }

    for (int i = 0; i < RWmaxmintexNum; i++)//最大最小用(内部計算用)
    {
        AutoTransitionTexture tex(m_heap, m_device, m_commandList, 1, 1, DXGI_FORMAT_R16G16B16A16_FLOAT);
        wchar_t name[30];
        swprintf(name, 30, L"m_RWmaxminInnerTex[%d]", i);
        tex.GetTextureData().texture->SetName(name);
        m_RWmaxminInnerTex.push_back(tex);
    }

    {
        AutoTransitionTexture tex(m_heap, m_device, m_commandList, m_texheight, 2, DXGI_FORMAT_R16G16B16A16_FLOAT);
        wchar_t name[30];
        swprintf(name, 30, L"m_RWlineInnerTex[%d]", 0);
        tex.GetTextureData().texture->SetName(name);
        m_RWlineInnerTex.push_back(tex);
    }

    {
        AutoTransitionTexture tex(m_heap, m_device, m_commandList, m_texwidth, m_texheight, DXGI_FORMAT_R16G16B16A16_FLOAT);
        wchar_t name[30] = L"m_RWGhostCachedTex";
        tex.GetTextureData().texture->SetName(name);
        m_GhostCachedTex = tex;
    }

    {
        AutoTransitionTexture tex(m_heap, m_device, m_commandList, m_texwidth, m_texheight, DXGI_FORMAT_R16G16B16A16_FLOAT);
        wchar_t name[30] = L"m_RWBurstCachedTex";
        tex.GetTextureData().texture->SetName(name);
        m_BurstCachedTex = tex;
    }


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ByteAddressBuffer a(m_device, m_heap);
    m_RandomIndexCounter.push_back(a);
    a.m_buffer->SetName(L"randomIndex");
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////以下はシェーダ側の設定

    {
        D3D12_DESCRIPTOR_RANGE T0{}, T1{ }, U0{}, U1{}, U2{}, U3{}, U4{};
        T0.RangeType = D3D12_DESCRIPTOR_RANGE_TYPE_SRV;
        T0.NumDescriptors = 1;//もともと
        T0.BaseShaderRegister = 0;

        //追加
        T1.RangeType = D3D12_DESCRIPTOR_RANGE_TYPE_SRV;
        T1.NumDescriptors = 1;//もともと
        T1.BaseShaderRegister = 1;


        U0.RangeType = D3D12_DESCRIPTOR_RANGE_TYPE_UAV;
        U0.NumDescriptors = 1;//もともと
        U0.BaseShaderRegister = 0;

        //追加
        U1.RangeType = D3D12_DESCRIPTOR_RANGE_TYPE_UAV;
        U1.NumDescriptors = 1;//もともと
        U1.BaseShaderRegister = 1;

        //追加
        U2.RangeType = D3D12_DESCRIPTOR_RANGE_TYPE_UAV;
        U2.NumDescriptors = 1;//もともと
        U2.BaseShaderRegister = 2;

        //追加
        U3.RangeType = D3D12_DESCRIPTOR_RANGE_TYPE_UAV;
        U3.NumDescriptors = 1;//もともと
        U3.BaseShaderRegister = 3;

        //追加
        U4.RangeType = D3D12_DESCRIPTOR_RANGE_TYPE_UAV;
        U4.NumDescriptors = 1;//もともと
        U4.BaseShaderRegister = 4;

        array<CD3DX12_ROOT_PARAMETER, 10> rootParams;//定数バッファ用にもう一つ追加したほうが良い
        rootParams[0].InitAsConstantBufferView(0);//ConstantBuffer<ComputeParameters> computeConstants : register(b0);
        rootParams[1].InitAsDescriptorTable(1, &T0);//Texture2D<float4> sourceImageR : register(t0);
        rootParams[2].InitAsDescriptorTable(1, &T1);//Texture2D<float4> sourceImageI : register(t1);
        rootParams[3].InitAsDescriptorTable(1, &U0); //RWTexture2D<float4> destinationImageR : register(u0);
        rootParams[4].InitAsDescriptorTable(1, &U1);//RWTexture2D<float4> destinationImageI : register(u1);
        rootParams[5].InitAsDescriptorTable(1, &U2); //RWTexture2D<float4> destinationImageR1 : register(u2);
        rootParams[6].InitAsDescriptorTable(1, &U3);//RWTexture2D<float4> destinationImageI1 : register(u3);
        rootParams[7].InitAsConstantBufferView(1);//ConstantBuffer<RandomTbl> randomTbl : register(b1);
        rootParams[8].InitAsDescriptorTable(1, &U4);//RWByteAddressBuffer randomTblIndex : register(u4);
        rootParams[9].InitAsConstantBufferView(2);//ConstantBuffer<GhostRotateMatrix> ghostGotateMat : register(b2);

        CD3DX12_ROOT_SIGNATURE_DESC rootSignatureDesc{};
        rootSignatureDesc.Init(
            UINT(rootParams.size()), rootParams.data(),
            0, nullptr,
            D3D12_ROOT_SIGNATURE_FLAG_ALLOW_INPUT_ASSEMBLER_INPUT_LAYOUT);

        ComPtr<ID3DBlob> signature, errBlob;
        hr = D3D12SerializeRootSignature(&rootSignatureDesc,
            D3D_ROOT_SIGNATURE_VERSION_1_0, &signature, &errBlob);
        if (FAILED(hr))
        {
            auto error = reinterpret_cast<const char*>(errBlob->GetBufferPointer());
            OutputDebugStringA(error);
        }
        hr = m_device->CreateRootSignature(0, signature->GetBufferPointer(), signature->GetBufferSize(), IID_PPV_ARGS(&m_csSignature));
        ThrowIfFailed(hr, "CreateRootSignature failed.");
    }

    vector<pair<wstring, string>> nameAtPipeline;
    getNameAtShaderAndPipeline(nameAtPipeline);

    std::vector<std::wstring> flags;
    std::vector<Shader::DefineMacro> defines;

    for (auto& p : nameAtPipeline)
    {
        //動的にヘッダ生成
        if (p.second == "fftCS_ROW")
        {
            createHeaderHLSL(m_texwidth, m_texheight, true, false);
        }
        else if (p.second == "fftCS_COL")
        {
            createHeaderHLSL(m_texwidth, m_texheight, false, false);
        }
        else if (p.second == "ifftCS_ROW")
        {
            createHeaderHLSL(m_texwidth, m_texheight, true, true);
        }
        else
        {
            createHeaderHLSL(m_texwidth, m_texheight, false, true);
        }

        Shader shader;
        shader.load(L"ComputeFilter.hlsl", Shader::Compute, p.first, flags, defines);//エントリポイント
        D3D12_COMPUTE_PIPELINE_STATE_DESC computeDesc{};
        computeDesc.CS = CD3DX12_SHADER_BYTECODE(shader.getCode().Get());
        computeDesc.pRootSignature = m_csSignature.Get();

        PipelineState pipeline;
        hr = m_device->CreateComputePipelineState(&computeDesc, IID_PPV_ARGS(&pipeline));
        ThrowIfFailed(hr, "CreateComputePipelineState failed.");
        m_pipelines[p.second] = pipeline;
    }
}

void ComputeFilterApp::getNameAtShaderAndPipeline(vector<pair<wstring, string>>& name)
{
    name.emplace_back(make_pair(L"mainSepia", "sepiaCS"));
    name.emplace_back(make_pair(L"mainSobel", "sobelCS"));
    name.emplace_back(make_pair(L"swap", "swapCS"));
    name.emplace_back(make_pair(L"mainMAXMIN", "maxminCS"));
    name.emplace_back(make_pair(L"mainNORMALIZE", "normalizeCS"));
    name.emplace_back(make_pair(L"mainMULTIPLY", "mulCS"));
    name.emplace_back(make_pair(L"mainFFT", "fftCS_ROW"));
    name.emplace_back(make_pair(L"mainFFT", "fftCS_COL"));
    name.emplace_back(make_pair(L"mainFFT", "ifftCS_ROW"));
    name.emplace_back(make_pair(L"mainFFT", "ifftCS_COL"));
    name.emplace_back(make_pair(L"mainDrawFRF", "frfCS0"));
    name.emplace_back(make_pair(L"mainColorScaling", "colorScaleCS"));
    name.emplace_back(make_pair(L"mainREAL", "realCS"));
    name.emplace_back(make_pair(L"mainIMAGE", "imageCS"));
    name.emplace_back(make_pair(L"mainAMPLITUDE", "ampCS"));
    name.emplace_back(make_pair(L"mainINTENSITY", "intenCS"));
    name.emplace_back(make_pair(L"mainPHASE", "phaseCS"));
    name.emplace_back(make_pair(L"mainDivByMaxAMP", "divByMaxAMPCS"));
    name.emplace_back(make_pair(L"mainREAL_IMAGE4Disp", "REALIMAGE4DISPCS"));
    name.emplace_back(make_pair(L"mainAMPLITUDE_INTENSITY4Disp", "AMPINT4DISPCS"));
    name.emplace_back(make_pair(L"mainPHASE4Disp", "PHASE4DISPCS"));
    name.emplace_back(make_pair(L"mainAdd", "AddCS"));
    name.emplace_back(make_pair(L"mainBinaryThreshold", "BTCS"));
    name.emplace_back(make_pair(L"mainDrawGaussian", "gaussianCS"));
    name.emplace_back(make_pair(L"mainDrawGaussianNoNormalize", "gaussianNoNormalizeCS"));
    name.emplace_back(make_pair(L"mainDrawGaussianNoNormalize2", "gaussianNoNormalize2CS"));
    name.emplace_back(make_pair(L"mainCopy", "copyCS"));
    name.emplace_back(make_pair(L"mainClear", "clearCS"));
    name.emplace_back(make_pair(L"mainRaiseBottom", "raiseCS"));
    name.emplace_back(make_pair(L"mainSpectrumScaling", "spectrumScalingCS"));
    name.emplace_back(make_pair(L"mainRaiseBottomRealImage", "raiseRICS"));
    name.emplace_back(make_pair(L"mainMAXMINfirst", "maxminfirstCS"));
    name.emplace_back(make_pair(L"mainMAXMINsecond", "maxminsecondCS"));
    name.emplace_back(make_pair(L"mainRadialBlur", "radialblurCS"));
    name.emplace_back(make_pair(L"mainDrawDecayLine", "decaylineCS"));
    name.emplace_back(make_pair(L"mainDrawQuadraticPhase", "quadraticCS0"));
    name.emplace_back(make_pair(L"mainDrawPolygon", "drawpolygonCS"));
    name.emplace_back(make_pair(L"mainDrawPolygonFixScale", "drawfixpolygonCS"));
    name.emplace_back(make_pair(L"mainDrawMovingPolygon", "drawmovingpolygonCS"));
    name.emplace_back(make_pair(L"mainSpectrumAmplitudeAdjustment", "amplitudeadjustmentCS"));
    name.emplace_back(make_pair(L"mainToneMapping", "tonemappingCS"));
    name.emplace_back(make_pair(L"mainDrawMovingDot", "drawMovingDotCS"));
    name.emplace_back(make_pair(L"mainDrawDotByRandomTbl", "drawDotByRandomTblCS"));
    name.emplace_back(make_pair(L"mainIncrementRandomTblIndex", "incrementRandomIndexCS"));
    name.emplace_back(make_pair(L"mainResetRandomTblIndex", "resetRandomIndexCS"));
    name.emplace_back(make_pair(L"mainScalingSizeByRandomTbl", "scalingSizeByRandomTblCS"));
    name.emplace_back(make_pair(L"mainShiftImageByRandomTbl", "shiftImageByRandomTblCS"));
    name.emplace_back(make_pair(L"mainShiftImageByRandomTblandAdd", "shiftImageByRandomTblansAddCS"));
    name.emplace_back(make_pair(L"mainShiftImageByTargetPos", "shiftImageByTargetPosCS"));
    name.emplace_back(make_pair(L"mainRotateByRandomTbl", "rotateImageByRandomTblCS"));
    name.emplace_back(make_pair(L"mainInverseRotateByRandomTbl", "invRotateImageByRandomTblCS"));
    name.emplace_back(make_pair(L"mainApplyVignettingByRandomTbl", "applyVignettingByRandomTblCS"));
    name.emplace_back(make_pair(L"mainDrawSmallCirc", "drawSmallCircCS"));
    name.emplace_back(make_pair(L"mainCaustic", "causticCS"));
    name.emplace_back(make_pair(L"mainBlur", "blurCS"));
    name.emplace_back(make_pair(L"mainFadeByGauss", "fadeByGaussCS"));
    name.emplace_back(make_pair(L"mainScaleByPos", "scaleByPosCS"));
}

ComputeFilterApp::TextureData ComputeFilterApp::LoadTextureFromFile(const std::wstring& name, ImageSize& size)//テクスチャのロード
{
    DirectX::TexMetadata metadata;
    DirectX::ScratchImage image;

    HRESULT hr;
    if (name.find(L".png") != std::wstring::npos)
    {
        hr = DirectX::LoadFromWICFile(name.c_str(), 0, &metadata, image);
    }
    if (name.find(L".dds") != std::wstring::npos)
    {
        hr = DirectX::LoadFromDDSFile(name.c_str(), 0, &metadata, image);
    }
    if (name.find(L".tga") != std::wstring::npos)
    {
        hr = DirectX::LoadFromTGAFile(name.c_str(), &metadata, image);
    }

    size.w = metadata.width;
    size.h = metadata.height;

    //if (m_texwidth < metadata.width)
        m_texwidth = metadata.width;
   // if (m_texheight < metadata.height)
        m_texheight = metadata.height;

    D3D12_RESOURCE_FLAGS resFlags = D3D12_RESOURCE_FLAG_NONE;
    resFlags |= D3D12_RESOURCE_FLAG_ALLOW_UNORDERED_ACCESS;
    ComPtr<ID3D12Resource> texture;
    CreateTextureEx(m_device.Get(), metadata, resFlags, false, &texture);

    Buffer srcBuffer;
    std::vector<D3D12_SUBRESOURCE_DATA> subresources;
    PrepareUpload(m_device.Get(), image.GetImages(), image.GetImageCount(), metadata, subresources);
    const auto totalBytes = GetRequiredIntermediateSize(texture.Get(), 0, UINT(subresources.size()));

    auto staging = CreateResource(CD3DX12_RESOURCE_DESC::Buffer(totalBytes), D3D12_RESOURCE_STATE_GENERIC_READ, nullptr, D3D12_HEAP_TYPE_UPLOAD);
    auto command = CreateCommandList();
    UpdateSubresources(command.Get(),
        texture.Get(), staging.Get(), 0, 0, UINT(subresources.size()), subresources.data());
    auto barrier = CD3DX12_RESOURCE_BARRIER::Transition(
        texture.Get(),
        D3D12_RESOURCE_STATE_COPY_DEST, D3D12_RESOURCE_STATE_PIXEL_SHADER_RESOURCE);
    command->ResourceBarrier(1, &barrier);
    FinishCommandList(command);

    TextureData texData;
    texture.As(&texData.texture);
    texData.handleRead = m_heap->Alloc();

    D3D12_SHADER_RESOURCE_VIEW_DESC srvDesc{};
    srvDesc.Format = metadata.format;
    srvDesc.TextureCube.MipLevels = UINT(metadata.mipLevels);
    srvDesc.TextureCube.MostDetailedMip = 0;
    srvDesc.TextureCube.ResourceMinLODClamp = 0;
    srvDesc.Shader4ComponentMapping = D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING;
    srvDesc.ViewDimension = D3D12_SRV_DIMENSION_TEXTURE2D;
    m_device->CreateShaderResourceView(texture.Get(), &srvDesc, texData.handleRead);

    return texData;
}