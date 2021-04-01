#include "../include/AutoTransitionTexture.h"

AutoTransitionTexture::AutoTransitionTexture(const std::shared_ptr<DescriptorManager> &heap, 
    const ComPtr<ID3D12Device> &device, 
    const ComPtr<ID3D12GraphicsCommandList> &commandlist, 
    UINT width, UINT height, DXGI_FORMAT format)
{
	m_commandList = commandlist;//コマンドリストの登録
    m_device = device;
    m_heap = heap;

    m_width = width;
    m_height = height;
    m_format = format;
    m_readWriteState = WRITE;
    m_sourceDestState = DEST;

	auto texDesc = CD3DX12_RESOURCE_DESC::Tex2D(format, width, height);
	texDesc.Flags |= D3D12_RESOURCE_FLAG_ALLOW_UNORDERED_ACCESS;

    //テクスチャデータの作成
	TextureData tex;

    auto cd3dx12_heap_properties = CD3DX12_HEAP_PROPERTIES(D3D12_HEAP_TYPE_DEFAULT);

    HRESULT hr = m_device->CreateCommittedResource(
        &cd3dx12_heap_properties,
        D3D12_HEAP_FLAG_NONE,
        &texDesc,
        D3D12_RESOURCE_STATE_UNORDERED_ACCESS,
        nullptr,
        IID_PPV_ARGS(&tex.texture)
    );
    ThrowIfFailed(hr, "CreateCommittedResource failed.");
    m_buffer = tex;

    //UAV・SRV設定
      //書き込み時の設定
    D3D12_UNORDERED_ACCESS_VIEW_DESC uavDesc{};
    uavDesc.ViewDimension = D3D12_UAV_DIMENSION_TEXTURE2D;
    uavDesc.Format = texDesc.Format;
    uavDesc.Texture2D.MipSlice = 0;
    uavDesc.Texture2D.PlaneSlice = 0;
    m_buffer.handleWrite = m_heap->Alloc();
    m_device->CreateUnorderedAccessView(
        m_buffer.texture.Get(),
        nullptr,
        &uavDesc,
        m_buffer.handleWrite
    );

    //読み出し時の設定
    D3D12_SHADER_RESOURCE_VIEW_DESC srvDesc{};
    srvDesc.Format = texDesc.Format;
    srvDesc.ViewDimension = D3D12_SRV_DIMENSION_TEXTURE2D;
    srvDesc.Texture2D.MipLevels = 1;
    srvDesc.Shader4ComponentMapping = D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING;
    m_buffer.handleRead = m_heap->Alloc();
    m_device->CreateShaderResourceView(
        m_buffer.texture.Get(),
        &srvDesc,
        m_buffer.handleRead
    );

    //バリアの作成
    auto barrierWtoR = CD3DX12_RESOURCE_BARRIER::Transition(
        m_buffer.texture.Get(),
        D3D12_RESOURCE_STATE_UNORDERED_ACCESS,
        D3D12_RESOURCE_STATE_PIXEL_SHADER_RESOURCE
    );
    m_barrierWtoR = barrierWtoR;

    auto barrierRtoW = CD3DX12_RESOURCE_BARRIER::Transition(
        m_buffer.texture.Get(),
        D3D12_RESOURCE_STATE_PIXEL_SHADER_RESOURCE,
        D3D12_RESOURCE_STATE_UNORDERED_ACCESS
    );
    m_barrierRtoW = barrierRtoW;

    auto barrierDtoS = CD3DX12_RESOURCE_BARRIER::Transition(
        m_buffer.texture.Get(),
        D3D12_RESOURCE_STATE_COPY_DEST,
        D3D12_RESOURCE_STATE_COPY_SOURCE
    );
    m_barrierDtoS = barrierDtoS;

    auto barrierStoD = CD3DX12_RESOURCE_BARRIER::Transition(
        m_buffer.texture.Get(),
        D3D12_RESOURCE_STATE_COPY_SOURCE,
        D3D12_RESOURCE_STATE_COPY_DEST
    );
    m_barrierStoD = barrierStoD;
}

DescriptorHandle AutoTransitionTexture::ReadState()
{
    if (m_readWriteState == WRITE)//W→R
    {
        m_readWriteState = READ;
        m_commandList->ResourceBarrier(1, &m_barrierWtoR);
    }
    return m_buffer.handleRead;
}

DescriptorHandle AutoTransitionTexture::WriteState()
{
    if (m_readWriteState == READ)//R→W
    {
        m_readWriteState = WRITE;
        m_commandList->ResourceBarrier(1, &m_barrierRtoW);
    }
    return m_buffer.handleWrite;
}

void AutoTransitionTexture::TransitionToSourceState()
{
    if (m_sourceDestState == DEST)//W→R
    {
        m_sourceDestState = SOURCE;
        m_commandList->ResourceBarrier(1, &m_barrierDtoS);
    }
}

void AutoTransitionTexture::TransitionToDestState()
{
    if (m_sourceDestState == SOURCE)//R→W
    {
        m_sourceDestState = DEST;
        m_commandList->ResourceBarrier(1, &m_barrierStoD);
    }
}

AutoTransitionTexture::SourceDestState AutoTransitionTexture::getSourceDestState()
{
    return m_sourceDestState;
}