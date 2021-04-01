#include "../include/ComputeFilterApp.h"
#include <random>

using namespace std;
using namespace DirectX;

void ComputeFilterApp::Render()//�`��̑S��
{
    m_frameIndex = m_swapchain->GetCurrentBackBufferIndex();
    m_commandAllocators[m_frameIndex]->Reset();
    m_commandList->Reset(
        m_commandAllocators[m_frameIndex].Get(), nullptr
    );

    // �X���b�v�`�F�C���\���\���烌���_�[�^�[�Q�b�g�`��\��
    auto barrierToRT = m_swapchain->GetBarrierToRenderTarget();
    m_commandList->ResourceBarrier(1, &barrierToRT);

    ID3D12DescriptorHeap* heaps[] = { m_heap->GetHeap().Get() };
    m_commandList->SetDescriptorHeaps(_countof(heaps), heaps);

    RenderToMain();

    RenderHUD();

    // �����_�[�^�[�Q�b�g����X���b�v�`�F�C���\���\��
    {
        auto barrierToPresent = m_swapchain->GetBarrierToPresent();
        CD3DX12_RESOURCE_BARRIER barriers[] = {
          barrierToPresent,
        };

        m_commandList->ResourceBarrier(_countof(barriers), barriers);
    }

    m_commandList->Close();

    ID3D12CommandList* lists[] = { m_commandList.Get() };
    m_commandQueue->ExecuteCommandLists(1, lists);

    //���̂�����
    if (m_recompile == true)
    {
        WaitForIdleGPU();
        PrepareComputeFilter();
        PrepareConstantBuffers();
        PrepareSimpleModel();
        m_recompile = false;
    }
    if (m_randomgenerate)
    {
        randomGenerate();
    }

   // m_swapchain->Present(1, 0);
    m_swapchain->Present(0, 0);
    m_swapchain->WaitPreviousFrame(m_commandQueue, m_frameIndex, GpuWaitTimeout);
}

void ComputeFilterApp::Flush()
{
    auto sceneCB = m_mainSceneCB[m_frameIndex];
    auto mtxProj = XMMatrixOrthographicLH(720.0f, 720.0f, 0.0f, 10.0f);
    //auto mtxProj = XMMatrixOrthographicLH(1280.0f, 720.0f, 0.0f, 10.0f);

    SceneParameters sceneParams{};
    XMStoreFloat4x4(&sceneParams.proj, XMMatrixTranspose(mtxProj));
    WriteToUploadHeapMemory(sceneCB.Get(), sizeof(sceneParams), &sceneParams);

    //m_commandList->ResourceBarrier(1, &m_RWfullsizeTexWtoR.at(0));

    m_commandList->SetPipelineState(m_pipelines["default"].Get());
    m_commandList->SetGraphicsRootSignature(m_rootSignature.Get());
    m_commandList->IASetPrimitiveTopology(D3D_PRIMITIVE_TOPOLOGY_TRIANGLESTRIP);

    //�����ɕ\��������
    m_commandList->IASetIndexBuffer(&m_quad.ibView);
    m_commandList->IASetVertexBuffers(0, 1, &m_quad.vbView);
    m_commandList->SetGraphicsRootConstantBufferView(0, sceneCB->GetGPUVirtualAddress());
    m_commandList->SetGraphicsRootDescriptorTable(1, m_RfullsizeTex.at(0).handleRead);//�\������摜
    m_commandList->DrawIndexedInstanced(4, 1, 0, 0, 0);

    //�E���ɕ\��������
    m_commandList->IASetIndexBuffer(&m_quad2.ibView);
    m_commandList->IASetVertexBuffers(0, 1, &m_quad2.vbView);
    m_commandList->SetGraphicsRootConstantBufferView(0, sceneCB->GetGPUVirtualAddress());
    m_commandList->SetGraphicsRootDescriptorTable(1, m_RWfullsizeTex.at(0).ReadState());//�\������摜
    m_commandList->DrawIndexedInstanced(4, 1, 0, 0, 0);

    //m_commandList->ResourceBarrier(1, &m_RWfullsizeTexRtoW.at(0));
}

void ComputeFilterApp::Flush(int index)
{
    auto sceneCB = m_mainSceneCB[m_frameIndex];
    auto mtxProj = XMMatrixOrthographicLH(720.0f, 720.0f, 0.0f, 10.0f);

    SceneParameters sceneParams{};
    XMStoreFloat4x4(&sceneParams.proj, XMMatrixTranspose(mtxProj));
    WriteToUploadHeapMemory(sceneCB.Get(), sizeof(sceneParams), &sceneParams);

    m_commandList->SetPipelineState(m_pipelines["default"].Get());
    m_commandList->SetGraphicsRootSignature(m_rootSignature.Get());
    m_commandList->IASetPrimitiveTopology(D3D_PRIMITIVE_TOPOLOGY_TRIANGLESTRIP);

    //�����ɕ\��������
    m_commandList->IASetIndexBuffer(&m_quad.ibView);
    m_commandList->IASetVertexBuffers(0, 1, &m_quad.vbView);
    m_commandList->SetGraphicsRootConstantBufferView(0, sceneCB->GetGPUVirtualAddress());
    m_commandList->SetGraphicsRootDescriptorTable(1, m_RfullsizeTex.at(0).handleRead);//�\������摜
    m_commandList->DrawIndexedInstanced(4, 1, 0, 0, 0);

    //�E���ɕ\��������
    m_commandList->IASetIndexBuffer(&m_quad2.ibView);
    m_commandList->IASetVertexBuffers(0, 1, &m_quad2.vbView);
    m_commandList->SetGraphicsRootConstantBufferView(0, sceneCB->GetGPUVirtualAddress());
    m_commandList->SetGraphicsRootDescriptorTable(1, m_RWfullsizeTex.at(index).ReadState());//�\������摜
    m_commandList->DrawIndexedInstanced(4, 1, 0, 0, 0);
}

void ComputeFilterApp::RenderToMain()//���C���Ƃ��ĕ`���摜(����̓L�����̉摜�����݂�����)��`�� �����ŃV�F�[�_�̊֐����Ă�
{
    auto rtv = m_swapchain->GetCurrentRTV();
    auto dsv = m_defaultDepthDSV;

    float m_clearColor[4] = { 0.1f,0.5f,0.75f,0 };
    m_commandList->ClearRenderTargetView(rtv, m_clearColor, 0, nullptr);

    m_commandList->ClearDepthStencilView(
        dsv, D3D12_CLEAR_FLAG_DEPTH, 1.0f, 0, 0, nullptr);

    auto d3d12_cpu_desc_handle_rtv = (D3D12_CPU_DESCRIPTOR_HANDLE)rtv;
    auto d3d12_cpu_desc_handle_dsv = (D3D12_CPU_DESCRIPTOR_HANDLE)dsv;

    m_commandList->OMSetRenderTargets(1, &d3d12_cpu_desc_handle_rtv,
        FALSE, &d3d12_cpu_desc_handle_dsv);

    auto viewport = CD3DX12_VIEWPORT(0.0f, 0.0f, float(m_width), float(m_height));
    auto scissorRect = CD3DX12_RECT(0, 0, LONG(m_width), LONG(m_height));
    m_commandList->RSSetViewports(1, &viewport);
    m_commandList->RSSetScissorRects(1, &scissorRect);

    m_commandList->SetComputeRootSignature(m_csSignature.Get());

    switch (m_mode)
    {
    case Mode_FFT:
        Flush(ExecuteFFTCommand());
        break;
    case Mode_IFFT:
        Flush(ExecuteIFFTCommand());
        break;
    case Mode_CONVOLUTION:
        Flush(ExecuteConvolutionCommand());
        break;
    case Mode_PROP_REAL:
    case Mode_PROP_IMAGE:
    case Mode_PROP_AMPLITUDE:
    case Mode_PROP_INTENSITY:
    case Mode_PROP_PHASE:
        Flush(ExecuteASMCommand());
        break;
    case Mode_KAMEHAMEHA:
        Flush(ExecuteKamehamehaCommand());
        break;
    case Mode_GAUSSIAN_BLUR:
        Flush(ExecuteGaussianBlurCommand());
        break;
    case Mode_GLARE:
        Flush(ExecuteGlareCommand());
        break;
    case Mode_GODRAY:
        Flush(ExecuteGodRayCommand());
        break;
    case Mode_GHOST:
        Flush(ExecuteLensFlareCommand());
        break;
    default:
        Flush(ExecuteOtherCommand());
        break;
    }
}