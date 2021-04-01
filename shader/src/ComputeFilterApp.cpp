#include "../include/ComputeFilterApp.h"
#include <random>

#include "imgui.h"
#include "examples/imgui_impl_dx12.h"
#include "examples/imgui_impl_win32.h"

#include <algorithm>

using namespace std;
using namespace DirectX;

const int ComputeFilterApp::m_ghostNum = 40;
//const int ComputeFilterApp::m_ghostNum = 4;

ComputeFilterApp::ComputeFilterApp()
{
    m_propdistance.resize(6);

    m_focus.resize(1);

    for (int i = 0; i < m_focus.size(); i++)
    {
        m_propdistance[i] = 10;
        m_focus[i] = 10;
    }

    m_mainComputeCB.resize(0);

    m_textureNames.resize(4);

    string ext = ".png";

    //m_textureNames.at(0) = "SKY" + ext;
    //m_textureNames.at(0).resize(MAX_PATH);
    //m_textureNames.at(1) = "eyerash" + ext;
    //m_textureNames.at(1).resize(MAX_PATH);
    //m_textureNames.at(2) = "zp" + ext;
    //m_textureNames.at(2).resize(MAX_PATH);
    //m_textureNames.at(3) = "dust" + ext;
    //m_textureNames.at(3).resize(MAX_PATH);

    m_textureNames.at(0) = "SKY1024" + ext;
    m_textureNames.at(0).resize(MAX_PATH);
    m_textureNames.at(1) = "eyerash1024" + ext;
    m_textureNames.at(1).resize(MAX_PATH);
    m_textureNames.at(2) = "zp1024" + ext;
    m_textureNames.at(2).resize(MAX_PATH);
    m_textureNames.at(3) = "dust1024" + ext;
    m_textureNames.at(3).resize(MAX_PATH);

    {
        randomGenerate();
    }
}

int ComputeFilterApp::inv2pow(int n) {
    int m = 0;
    m = static_cast<int>(log(static_cast<double>(n)) / log(2.0));
    return m;
}

void ComputeFilterApp::Cleanup()
{
  
}

void ComputeFilterApp::PrepareSimpleModel()//ï`âÊÇ∑ÇÈÉÇÉfÉãÇçÏê¨
{
    using VertexData = std::vector<Vertex>;
    using IndexData = std::vector<UINT>;
    VertexData quadVertices;
    IndexData  quadIndices;

    float offset = 10.0f;

    float w = 135.0f * 2 * m_texwidth / m_texheight * 0.5  * 2;
    float h = 135.0f * 2 * 0.5 * 2;

    /*  quadVertices = {
        { XMFLOAT3(-480.0f - offset,  135.0f, 0.0f), XMFLOAT2(0.0f, 0.0f) },
        { XMFLOAT3(0.0f - offset,  135.0f, 0.0f), XMFLOAT2(1.0f, 0.0f) },
        { XMFLOAT3(-480.0f - offset, -135.0f, 0.0f), XMFLOAT2(0.0f, 1.0f) },
        { XMFLOAT3(0.0f - offset, -135.0f, 0.0f), XMFLOAT2(1.0f, 1.0f) },
      };*/
    quadVertices = {
     { XMFLOAT3(-w - offset,  h, 0.0f), XMFLOAT2(0.0f, 0.0f) },
     { XMFLOAT3(0.0f - offset,  h, 0.0f), XMFLOAT2(1.0f, 0.0f) },
     { XMFLOAT3(-w - offset, -h, 0.0f), XMFLOAT2(0.0f, 1.0f) },
     { XMFLOAT3(0.0f - offset, -h, 0.0f), XMFLOAT2(1.0f, 1.0f) },
    };

    quadVertices = {
 { XMFLOAT3(0,  0, 0.0f), XMFLOAT2(0.0f, 0.0f) },
 { XMFLOAT3(0.0f - 0,  0, 0.0f), XMFLOAT2(1.0f, 0.0f) },
 { XMFLOAT3(-0 - 0, -0, 0.0f), XMFLOAT2(0.0f, 1.0f) },
 { XMFLOAT3(0.0f - 0, -0, 0.0f), XMFLOAT2(1.0f, 1.0f) },
    };
    quadIndices = { 0, 1, 2, 3, };

    m_quad = CreateSimpleModel(quadVertices, quadIndices);

    /* quadVertices = {
       { XMFLOAT3(0.0f + offset,  135.0f, 0.0f), XMFLOAT2(0.0f, 0.0f) },
       { XMFLOAT3(480.0f + offset,  135.0f, 0.0f), XMFLOAT2(1.0f, 0.0f) },
       { XMFLOAT3(0.0f + offset, -135.0f, 0.0f), XMFLOAT2(0.0f, 1.0f) },
       { XMFLOAT3(480.0f + offset, -135.0f, 0.0f), XMFLOAT2(1.0f, 1.0f) },
     };*/
    quadVertices = {
      { XMFLOAT3(-w/2 + 0,  h, 0.0f), XMFLOAT2(0.0f, 0.0f) },
      { XMFLOAT3(w/2 + 0,  h, 0.0f), XMFLOAT2(1.0f, 0.0f) },
      { XMFLOAT3(-w/2 + 0, -h, 0.0f), XMFLOAT2(0.0f, 1.0f) },
      { XMFLOAT3(w/2 + 0, -h, 0.0f), XMFLOAT2(1.0f, 1.0f) },
    };
    m_quad2 = CreateSimpleModel(quadVertices, quadIndices);
}

void ComputeFilterApp::CopyTextute(AutoTransitionTexture& Dst, AutoTransitionTexture& Src)
{
    //if (Dst.getSourceDestState() == 0)//SOURCEèÛë‘  Ç»ÇÁDESTÇ÷
    //{
    //    Dst.TransitionToDestState();
    //}

    Dst.TransitionToDestState();

    //if (Src.getSourceDestState() == 1)//DESTèÛë‘  Ç»ÇÁSOURCEÇ÷
    //{
    //    Src.TransitionToSourceState();
    //}
    Src.TransitionToSourceState();

    m_commandList.Get()->CopyResource(Dst.GetTextureData().texture.Get(), Src.GetTextureData().texture.Get());
}

//ëÄçÏån
void ComputeFilterApp::OnMouseButtonDown(UINT msg)
{
    auto io = ImGui::GetIO();
    if (io.WantCaptureMouse)
    {
        return;
    }

    m_IsDragged = true;
    m_ButtonType = msg;
}
void ComputeFilterApp::OnMouseButtonUp(UINT msg) {
    m_IsDragged = false;
}

void ComputeFilterApp::OnMouseMove(UINT msg, int dx, int dy)
{
    auto io = ImGui::GetIO();
    if (io.WantCaptureMouse)
    {
        return;
    }
   // m_camera.OnMouseMove(dx, dy);

    if (!m_IsDragged)
    {
        return;
    }

    if (m_ButtonType == 0)
    {
        m_posx += dx * 0.001f;
        m_posy += dy * 0.001f * m_width / m_height;

        m_posx = clamp(m_posx, -0.5f, 1.5f);
        m_posy = clamp(m_posy, -0.5f, 1.5f);
    }
}