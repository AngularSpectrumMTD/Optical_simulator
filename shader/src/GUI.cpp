#define _CRT_SECURE_NO_WARNINGS

#include "../include/ComputeFilterApp.h"

#include "imgui.h"
#include "examples/imgui_impl_dx12.h"
#include "examples/imgui_impl_win32.h"

void ComputeFilterApp::RenderHUD()//UI‚Ì•`‰æ
{
	ImGui_ImplDX12_NewFrame();
	ImGui_ImplWin32_NewFrame();
	ImGui::NewFrame();

	auto framerate = ImGui::GetIO().Framerate;
	ImGui::Begin("Information");
	ImGui::Text("Framerate Entire: %.3f fps", framerate);
	ImGui::Text("DeltaTime Entire: %.3f ms", 1000.0f / framerate);
	if (ImGui::Button("Shader Recompile")) {
		m_recompile = true;
		m_ghostKernelRegenerate = true;
		m_burstKernelRegenerate = true;
	}
	ImGui::Combo("Filter", (int*)&m_mode,
		"Sepia Filter\0Sobel Filter\0Swap\0FFT\0IFFT\0Convolution\0AsmProp_Real\0AsmProp_Image"
		"\0AsmProp_Amplitude\0AsmProp_Intensity\0AsmProp_Phase\0Kamehameha\0Gaussian_Blur\0Glare\0GodRay\0Ghost\0\0");
	ImGui::Spacing();

	if (ImGui::CollapsingHeader("Images[.png]"))
	{
		ImGui::InputText("input image", m_textureNames[0].data(), MAX_PATH);
		ImGui::InputText("eyerash image", m_textureNames[1].data(), MAX_PATH);
		ImGui::InputText("zoneplate image", m_textureNames[2].data(), MAX_PATH);
		ImGui::InputText("dust image", m_textureNames[3].data(), MAX_PATH);
	}

    ImGui::Spacing();

	if (!(m_mode == Mode_Sepia
		|| m_mode == Mode_Sobel
		|| m_mode == Mode_swap
		|| m_mode == Mode_FFT
		|| m_mode == Mode_IFFT
		|| m_mode == Mode_CONVOLUTION))
	{
		displayParameters();
	}

	ImGui::End();

	switch (m_mode)
	{
	case Mode_CONVOLUTION:
		displayRImageImGui(1, "ConvolutionKernel");
		break;
	case Mode_KAMEHAMEHA:
		displayRWImageImGui(0, "OverlaidImage");
		displayRWImageImGui(1, "ApertureImage");
		break;
	case Mode_GAUSSIAN_BLUR:
		displayRWImageImGui(0, "OverlaidImage");
		break;
	case Mode_GLARE:
	{
		switch (m_glaremode)
		{
		case GlareMode_EyeRush:
			displayRImageImGui(1, "EyeRushImage");
			break;
		case GlareMode_Cross:
			displayRWImageImGui(2, "ApertureBladesImage");
			break;
		case GlareMode_Halo:
			displayRImageImGui(2, "ZonePlateImage");
			break;
		}
		displayRWImageImGui(0, "OverlaidImage");
		displayRWImageImGui(1, "GlareImage");
	}
	break;
	case Mode_GODRAY:
		displayRWImageImGui(0, "OverlaidImage");
		break;
	case Mode_GHOST:
		displayRWImageImGui(0, "OverrideImage");
		displayRWImageImGui(1, "ApertureMaskedImage");
		displayRWImageImGui(2, "ApertureImage");
	default:
		break;
	}

	ImGui::Render();
	ImGui_ImplDX12_RenderDrawData(ImGui::GetDrawData(), m_commandList.Get());
}

void ComputeFilterApp::displayParameters()
{
	if (ImGui::CollapsingHeader("Parameters"))
	{
		switch (m_mode)
		{
		case Mode_PROP_REAL:
		case Mode_PROP_IMAGE:
		case Mode_PROP_AMPLITUDE:
		case Mode_PROP_INTENSITY:
		case Mode_PROP_PHASE:
			ImGui::RadioButton("Input Is Image", &m_inputfrag, InputFrag_Image);
			ImGui::RadioButton("Input Is Aperture", &m_inputfrag, InputFrag_Aperture);
			ImGui::SliderInt("apertureNum", &m_N, 3, 20);
			ImGui::SliderFloat("apertureradius", &m_r, 0.0f, 1.0f);
			ImGui::SliderFloat("intervalx [um]", &m_interval_x, 0.1f, 10.0f);
			ImGui::SliderFloat("intensity", &m_glareintensity, 1.0f, 50.0f);
			ImGui::SliderFloat("propdistance [mm]", &m_propdistance[0], 0.0f, 40.0f);
			ImGui::RadioButton("FocusON", &m_focusfrag, FocusFrag_ON);
			ImGui::RadioButton("FocusOFF", &m_focusfrag, FocusFrag_OFF);
			ImGui::SliderFloat("focus [mm]", &m_focus[0], 0.0f, 40.0f);
			break;
		case Mode_KAMEHAMEHA:
			ImGui::SliderFloat("intervalx [um]", &m_interval_x, 0.1f, 10.0f);
			ImGui::SliderFloat("propdistance [mm]", &m_propdistance[0], 0.0f, 10.0f);
			ImGui::SliderInt("INapertureNum", &m_N1, 3, 20);
			ImGui::SliderFloat("INapertureradius", &m_r1, 0.0f, 1.0f);
			ImGui::SliderFloat("posx", &m_posx, 0.0f, 0.99f);
			ImGui::SliderFloat("posy", &m_posy, 0.0f, 0.99f);
			ImGui::RadioButton("FocusON", &m_focusfrag, FocusFrag_ON);
			ImGui::RadioButton("FocusOFF", &m_focusfrag, FocusFrag_OFF);
			ImGui::SliderFloat("focus [mm]", &m_focus[0], 0.0f, 10.0f);
			break;
		case Mode_GAUSSIAN_BLUR:
			ImGui::SliderFloat("intervalx [um]", &m_interval_x, 0.1f, 10.0f);
			ImGui::SliderFloat("sigma", &m_gausssigma, 0.1f, 10.0f);
			ImGui::SliderFloat("threshold", &m_threshold, 0.85f, 1.00f);
			break;
		case Mode_GLARE:
			ImGui::SliderInt("lambdasamplenum", &m_glarelambdasamplenum, 1, 20);
			ImGui::RadioButton("EyeRush", &m_glaremode, GlareMode_EyeRush);
			ImGui::RadioButton("Cross", &m_glaremode, GlareMode_Cross);
			ImGui::SliderInt("apertureNum", &m_N, 3, 20);
			ImGui::SliderFloat("apertureradius", &m_r, 0.0f, 1.0f);
			ImGui::RadioButton("ApertureDust_OFF", &m_aperturedust, ApertureDust_OFF);
			ImGui::RadioButton("ApertureDust_ON", &m_aperturedust, ApertureDust_ON);
			ImGui::RadioButton("Halo", &m_glaremode, GlareMode_Halo);
			ImGui::SliderFloat("intensity", &m_glareintensity, 0.1f, 50.0f);
			ImGui::SliderFloat("threshold", &m_threshold, 0.85f, 1.00f);
			break;
		case Mode_GODRAY:
			ImGui::SliderFloat("RayLength", &m_raylength, 1.0f, 10.0f);
			ImGui::SliderFloat("posx", &m_posx, 0.0f, 0.99f);
			ImGui::SliderFloat("posy", &m_posy, 0.0f, 0.99f);
			ImGui::SliderFloat("threshold", &m_threshold, 0.85f, 1.00f);
			ImGui::SliderInt("quality", &m_godrayquarity, 1, 10);
			ImGui::SliderFloat("sigma", &m_gausssigma, 0.1f, 10.0f);
			break;
		case Mode_GHOST:

			ImGui::RadioButton("Input Is Image", &m_inputfrag, InputFrag_Image);
			ImGui::RadioButton("Input Is Non", &m_inputfrag, InputFrag_Aperture);

			if (ImGui::SliderFloat("randSigmaOpposite(Near)", &m_randSigmaOppositeNear, 0.01f, 30.0f))
			{
				m_randomgenerate = true;
			}
			if (ImGui::SliderFloat("randSigmaOpposite(Far)", &m_randSigmaOppositeFar, 0.01f, 30.0f))
			{
				m_randomgenerate = true;
			}
			if (ImGui::SliderFloat("randSigma", &m_randSigma, 0.01f, 30.0f))
			{
				m_randomgenerate = true;
			}
			if (ImGui::Button("Random Regenerate")) {
				m_randomgenerate = true;
			}

			ImGui::RadioButton("KerareON(More Accurate)", &m_kerare_ON, 1);
			ImGui::RadioButton("KerareOFF", &m_kerare_ON, 0);

			ImGui::ColorPicker3("ghostBaseColor", &m_baseColor.w_x);
			//m_baseColor /= 255.99f;

			ImGui::SliderFloat("kerareRange", &m_gausssigma2, 0.0f, 50.0f);

			//ImGui::SliderFloat("ghostScale", &m_ghostScale, 0.1f, 50.0f);
			
			/*


			ImGui::SliderFloat("posx", &m_posx, 0.0f, 0.99f);
			ImGui::SliderFloat("posy", &m_posy, 0.0f, 0.99f);*/

			if (ImGui::SliderFloat("burstIntensity", &m_glareintensity, 0.1f, 200.0f))
			{
				m_burstKernelRegenerate = true;
			}

			if (ImGui::SliderFloat("propdistance [mm]", &m_propdistance[0], 0.0f, 10.0f))
			{
				m_ghostKernelRegenerate = true;
			}
			if (ImGui::SliderInt("apertureNum", &m_N, 5, 20))
			{
				m_ghostKernelRegenerate = true;
				m_burstKernelRegenerate = true;
			}
			if (ImGui::SliderFloat("apertureradius", &m_r, 0.1f, 1.0f))
			{
				m_ghostKernelRegenerate = true;
				m_burstKernelRegenerate = true;
			}
			if (ImGui::SliderFloat("rotAngle", &m_rotAngle, 0.0f, 360.0f))
			{
				m_ghostKernelRegenerate = true;
				m_burstKernelRegenerate = true;
			}

		
			break;
		default:
			break;
		}
	}
}

void ComputeFilterApp::displayRImageImGui(UINT index, const char* tag)
{
	if (tag == nullptr)
	{
		tag = "image";
	}
	ImGui::SetNextWindowSize(ImVec2(m_texwidth / 2, m_texheight / 2), ImGuiCond_Once);
	ImGui::Begin(tag, NULL, ImGuiWindowFlags_NoResize);
	{
		ImVec2 avail_size = ImGui::GetContentRegionAvail();
		ImGui::Image((void*)D3D12_GPU_DESCRIPTOR_HANDLE(m_RfullsizeTex.at(index).handleRead).ptr, avail_size);
	}
	ImGui::End();
}

void ComputeFilterApp::displayRWImageImGui(UINT index, const char* tag)
{
	if (tag == nullptr)
	{
		tag = "image";
	}
	ImGui::SetNextWindowSize(ImVec2(m_texwidth / 2, m_texheight / 2), ImGuiCond_Once);
	ImGui::Begin(tag, NULL, ImGuiWindowFlags_NoResize);
	{
		ImVec2 avail_size = ImGui::GetContentRegionAvail();
		ImGui::Image((void*)D3D12_GPU_DESCRIPTOR_HANDLE(m_RWdisplayTex.at(index).GetTextureData().handleRead).ptr, avail_size);
	}
	ImGui::End();
}