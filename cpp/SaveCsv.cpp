#include"WaveFront.h"

void WaveFront::SaveAsCsv(const char* fname, Axis axis, int ij)
{
	FILE* fp = fopen(fname, "w+t");

	if (axis == X_AXIS)
	{
		if (ij < 0 || ij >= _ny)
		{
			printf("Index is invalid");
			return;
		}

		fprintf(fp, "%s,  %s,  %s,  %s,  %s,  %s\n",
			"X", "Real", "Imaginary", "Amplitude", "Phase", "Intensity");
		for (int i = 0; i < _nx; i++)
		{
			fprintf(fp, "%g,  %g,  %g,  %g,  %g,  %g\n", itox(i),
				GetReal(i, ij), GetImage(i, ij),
				GetAmplitude(i, ij), (double)GetPhase(i, ij), GetIntensity(i, ij));
		}
	}
	else
	{
		if (ij < 0 || ij >= _nx)
		{
			printf("Index is invalid");
			return;
		}

		fprintf(fp, "%s,  %s,  %s,  %s,  %s,  %s\n",
			"Y", "Real", "Imaginary", "Amplitude", "Phase", "Intensity");
		for (int j = 0; j < _ny; j++)
		{
			fprintf(fp, "%g,  %g,  %g,  %g,  %g,  %g\n", jtoy(j),
				GetReal(ij, j), GetImage(ij, j),
				GetAmplitude(ij, j), (double)GetPhase(ij, j), GetIntensity(ij, j));
		}
	}

	fclose(fp);

	return;
}