#include"../include/WaveFront.h"
using namespace std;
void WaveFront::SaveAsWaveFront(const char* filename)
{
	ofstream file(filename);
	if (!file) {
		cerr << "Error: file not opened." << endl;
		exit(0);
	}

	file << GetNx() << " ";
	file << GetNy() << " ";
	file << GetPx() << " ";
	file << GetPy() << " ";
	file << GetLambda() << " ";
	file << "\n";

	file << GetOrigin().w_x << " ";
	file << GetOrigin().w_y << " ";
	file << GetOrigin().w_z << " ";
	file << "\n";

	int i = 0, j = 0;
	for (j = 0; j < w_ny; j++)
	{
		for (i = 0; i < w_nx; i++)
		{
			file << GetPixel(i, j) << "\t";
		}
		file << "\n";
	}
	file << "\n";
	file.close();
}

WaveFront& WaveFront::LoadAsWaveFront(const char* filename)
{
	vector<vector<complex<double>>> matrix;
	ifstream file(filename);
	if (!file) {
		cerr << "Error: file not opened." << endl;
		exit(0);
	}

	string line;
	getline(file, line);
	const char* buf;
	buf = line.c_str();

	WaveFront& ret = *this;

	int nx, ny;
	double px, py, lam;
	sscanf_s(buf, "%d %d %lf %lf %lf", &nx, &ny, &px, &py, &lam);

	ret.SetNx(nx); ret.SetNy(ny); ret.SetPx(px); ret.SetPy(py); ret.SetLambda(lam);

	ret.Init();

	getline(file, line);
	buf = line.c_str();
	float x, y, z;
	sscanf_s(buf, "%g %g %g", &x, &y, &z);

	vec3 origin {x, y, z};
	ret.SetOrigin(origin);

	int i = 0, j = 0;
	while (getline(file, line)) {
		istringstream stream(line);
		complex<double> c;
		vector<complex<double>> row;
		while (stream >> c) {
			//row.push_back(c);
			ret.SetPixel(i,j,c);
			i++;
		}
		i = 0;
		//matrix.push_back(row);
		j++;
	}
	file.close();
	return ret;
}