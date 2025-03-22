/*
 *  This code is freely available under the following conditions:
 *
 *  1) The code is to be used only for non-commercial purposes.
 *  2) No changes and modifications to the code without prior permission of the developer.
 *  3) No forwarding the code to a third party without prior permission of the developer.
 *
 *                      FDEM_BfiPF
 * The file contains basic utility functions and simple operations
 *
 *  Written by D.Sc. Vagin D.V.
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 December 10, 2024 
*/

#include <stdio.h>
#include <stdlib.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <Windows.h>
#include <direct.h>

using namespace std;

bool isFileExists(char* fname)
{
	ifstream inf;
	bool flag;
	flag = false;
	inf.open(fname);
	if (inf)
	{
		flag = true;
		inf.close();
	}
	inf.clear();
	return flag;
}
//------------------------------------------------------     
// Calculation settings
//------------------------------------------------------     
struct Settings
{
public:
	double steps3D[2];    // Mesh initial steps
	double sparse3D[2];   // Sparce coefficient
	double farBound3D[2]; // Distance from receivers bounding box+gap to the boundary of the calculation domain
	double gapFromReceivers3D[2]; // Gap between receivers bounding box and sparced mesh area
	int receiversRefinementsCount; // Mesh refinements around receivers

	int numberOfABDipoles;

	double steps2D[2][2];    // Mesh initial steps
	double sparse2D[2][3];   // Sparce coefficient
	double farBound2D; // Distance from receivers bounding

	double Nu;	// Frequency

	int threadsCount; // Threads count for parallel calculation

	int met2d;	// Method for 2D (0 - U task, 1 - AV task)

	int met_slae_3d;	// SLAE solver for 3D (0 - pardiso, 1 - iteration)
	double eps3d;		// SLAE maximum residual for 3D if solver "iteration"
	int maxiter3d;		// SLAE maximum iteration for 3D if solver "iteration"

	Settings()
	{
		numberOfABDipoles = 100;
		Nu = 0.0;
		threadsCount = 1;
		met2d = 0;
		met_slae_3d = 0;
		eps3d = 1e-6;
		maxiter3d = 1000;
	}
};
//------------------------------------------------------     
// Electrode position
//------------------------------------------------------     
struct Electrode
{
public:

	double coordinates[3]; // x, y, z
};

FILE *log_file = NULL;
//------------------------------------------------------     
// Write message to log
//------------------------------------------------------     
void write_to_log(const char *str)
{
	printf("%s", str);
	if (log_file != NULL)
	{
		fprintf(log_file, "%s", str);
		fflush(log_file);
	}
}
//------------------------------------------------------     
// Open file for reading
//------------------------------------------------------     
int open_file_r(const char *file_name, FILE **file_stream)
{
	char buf[2048];

	if (!((*file_stream) = fopen(file_name, "r")))
	{
		sprintf(buf, "Error : Could not read file '%s'\n", file_name);
		write_to_log(buf);
		return 1;
	}

	return 0;
}
//------------------------------------------------------     
// Open file for writing
//------------------------------------------------------      
int open_file_w(const char *file_name, FILE **file_stream)
{
	char buf[2048];

	if (!((*file_stream) = fopen(file_name, "w")))
	{
		sprintf(buf, "Error : Could not write file '%s'\n", file_name);
		write_to_log(buf);
		return 1;
	}

	return 0;
}
//------------------------------------------------------     
// Open log file
//------------------------------------------------------     
int open_log(const char *file_name)
{
	if (open_file_w(file_name, &log_file) != 0)
	{
		printf("Error : could not open file '%s'\n", file_name);
		return 1;
	}

	return 0;
}

//------------------------------------------------------     
// Run Windows executable file
//------------------------------------------------------     
int ExecuteExe(const char *cmdline, const char *workdir)
{
	int retCode;
	STARTUPINFOA si;
	ZeroMemory(&si, sizeof(si));
	si.cb = sizeof(si);

	PROCESS_INFORMATION pi;
	ZeroMemory(&pi, sizeof(pi));

	CreateProcessA(NULL, (LPSTR)cmdline, NULL, NULL, FALSE, CREATE_NO_WINDOW, NULL, workdir, &si, &pi);
	WaitForSingleObject(pi.hProcess, INFINITE);
	GetExitCodeProcess(pi.hProcess, (LPDWORD)&retCode);

	CloseHandle(pi.hProcess);
	CloseHandle(pi.hThread);
	return retCode;
}

//------------------------------------------------------     
// Write single string file
//------------------------------------------------------     
int WriteSingleStringFile(const char *file_name, const char *str)
{
	FILE *file_out = NULL;
	if (open_file_w(file_name, &file_out) != 0) return 1;

	fprintf(file_out, "%s\n", str);

	fclose(file_out);
	return 0;

}
//------------------------------------------------------     
// Write single int file
//------------------------------------------------------     
int WriteSingleIntFile(const char *file_name, const int val)
{
	FILE *file_out = NULL;
	if (open_file_w(file_name, &file_out) != 0) return 1;

	fprintf(file_out, "%d\n", val);

	fclose(file_out);
	return 0;

}
//------------------------------------------------------     
// Write single string file
//------------------------------------------------------     
int WriteLcTxt(const char *file_name)
{
	FILE *file_out = NULL;
	if (open_file_w(file_name, &file_out) != 0) return 1;

	fprintf(file_out, "3\n");
	fprintf(file_out, "1e-3\n");
	fprintf(file_out, "1e-2\n");
	fprintf(file_out, "\n");
	fprintf(file_out, "3\n");
	fprintf(file_out, "1e-3\n");
	fprintf(file_out, "1e-2\n");
	fprintf(file_out, " \n");
	fprintf(file_out, "\n");
	fprintf(file_out, "\n");
	fprintf(file_out, "0\n");
	fprintf(file_out, "1e-3\n");
	fprintf(file_out, "1e-3\n");
	fprintf(file_out, "0.5\n");
	fprintf(file_out, "0.5\n");
	fprintf(file_out, "0.0\n");
	fprintf(file_out, "\n");
	fprintf(file_out, "0\n");
	fprintf(file_out, "1e-3\n");
	fprintf(file_out, "1e-3\n");
	fprintf(file_out, "0.5\n");
	fprintf(file_out, "0.5\n");
	fprintf(file_out, "0.5\n");
	fprintf(file_out, "\n");
	fprintf(file_out, "1\n");
	fprintf(file_out, "1e+30\n");
	fprintf(file_out, "1e-3\n");

	fclose(file_out);
	return 0;

}
//------------------------------------------------------     
// Write mesh settings file
//------------------------------------------------------     
int WriteMeshSettings(const char *file_name, Settings &settings)
{
	FILE *file_out = NULL;
	if (open_file_w(file_name, &file_out) != 0) return 1;
	char c;

	fprintf(file_out, "%.13e\t%.13e\n", settings.steps3D[0], settings.steps3D[1]);
	fprintf(file_out, "%.13e\t%.13e\n", settings.sparse3D[0], settings.sparse3D[1]);
	fprintf(file_out, "%.13e\t%.13e\n", settings.gapFromReceivers3D[0], settings.gapFromReceivers3D[1]);
	fprintf(file_out, "%.13e\t%.13e\n", settings.farBound3D[0], settings.farBound3D[1]);
	fprintf(file_out, "%d\n", settings.receiversRefinementsCount);

	fclose(file_out);
	return 0;
}
//------------------------------------------------------     
// Write mesh2d settings file
//------------------------------------------------------     
int WriteMeshSettings2D(const char* file_name, Settings& settings)
{
	FILE* file_out = NULL;
	if (open_file_w(file_name, &file_out) != 0) return 1;
	char c;

	fprintf(file_out, "%.13e\n", settings.farBound2D);
	fprintf(file_out, "%.13e\t%.13e\t", settings.steps2D[0][0], settings.steps2D[0][1]);
	fprintf(file_out, "%.13e\t%.13e\t%.13e\n", settings.sparse2D[0][0], settings.sparse2D[0][1], settings.sparse2D[0][2]);
	fprintf(file_out, "%.13e\t%.13e\t", settings.steps2D[1][0], settings.steps2D[1][1]);
	fprintf(file_out, "%.13e\t%.13e\t%.13e\n", settings.sparse2D[1][0], settings.sparse2D[1][1], settings.sparse2D[1][2]);

	fclose(file_out);
	return 0;
}
//------------------------------------------------------     
// Read application settings file
//------------------------------------------------------     
int ReadSettings(const char *file_name, Settings &settings)
{
	FILE *file_in = NULL;
	if (open_file_r(file_name, &file_in) != 0) return 1;

	fscanf(file_in, "%lf%lf", settings.steps3D, settings.steps3D + 1);
	fscanf(file_in, "%lf%lf", settings.sparse3D, settings.sparse3D + 1);
	fscanf(file_in, "%lf%lf", settings.gapFromReceivers3D, settings.gapFromReceivers3D + 1);
	fscanf(file_in, "%lf%lf", settings.farBound3D, settings.farBound3D + 1);
	fscanf(file_in, "%d", &settings.receiversRefinementsCount);

	//fscanf(file_in, "%d", &settings.numberOfMNPoints);
	fscanf(file_in, "%d", &settings.numberOfABDipoles);	

	fscanf(file_in, "%lf", &(settings.farBound2D));
	fscanf(file_in, "%lf%lf", &(settings.steps2D[0][0]), &(settings.steps2D[0][1]));
	fscanf(file_in, "%lf%lf%lf", &(settings.sparse2D[0][0]), &(settings.sparse2D[0][1]), &(settings.sparse2D[0][2]));
	fscanf(file_in, "%lf%lf", &(settings.steps2D[1][0]), &(settings.steps2D[1][1]));
	fscanf(file_in, "%lf%lf%lf", &(settings.sparse2D[1][0]), &(settings.sparse2D[1][1]), &(settings.sparse2D[1][2]));

	fscanf(file_in, "%lf", &settings.Nu);

	fscanf(file_in, "%d", &settings.threadsCount);

	fscanf(file_in, "%d", &settings.met2d);

	fscanf(file_in, "%d", &settings.met_slae_3d);
	fscanf(file_in, "%lf", &settings.eps3d);
	fscanf(file_in, "%d", &settings.maxiter3d);

	fclose(file_in);

	return 0;
}
//------------------------------------------------------     
// Prepare directory for meshing
//------------------------------------------------------     
int PrepareDirectoryForMeshing(const char *path, Settings &settings)
{
	char fileName[2048];
	ifstream inf;
	ofstream ofp;
	double z, ax, ay, bx, by;
	int i,nLay;
	vector<double> LayZ, LaySr, LaySz;

	sprintf(fileName, "%s/objects", path); CopyFileA("objects", fileName, false);
	//sprintf(fileName, "%s/z_sig_2d", path); CopyFileA("layers", fileName, false);
	
	inf.open("layers");
	if (!inf)
	{
		return 1;
	}
	inf >> nLay;
	LayZ.resize(nLay);
	LaySr.resize(nLay);
	LaySz.resize(nLay);
	for (i = 0; i < nLay; i++)
	{
		inf >> LayZ[i] >> LaySr[i] >> LaySz[i];
	}
	inf.close();
	inf.clear();

	sprintf(fileName, "%s/z_sig_2d", path);
	ofp.open(fileName);
	ofp << scientific << setprecision(14);
	ofp << nLay << '\n';
	for (i = 0; i < nLay; i++)
	{
		ofp << LayZ[i] <<' '<< LaySr[i] << '\n';
	}
	ofp.close();
	ofp.clear();

	sprintf(fileName, "%s/mlayersAdditional", path);
	ofp.open(fileName);
	ofp << scientific << setprecision(14);
	ofp << nLay << '\n';
	for (i=nLay-1;i>=0;i--)
	{
		ofp << 1.0 / LaySr[i] << ' ' << 1.0 / LaySr[i] << ' ' << 1.0 / LaySz[i] << ' ' << 1.0<< ' ' << 0.0 << '\n';
	}
	ofp.close();
	ofp.clear();

	inf.open("gen");
	if (!inf)
	{
		return 1;
	}
	inf >> z >> ax >> ay >> bx >> by;
	inf.close();
	inf.clear();

	sprintf(fileName, "%s/sours", path);
	ofp.open(fileName);
	ofp << scientific << setprecision(14);
	ofp << ax << ' ' << ay << ' ' << z << ' ' << bx << ' ' << by << ' ' << z << '\n';
	ofp.close();
	ofp.clear();

	sprintf(fileName, "%s/xyzVectorB", path); CopyFileA("recB", fileName, false);
	sprintf(fileName, "%s/xyzVectorE", path); CopyFileA("recE", fileName, false);

	sprintf(fileName, "%s/RegularMeshBuilderSettings.cfg", path);
	if (WriteMeshSettings(fileName, settings) != 0)
		return 1;

	sprintf(fileName, "%s/RegularMeshBuilderSettings2D.cfg", path);
	if (WriteMeshSettings2D(fileName, settings) != 0)
		return 1;

	return 0;
}
//------------------------------------------------------     
// Prepare calculations directory
//------------------------------------------------------     
int PrepareDirectoryCalculation(const char *path, const char *modulesPath, Settings &settings)
{
	char fileName[2048];
	int nRec;
	ifstream inf;
	ofstream ofp;

	sprintf(fileName, "%s/ifreq", path);
	if (WriteSingleStringFile(fileName, "1") != 0)
		return 1;

	sprintf(fileName, "%s/dminsrc", path);
	if (WriteSingleStringFile(fileName, "0.1") != 0)
		return 1;

	sprintf(fileName, "%s/n_neib", path);
	if (WriteSingleStringFile(fileName, "3") != 0)
		return 1;

	sprintf(fileName, "%s/scalingZ", path);
	if (WriteSingleStringFile(fileName, "1") != 0)
		return 1;

	sprintf(fileName, "%s/numberofabdipoles", path);
	if (!isFileExists(fileName))
		if (WriteSingleStringFile(fileName, std::to_string(settings.numberOfABDipoles).c_str()) != 0)
			return 1;

	sprintf(fileName, "%s/nu", path);
	if (WriteSingleStringFile(fileName, std::to_string(settings.Nu).c_str()) != 0)
		return 1;

	sprintf(fileName, "%s/met2d", path);
	if (WriteSingleStringFile(fileName, std::to_string(settings.met2d).c_str()) != 0)
		return 1;

	sprintf(fileName, "%s/met_slae_3d", path);
	ofp.open(fileName);
	ofp << settings.met_slae_3d << '\n';
	ofp << settings.eps3d << '\n';
	ofp << settings.maxiter3d << '\n';
	ofp.close();
	ofp.clear();

	inf.open("recB");
	if (!inf)
	{
		return 1;
	}
	inf >> nRec;
	inf.close();
	inf.clear();

	sprintf(fileName, "%s/recvsb", path);
	ofp.open(fileName);
	ofp << nRec << '\n';
	ofp.close();
	ofp.clear();

	inf.open("recE");
	if (!inf)
	{
		return 1;
	}
	inf >> nRec;
	inf.close();
	inf.clear();

	sprintf(fileName, "%s/recvse", path);
	ofp.open(fileName);
	ofp << nRec << '\n';
	ofp.close();
	ofp.clear();

	sprintf(fileName, "%s/group", path);
	ofp.open(fileName);
	ofp << 1 << '\n';
	ofp << 1 << ' ' << 1 << ' ' << 1 << '\n';
	ofp.close();
	ofp.clear();

	sprintf(fileName, "%s/nthreads.txt", path);
	if (WriteSingleIntFile(fileName, settings.threadsCount) != 0)
		return 1;

	return 0;
}
//------------------------------------------------------     
// Prepare calculation directories
//------------------------------------------------------     
int PrepareDirectory(Settings &settings, const char *calculationPath, const char *modulesPath, const char* ResultsPath)
{
	_mkdir(calculationPath);
	_mkdir(ResultsPath);

	PrepareDirectoryForMeshing(calculationPath, settings);
	PrepareDirectoryCalculation(calculationPath, modulesPath, settings);

	return 0;
}
//------------------------------------------------------     
// Run calculation
//------------------------------------------------------     
int RunCalculation(const char *calculationPath, const char *modulesPath)
{
	char fileName[2048];

	write_to_log("Starting CalcFreq.exe\n");
	sprintf(fileName, "%s/CalcFreq.exe", modulesPath);
	if (ExecuteExe(fileName, calculationPath) != 0) return 1;

	return 0;
}
//------------------------------------------------------     
// Calling common procedures
//------------------------------------------------------     
int MainProcedure()
{
	Settings settings;

	const char *calculationPath = "Calculations";
	const char *modulesPath = "Modules";

	write_to_log("Reading settings.cfg\n");
	if (ReadSettings("settings.cfg", settings) != 0)
		return 1;

	write_to_log("Preparing directories\n");
	if (PrepareDirectory(settings, calculationPath, modulesPath,"Results") != 0)
		return 1;

	write_to_log("Running calculation\n");
	if (RunCalculation(calculationPath, modulesPath) != 0)
		return 1;

	return 0;
}

int main()
{
	char buf[2048];

	if (open_log("CalcStarter.log") != 0)
	{
		printf("Cound not open CalcStarter.log\n");
		return 1;
	}

	int status = MainProcedure();
	sprintf(buf, "Status = %d\n", status);
	write_to_log(buf);

	fprintf(log_file, "All done\n");
	fclose(log_file);
	return status;
}
