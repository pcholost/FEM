#pragma once
#include<fstream>
#include<vector>

class File
{
	std::fstream plik;
	static double gridHight;
	static double gridLength;
	static int nodeNumberHight;
	static int nodeNumberLength;
	static double initialTemperature;
	static double ambientTemperature;
	static double conductivity;
	static int pointCount;
	static double density;
	static double specificHeat;
	static double convection;
	static double dTau;
	static double simulationTime;
	double *tabX;
	double *tabY;
	static double integrationPoint;

public:
	File();
	~File();
	static double getGridHight();
	static double getGridLength();
	static int getNodeNumberHight();
	static int getNodeNumberLength();
	static double getInitialTemperature();
	static double getAmbientTemperature();
	static double getConductivity();
	static double getIntegrationPoint();
	static int getPointCount();
	static double getDensity();
	static double getSpecificHeat();
	static double getConvection();
	static double getdTau();
	static double getSimulationTime();
	double* getTabX();
	double* getTabY();

	void showTabXY();
};