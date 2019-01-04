#include<string>
#include<iostream>
#include"File.h"

using namespace std;

File::File()
{
	plik.open("plik.csv", ios::in);
	if (!plik.good())
	{
		cout << "Nie udalo sie otworzyc pliku" << endl;
	}
	else
	{
		string linia;
		getline(plik, linia, ';');
		getline(plik, linia, '\n');
		gridHight = atof(linia.c_str());

		getline(plik, linia, ';');
		getline(plik, linia, '\n');
		gridLength = atof(linia.c_str());

		getline(plik, linia, ';');
		getline(plik, linia, '\n');
		nodeNumberHight = atoi(linia.c_str());

		getline(plik, linia, ';');
		getline(plik, linia, '\n');
		nodeNumberLength = atoi(linia.c_str());

		getline(plik, linia, ';');
		getline(plik, linia, '\n');
		initialTemperature = atof(linia.c_str());

		getline(plik, linia, ';');
		getline(plik, linia, '\n');
		ambientTemperature = atof(linia.c_str());

		getline(plik, linia, ';');
		getline(plik, linia, '\n');
		conductivity = atof(linia.c_str());

		getline(plik, linia, ';');
		getline(plik, linia, '\n');
		integrationPoint = atof(linia.c_str());

		getline(plik, linia, ';');
		getline(plik, linia, '\n');
		pointCount = atoi(linia.c_str());

		getline(plik, linia, ';');
		getline(plik, linia, '\n');
		density = atof(linia.c_str());

		getline(plik, linia, ';');
		getline(plik, linia, '\n');
		specificHeat = atof(linia.c_str());

		getline(plik, linia, ';');
		getline(plik, linia, '\n');
		convection = atof(linia.c_str());

		getline(plik, linia, ';');
		getline(plik, linia, '\n');
		dTau = atof(linia.c_str());

		getline(plik, linia, ';');
		getline(plik, linia, '\n');
		simulationTime = atof(linia.c_str());

		tabX = new double[pointCount];
		tabY = new double[pointCount];

		for (int i = 0; i < pointCount; i++)
		{
			getline(plik, linia, ';');
			tabX[i] = atof(linia.c_str());
		}

		for (int i = 0; i < pointCount; i++)
		{
			tabY[i] = atof(linia.c_str());
			getline(plik, linia, ';');
		}

		plik.close();
	}

}


File::~File()
{

}

double File::getGridHight()
{
	return gridHight;
}

double File::getGridLength()
{
	return gridLength;
}

int File::getNodeNumberHight()
{
	return nodeNumberHight;
}

int File::getNodeNumberLength()
{
	return nodeNumberLength;
}

double File::getConductivity()
{
	return conductivity;
}

double File::getInitialTemperature()
{
	return initialTemperature;
}

double File::getAmbientTemperature()
{
	return ambientTemperature;
}

double File::getIntegrationPoint()
{
	return integrationPoint;
}

int File::getPointCount()
{
	return pointCount;
}

double File::getDensity()
{
	return density;
}

double File::getSpecificHeat()
{
	return specificHeat;
}

double File::getConvection()
{
	return convection;
}

double File::getdTau()
{
	return dTau;
}

double File::getSimulationTime()
{
	return simulationTime;
}

void File::showTabXY()
{
	for (int i = 0; i < pointCount; i++)
	{
		cout << "X" << i << ": " << tabX[i] << endl;
		cout << "Y" << i << ": " << tabY[i] << endl;
	}
}

double* File::getTabX()
{
	return this->tabX;
}

double* File::getTabY()
{
	return this->tabY;
}


double File::gridHight;
double File::gridLength;
int File::nodeNumberHight;
int File::nodeNumberLength;
double File::initialTemperature;
double File::ambientTemperature;
double File::conductivity;
double File::integrationPoint;
int File::pointCount;
double File::density;
double File::specificHeat;
double File::convection;
double File::dTau;
double File::simulationTime;