#include<iostream>
#include<string>
#include"Element.h"
#include"File.h"

using namespace std;

Element::Element()
{
	id = new int[4];
	nodeTab = new Node[4];
	initTabs();
}

Element::Element(double conduct)
{
	this->conductivity = conduct;
	id = new int[4];
	initTabs();
	nodeTab = new Node[4];
}

Element::~Element()
{

}

void Element::showMatrixes()
{
	cout << "Macierz C: " << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cout << matrixC[i][j] << "  ";
		}
		cout << endl;
	}

	cout << "Macierz H: " << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cout << matrixH[i][j] << "  ";
		}
		cout << endl;
	}

	cout << "Macierz HBC: " << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cout << matrixHBC[i][j] << "  ";
		}
		cout << endl;
	}

	cout << "Wektor P: " << endl;
	for (int i = 0; i < 4; i++)
	{
		cout << vectorP[i] << " ";
	}
	cout << endl;
}

void Element::initTabs()
{
	matrixH = new double*[4];
	for (int i = 0; i < 4; i++)
		matrixH[i] = new double[4];

	matrixC = new double*[4];
	for (int i = 0; i < 4; i++)
		matrixC[i] = new double[4];

	vectorP = new double[4];
}
void Element::showId()
{
	for (int i = 0; i < 4; i++)
	{
		cout << id[i] << endl;
	}
}


void Element::showNodes()
{
	for (int i = 0; i < 4; i++)
	{
		cout << nodeTab[i].x << " " << nodeTab[i].y << " " << nodeTab[i].BC << endl;
	}
}

