#pragma once
#include"Node.h"
#include"File.h"

class Element
{
public:
	File plik;
	int *id;
	double conductivity;
	double *vectorP;
	double **matrixC;
	double **matrixH;
	double **matrixHBC;
	Node *nodeTab;
	Node *tabNode;

	Element();
	Element(double);
	~Element();

	void showId();
	void showNodes();
	void showMatrixes();
	void initTabs();

};