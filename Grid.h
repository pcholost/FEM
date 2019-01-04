#pragma once
#include"Node.h"
#include"Element.h"
#include"File.h"
#include"LocalExecute.h"
#include <vector>
#include <iostream>
#include <memory>
#include <string>
#include <algorithm>
#include <iterator>

using namespace std;

class Grid
{
public:
	double **globalMatrixH;
	double **globalMatrixC;
	double *globalVectorP;
	double **globalMatrixHBC;
	double *vectorP;
	vector<long double> vectorOfTemperatures;
	vector<Node*> nodeTab;
	double **matrixHAndVectorP;

	vector<Element*> elementTab;

	File file;
	LocalExecute local;

	double a = file.getNodeNumberLength()*file.getNodeNumberLength();

	int modifiedIndexI, modifiedIndexJ;

	void setNodesCordinates();
	void setNodeTab();
	void setElementTab();
	void setElementID();
	void showElementID();
	void showNodesCordinates();
	void setNodesInElement();
	void setBCInNodesElement();
	void setValueOfGlobalMatrix(Element*);
	void setValueOfVectorP(Element*);
	void setValueofMatrixHAndVectorP();
	void calculateTemperatures();
	void upgradeVectorP();
	bool gaussMethod(int);
	void showGlobalVectorP();
	void showGlobalMatrixH();
	void showGlobalMatrixC();
	void setMatrixs();
	void modifyIndexes(int, int, Element*);
	Grid();
	~Grid();
};