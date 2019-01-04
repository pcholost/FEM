#include"Grid.h"
#include"File.h"
#include<iomanip>
#include<iostream>
#include<cmath>


Grid::Grid()
{
	int nodeSize = static_cast<int>(file.getNodeNumberHight()) * static_cast<int>(file.getNodeNumberLength());
	vectorOfTemperatures.assign(file.getNodeNumberLength() * file.getNodeNumberHight(), file.getInitialTemperature());
	setNodeTab();
	setElementTab();

	globalMatrixH = new double*[nodeSize];
	for (int i = 0; i < nodeSize; i++)
		globalMatrixH[i] = new double[nodeSize];

	globalMatrixC = new double*[nodeSize];
	for (int i = 0; i < nodeSize; i++)
		globalMatrixC[i] = new double[nodeSize];

	globalVectorP = new double[nodeSize];
	vectorP = new double[nodeSize];

	matrixHAndVectorP = new double*[nodeSize];
	for (int i = 0; i < nodeSize; i++)
		matrixHAndVectorP[i] = new double[nodeSize];

	globalMatrixHBC = new double*[nodeSize];
	for (int i = 0; i < nodeSize; i++)
		globalMatrixHBC[i] = new double[nodeSize];

	for (int i = 0; i < file.getNodeNumberHight() * file.getNodeNumberHight(); i++)
	{
		for (int j = 0; j < file.getNodeNumberHight() * file.getNodeNumberHight(); j++)
		{
			globalMatrixC[i][j] = 0;
			globalMatrixH[i][j] = 0;
			globalMatrixHBC[i][j] = 0;
		}
		globalVectorP[i] = 0;
	}

	setNodesCordinates();
	setElementID();
	setNodesInElement();
	setBCInNodesElement();
	setMatrixs();

}

Grid::~Grid()
{

}

void Grid::setNodeTab()
{
	for (int i = 0; i < file.getNodeNumberLength() * file.getNodeNumberHight(); i++)
		nodeTab.emplace_back(new Node());
}

void Grid::setElementTab()
{
	for (int i = 0; i < (file.getNodeNumberLength() - 1) * (file.getNodeNumberHight() - 1); i++)
		elementTab.emplace_back(new Element(file.getConductivity()));
}

void Grid::setNodesCordinates()
{
	double dh = file.getGridHight() / (file.getNodeNumberHight() - 1);
	double dl = file.getGridLength() / (file.getNodeNumberLength() - 1);
	int nhVal = static_cast<int>(file.getNodeNumberHight());

	nodeTab[0]->x = 0;
	nodeTab[0]->y = 0;
	int nodeVal = 0;
	for (int i = 1; i < (file.getNodeNumberHight() *file.getNodeNumberLength()); i++)
	{
		if (i - nhVal == nodeVal)
		{
			nodeVal = i;
			nodeTab[i]->x = nodeTab[i - 1]->x + dl;
			nodeTab[i]->y = 0;
		}

		else {
			nodeTab[i]->x = nodeTab[i - 1]->x;
			nodeTab[i]->y = nodeTab[i - 1]->y + dh;
		}
	}
}

void Grid::setElementID()
{
	int counter = 0;
	int loop = (file.getNodeNumberHight() - 1) * (file.getNodeNumberLength() - 1);
	int nhVal = static_cast<int>(file.getNodeNumberHight());
	int nlVal = static_cast<int>(file.getNodeNumberLength());
	for (int i = 1; i <= file.getNodeNumberHight() * file.getNodeNumberLength(); i++)
	{
		if (i%nhVal == 0)continue;

		if (i > nhVal * nlVal - nhVal)break;

		elementTab[counter]->id[0] = i;
		elementTab[counter]->id[1] = i + file.getNodeNumberHight();
		elementTab[counter]->id[2] = elementTab[counter]->id[1] + 1;
		elementTab[counter]->id[3] = elementTab[counter]->id[0] + 1;

		counter++;
	}
}

void Grid::showElementID()
{
	int i = 0;

	for (vector<Element*>::iterator iter = elementTab.begin(); iter != elementTab.end(); iter++)
	{
		cout << "Element nr: " << ++i << endl;
		(*iter)->showId();
	}
}

void Grid::showNodesCordinates()
{
	int i = 0;
	for (vector<Node*>::iterator iter = nodeTab.begin(); iter != nodeTab.end(); iter++)
	{
		cout << "Wezel nr: " << i << " " << "Wsp. x: " << (*iter)->x << " " << "Wsp. y: " << (*iter)->y << endl;
		i++;
	}
}

void Grid::setNodesInElement()
{
	int counter = 0;
	int nhVal = static_cast<int>(file.getNodeNumberHight());
	int nlVal = static_cast<int>(file.getNodeNumberLength());

	for (int i = 1; i <= nodeTab.size(); i++)
	{
		if (i % nhVal == 0) continue;

		if (i > nhVal * nlVal - nhVal)break;

		elementTab[counter]->nodeTab[0].x = nodeTab[i - 1]->x;
		elementTab[counter]->nodeTab[0].y = nodeTab[i - 1]->y;

		elementTab[counter]->nodeTab[1].x = nodeTab[(i + file.getNodeNumberHight()) - 1]->x;
		elementTab[counter]->nodeTab[1].y = nodeTab[(i + file.getNodeNumberHight()) - 1]->y;

		elementTab[counter]->nodeTab[2].x = nodeTab[(i + file.getNodeNumberHight())]->x;
		elementTab[counter]->nodeTab[2].y = nodeTab[(i + file.getNodeNumberHight())]->y;

		elementTab[counter]->nodeTab[3].x = nodeTab[i]->x;
		elementTab[counter]->nodeTab[3].y = nodeTab[i]->y;

		counter++;
	}
}

void Grid::setBCInNodesElement()
{
	int counter = 0;
	for (int i = 1; i <= file.getNodeNumberHight() - 1; i++)
	{
		for (int j = 1; j <= file.getNodeNumberHight() - 1; j++)
		{
			if (i == 1)
			{
				if (j == 1)
				{
					elementTab[counter]->nodeTab[0].BC = true;
					elementTab[counter]->nodeTab[1].BC = true;
					elementTab[counter]->nodeTab[2].BC = false;
					elementTab[counter]->nodeTab[3].BC = true;
				}
				else if (j > 1 && j < file.getNodeNumberHight() - 1)
				{
					elementTab[counter]->nodeTab[0].BC = true;
					elementTab[counter]->nodeTab[1].BC = false;
					elementTab[counter]->nodeTab[2].BC = false;
					elementTab[counter]->nodeTab[3].BC = true;
				}
				else
				{
					elementTab[counter]->nodeTab[0].BC = true;
					elementTab[counter]->nodeTab[1].BC = false;
					elementTab[counter]->nodeTab[2].BC = true;
					elementTab[counter]->nodeTab[3].BC = true;
				}
			}
			else if (i > 1 && i < file.getNodeNumberHight() - 1)
			{
				if (j == 1)
				{
					elementTab[counter]->nodeTab[0].BC = true;
					elementTab[counter]->nodeTab[1].BC = true;
					elementTab[counter]->nodeTab[2].BC = false;
					elementTab[counter]->nodeTab[3].BC = false;
				}
				else if (j > 1 && j < file.getNodeNumberHight() - 1)
				{
					elementTab[counter]->nodeTab[0].BC = false;
					elementTab[counter]->nodeTab[1].BC = false;
					elementTab[counter]->nodeTab[2].BC = false;
					elementTab[counter]->nodeTab[3].BC = false;
				}
				else
				{
					elementTab[counter]->nodeTab[0].BC = false;
					elementTab[counter]->nodeTab[1].BC = false;
					elementTab[counter]->nodeTab[2].BC = true;
					elementTab[counter]->nodeTab[3].BC = true;
				}
			}
			else
			{
				if (j == 1)
				{
					elementTab[counter]->nodeTab[0].BC = true;
					elementTab[counter]->nodeTab[1].BC = true;
					elementTab[counter]->nodeTab[2].BC = true;
					elementTab[counter]->nodeTab[3].BC = false;
				}
				else if (j > 1 && j < file.getNodeNumberHight() - 1)
				{
					elementTab[counter]->nodeTab[0].BC = false;
					elementTab[counter]->nodeTab[1].BC = true;
					elementTab[counter]->nodeTab[2].BC = true;
					elementTab[counter]->nodeTab[3].BC = false;
				}
				else
				{
					elementTab[counter]->nodeTab[0].BC = false;
					elementTab[counter]->nodeTab[1].BC = true;
					elementTab[counter]->nodeTab[2].BC = true;
					elementTab[counter]->nodeTab[3].BC = true;
				}
			}
			counter++;
		}
	}
}

void Grid::setValueofMatrixHAndVectorP()
{
	int counter = 0;
	for (int i = 0; i < file.getNodeNumberHight() * file.getNodeNumberHight(); i++)
	{
		for (int j = 0; j < file.getNodeNumberHight() * file.getNodeNumberHight() + 1; j++)
		{
			if (j == file.getNodeNumberHight() * file.getNodeNumberHight())
			{
				matrixHAndVectorP[i][j] = globalVectorP[i];
				continue;
			}
			matrixHAndVectorP[i][j] = globalMatrixHBC[i][j];
		}
	}
}

//zmienianie indeksow na etapie agregacji
void Grid::modifyIndexes(int i, int j, Element *element)
{
	switch (i + 1)
	{
	case 1:
		modifiedIndexI = element->id[0] - 1;
		break;
	case 2:
		modifiedIndexI = element->id[1] - 1;
		break;
	case 3:
		modifiedIndexI = element->id[2] - 1;
		break;
	case 4:
		modifiedIndexI = element->id[3] - 1;
		break;
	default:
		break;
	}

	switch (j + 1)
	{
	case 1:
		modifiedIndexJ = element->id[0] - 1;
		break;
	case 2:
		modifiedIndexJ = element->id[1] - 1;
		break;
	case 3:
		modifiedIndexJ = element->id[2] - 1;
		break;
	case 4:
		modifiedIndexJ = element->id[3] - 1;
		break;
	default:
		break;
	}
}

//agregacja macierzy H i C
void Grid::setValueOfGlobalMatrix(Element *element)
{
	double val, val1, val2;

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			val = element->matrixHBC[i][j];
			val1 = element->matrixC[i][j];
			val2 = element->matrixH[i][j];
			modifyIndexes(i, j, element);
			//cout << "miejsce w macierzy globalnej" << modifiedIndexI << " " << modifiedIndexJ << endl;
			globalMatrixC[modifiedIndexI][modifiedIndexJ] += val1;
			globalMatrixHBC[modifiedIndexI][modifiedIndexJ] += val;
			globalMatrixH[modifiedIndexI][modifiedIndexJ] += val2;
		}
	}
}

void Grid::upgradeVectorP()
{
	double val1 = 0;

	for (int j = 0; j < file.getNodeNumberLength() * file.getNodeNumberHight(); j++)
	{
		globalVectorP[j] = 0;
		for (int i = 0; i < file.getNodeNumberLength() * file.getNodeNumberHight(); i++)
		{
			val1 += ((globalMatrixC[j][i] / file.getdTau()) * vectorOfTemperatures[i]);
		}
		globalVectorP[j] = vectorP[j] + val1;
		val1 = 0;
	}
}

void Grid::setMatrixs()
{
	vector<LocalExecute*> matrixVec;
	double val1 = 0;
	int i = 0;
	for (vector<Element*>::iterator iter = elementTab.begin(); iter != elementTab.end(); iter++)
	{
		/*cout << "Element nr " << ++i << endl;
		(*iter)->showId();*/
		matrixVec.emplace_back(new LocalExecute(*(iter)));

		setValueOfGlobalMatrix(*(iter));
		setValueOfVectorP(*(iter));
	}
	setValueofMatrixHAndVectorP();
	//matrixHAndVectorP = move(globalMatrixH);
	for (int j = 0; j < file.getNodeNumberLength() * file.getNodeNumberHight(); ++j)
	{
		for (int i = 0; i < file.getNodeNumberLength() * file.getNodeNumberHight(); ++i)
		{
			globalMatrixHBC[j][i] += globalMatrixC[j][i] / file.getdTau();
		}
	}
	/*cout << "globalna macierz P bez dodawania\n";
	showGlobalVectorP();
	cout << endl;*/
	for (int k = 0; k < file.getNodeNumberLength() * file.getNodeNumberHight(); k++)
	{
		vectorP[k] = globalVectorP[k];
	}
	for (int j = 0; j < file.getNodeNumberLength() * file.getNodeNumberHight(); j++)
	{
		for (int i = 0; i < file.getNodeNumberLength() * file.getNodeNumberHight(); i++)
		{
			val1 += (globalMatrixC[j][i] / file.getdTau() * vectorOfTemperatures[i]);
		}
		/* cout << val1 << "  ";
		cout << endl;*/
		globalVectorP[j] = globalVectorP[j] + val1;
		val1 = 0;
	}
}

void Grid::calculateTemperatures()
{
	int n = static_cast<int>(file.getNodeNumberHight() * file.getNodeNumberHight());
	long double max = 0;
	long double min = 0;
	int simulationTime = static_cast<int>(file.getSimulationTime());
	int dTau = static_cast<int>(file.getdTau());
	for (int j = 0; j < simulationTime; j = j + dTau)
	{
		setValueofMatrixHAndVectorP();
		//cout << "wektor p\n";
		//showGlobalVectorP();
		gaussMethod(n);
		cout << "Temperatury: " << endl;
		std::copy(vectorOfTemperatures.begin(), vectorOfTemperatures.end(),
			std::ostream_iterator<long double>(std::cout, " "));
		cout << endl;
		upgradeVectorP();
		//setValueOfVectorP();
	}
}

bool Grid::gaussMethod(int n)
{
	const double eps = 1e-12;
	int i, j, k;
	double m, s;
	fill(vectorOfTemperatures.begin(), vectorOfTemperatures.end(), 0);
	//eliminacja wspolczynnikow

	for (i = 0; i < n - 1; i++)
	{
		for (j = i + 1; j < n; j++)
		{
			if (fabs(matrixHAndVectorP[i][i]) < eps)
				return false;
			m = -matrixHAndVectorP[j][i] / matrixHAndVectorP[i][i];
			for (k = i + 1; k <= n; k++)
				matrixHAndVectorP[j][k] += m * matrixHAndVectorP[i][k];
		}
	}

	//wyliczanie niewiadomych

	for (i = n - 1; i >= 0; i--)
	{
		s = matrixHAndVectorP[i][n];
		for (j = n - 1; j >= i + 1; j--)
			s -= matrixHAndVectorP[i][j] * vectorOfTemperatures[j];
		if (fabs(matrixHAndVectorP[i][i]) < eps)
			return false;
		vectorOfTemperatures[i] = s / matrixHAndVectorP[i][i];
	}
	return true;

}

void Grid::setValueOfVectorP(Element *element)
{
	double val, val1 = 0;
	for (int i = 0; i < 4; i++)
	{
		val = element->vectorP[i];
		modifyIndexes(i, 0, element);
		globalVectorP[modifiedIndexI] += val;
	}
}

void Grid::showGlobalVectorP()
{
	cout << "Wektor P globalny: " << endl;
	for (int i = 0; i < file.getNodeNumberLength() * file.getNodeNumberHight(); i++)
	{
		cout << globalVectorP[i] << endl;
	}
}

void Grid::showGlobalMatrixH()
{
	cout << "Macierz H globalna: " << endl;
	for (int i = 0; i < file.getNodeNumberLength() * file.getNodeNumberHight(); i++)
	{
		for (int j = 0; j < file.getNodeNumberLength() * file.getNodeNumberHight(); j++)
		{
			cout << globalMatrixH[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
	cout << endl;
	cout << "Macierz H i P: " << endl;
	for (int i = 0; i < file.getNodeNumberLength() * file.getNodeNumberHight(); i++)
	{
		for (int j = 0; j < file.getNodeNumberLength() * file.getNodeNumberHight() - 1; j++)
		{
			cout << matrixHAndVectorP[i][j] << " ";
		}
		cout << endl;
	}
}

void Grid::showGlobalMatrixC()
{
	cout << "Macierz C globalna: " << endl;
	for (int i = 0; i < file.getNodeNumberLength() * file.getNodeNumberHight(); i++)
	{
		for (int j = 0; j < file.getNodeNumberLength() * file.getNodeNumberHight(); j++)
		{
			cout << globalMatrixC[i][j] << " ";
		}
		cout << endl;
	}
}

