#include"LocalExecute.h"
#include"File.h"
#include<iostream>

using namespace std;

LocalExecute::LocalExecute()
{

}

LocalExecute::~LocalExecute()
{

}

LocalExecute::LocalExecute(Element* element)
{
	//w konstuktorze wykonujemy wszystkie potrzebne funkcje dla kazdego z elementow
	this->element = element;
	funkcjeKsztaltu();
	matrixJ();
	matrixH();
	matrixC();
	vectorP();
	//przekazanie tablic do klasy element, z ktorej korzysta grid
	element->matrixC = move(tabMacierzC);
	element->matrixHBC = move(tabMacierzHBC);
	element->vectorP = move(tabWektorP);
	element->matrixH = move(tabMacierzH);
}

//integrationPoint wczytany jest z pliku i wykorzystywany do wyliczenia Ksi i Eta

void LocalExecute::setTabKsi(int size, double integrationPoint)
{
	this->tabKsi = new double[size];
	this->tabKsi[0] = -integrationPoint;
	this->tabKsi[1] = integrationPoint;
	this->tabKsi[2] = integrationPoint;
	this->tabKsi[3] = -integrationPoint;
}

void LocalExecute::setTabEta(int size, double integrationPoint)
{
	this->tabEta = new double[size];
	this->tabEta[0] = -integrationPoint;
	this->tabEta[1] = -integrationPoint;
	this->tabEta[2] = integrationPoint;
	this->tabEta[3] = integrationPoint;
}

double* LocalExecute::getTabEta()
{
	return tabEta;
}

double* LocalExecute::getTabKsi()
{
	return tabKsi;
}

//funkcje pomocnicze do liczenia f.ksztaltu

double LocalExecute::N1(double tabKsi, double tabEta)
{
	return 0.25*(1 - tabKsi)*(1 - tabEta);
}
double LocalExecute::N2(double tabKsi, double tabEta)
{
	return 0.25*(1 + tabKsi)*(1 - tabEta);
}
double LocalExecute::N3(double tabKsi, double tabEta)
{
	return 0.25*(1 + tabKsi)*(1 + tabEta);
}
double LocalExecute::N4(double tabKsi, double tabEta)
{
	return 0.25*(1 - tabKsi)*(1 + tabEta);
}

double** LocalExecute::N(int size) //funkcje ksztaltu
{

	tabN = new double*[size];
	for (int i = 0; i < size; i++)
		tabN[i] = new double[size];
	for (int i = 0; i < size; i++)
	{
		tabN[i][0] = N1(this->tabKsi[i], this->tabEta[i]);
		tabN[i][1] = N2(this->tabKsi[i], this->tabEta[i]);
		tabN[i][2] = N3(this->tabKsi[i], this->tabEta[i]);
		tabN[i][3] = N4(this->tabKsi[i], this->tabEta[i]);
	}

	return tabN;
}

//funckje pomocniczne do liczenia pochodnych funkcji ksztaltu

double LocalExecute::FKszt1Eta(double tabKsi)
{
	return -0.25*(1 - tabKsi);
}

double LocalExecute::FKszt2Eta(double tabKsi)
{
	return -0.25*(1 + tabKsi);
}

double LocalExecute::FKszt3Eta(double tabKsi)
{
	return 0.25*(1 + tabKsi);
}

double LocalExecute::FKszt4Eta(double tabKsi)
{
	return 0.25*(1 - tabKsi);
}

double LocalExecute::FKszt1Ksi(double tabEta)
{
	return -0.25*(1 - tabEta);
}

double LocalExecute::FKszt2Ksi(double tabEta)
{
	return 0.25*(1 - tabEta);
}

double LocalExecute::FKszt3Ksi(double tabEta)
{
	return 0.25*(1 + tabEta);
}

double LocalExecute::FKszt4Ksi(double tabEta)
{
	return -0.25*(1 + tabEta);
}

double** LocalExecute::pochodnaFunkcjiKsztaltuEta(int size)
{
	tabFunkcjiKsztaltuEta = new double*[size];
	for (int i = 0; i < size; i++)
	{
		tabFunkcjiKsztaltuEta[i] = new double[size];
	}
	for (int i = 0; i<size; i++)
	{
		tabFunkcjiKsztaltuEta[i][0] = FKszt1Eta(tabKsi[i]);
		tabFunkcjiKsztaltuEta[i][1] = FKszt2Eta(tabKsi[i]);
		tabFunkcjiKsztaltuEta[i][2] = FKszt3Eta(tabKsi[i]);
		tabFunkcjiKsztaltuEta[i][3] = FKszt4Eta(tabKsi[i]);
	}
	return tabFunkcjiKsztaltuEta;
}

double** LocalExecute::pochodnaFunkcjiKsztaltuKsi(int size)
{
	tabFunkcjiKsztaltuKsi = new double*[size];
	for (int i = 0; i < size; i++)
	{
		tabFunkcjiKsztaltuKsi[i] = new double[size];
	}
	for (int i = 0; i<size; i++)
	{
		tabFunkcjiKsztaltuKsi[i][0] = FKszt1Ksi(tabEta[i]);
		tabFunkcjiKsztaltuKsi[i][1] = FKszt2Ksi(tabEta[i]);
		tabFunkcjiKsztaltuKsi[i][2] = FKszt3Ksi(tabEta[i]);
		tabFunkcjiKsztaltuKsi[i][3] = FKszt4Ksi(tabEta[i]);
	}
	return tabFunkcjiKsztaltuKsi;
}

//funkcje pomocniczne do liczenia jacobianu przeksztalcenia

double LocalExecute::calculateJacobianX(double *tab)
{
	double suma = 0;
	for (int i = 0; i < File::getPointCount(); i++)
	{
		suma += tab[i] * (element->nodeTab[i].x);
	}
	return suma;
}

double LocalExecute::calculateJacobianY(double *tab)
{
	double suma = 0;
	for (int i = 0; i < File::getPointCount(); i++)
	{
		suma += tab[i] * (element->nodeTab[i].y);
	}
	return suma;
}

//mnozymy krawedzie elementow, x i y, razy pochodne funkcji ksztaltu otrzymane przed chwila
//dzieki jacobianowi przeksztalcenia mozemy przez pomiedzy ukladem lokalnym a globalnym
double** LocalExecute::obliczJacobian(int size)
{
	jacobian = new double*[size];
	for (int i = 0; i < size; i++)
	{
		jacobian[i] = new double[size];
	}

	for (int i = 0; i < size; i++)
	{
		jacobian[0][i] = calculateJacobianX(tabFunkcjiKsztaltuKsi[0]);
		jacobian[1][i] = calculateJacobianY(tabFunkcjiKsztaltuKsi[1]);
		jacobian[2][i] = calculateJacobianX(tabFunkcjiKsztaltuEta[2]);
		jacobian[3][i] = calculateJacobianY(tabFunkcjiKsztaltuEta[3]);
	}

	return jacobian;
}

//detJ mowi nam jaki blad wystepuje pomiedzy przejsciem miedzy ukladami

double LocalExecute::detJ(double *a)
{
	return a[0] * a[3] - a[1] * a[2];
}

double* LocalExecute::detJacobian(int size)
{
	tabDetJ = new double[size];
	tabDetJ[0] = (jacobian[0][0] * jacobian[3][0] - (jacobian[2][0] * jacobian[1][0]));
	tabDetJ[1] = (jacobian[0][1] * jacobian[3][1] - (jacobian[2][1] * jacobian[1][1]));
	tabDetJ[2] = (jacobian[0][2] * jacobian[3][2] - (jacobian[2][2] * jacobian[1][2]));
	tabDetJ[3] = (jacobian[0][3] * jacobian[3][3] - (jacobian[2][3] * jacobian[1][3]));

	return tabDetJ;
}

//dzielimy jacobian przez detj i otrzymujemy macierz Jacobiego

double** LocalExecute::obliczJacobianPrzezDet(int size)
{
	jacobianPrzezDet = new double*[size];
	for (int i = 0; i < size; i++)
	{
		jacobianPrzezDet[i] = new double[size];

	}
	for (int i = 0; i < size; i++)
	{

		jacobianPrzezDet[i][0] = jacobian[3][i] / tabDetJ[i];
		jacobianPrzezDet[i][1] = -jacobian[1][i] / tabDetJ[i];
		jacobianPrzezDet[i][2] = jacobian[2][i] / tabDetJ[i];
		jacobianPrzezDet[i][3] = jacobian[0][i] / tabDetJ[i];
	}

	/*cout << "Macierz Jacobianu: " << endl;
	for (int i = 0; i < size; i++)
	{
	for (int j = 0; j < size; j++)
	{
	cout << jacobianPrzezDet[i][j] << "  ";
	}
	cout << endl;
	}*/


	return jacobianPrzezDet;
}

double** LocalExecute::pochodnaJacobianuPoDx(int size)
{
	jacobianDx = new double*[size];
	for (int i = 0; i < size; i++)
	{
		jacobianDx[i] = new double[size];
	}
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			jacobianDx[i][j] = jacobianPrzezDet[i][0] * tabFunkcjiKsztaltuKsi[i][j] + jacobianPrzezDet[i][1] * tabFunkcjiKsztaltuEta[i][j];
			jacobianDx[i][j] = jacobianPrzezDet[i][0] * tabFunkcjiKsztaltuKsi[i][j] + jacobianPrzezDet[i][1] * tabFunkcjiKsztaltuEta[i][j];
			jacobianDx[i][j] = jacobianPrzezDet[i][0] * tabFunkcjiKsztaltuKsi[i][j] + jacobianPrzezDet[i][1] * tabFunkcjiKsztaltuEta[i][j];
			jacobianDx[i][j] = jacobianPrzezDet[i][0] * tabFunkcjiKsztaltuKsi[i][j] + jacobianPrzezDet[i][1] * tabFunkcjiKsztaltuEta[i][j];
		}
	}

	return jacobianDx;
}


double** LocalExecute::pochodnaJacobianuPoDy(int size)
{
	jacobianDy = new double*[size];
	for (int i = 0; i < size; i++)
	{
		jacobianDy[i] = new double[size];
	}
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			jacobianDy[i][j] = jacobianPrzezDet[i][2] * tabFunkcjiKsztaltuKsi[i][j] + jacobianPrzezDet[i][3] * tabFunkcjiKsztaltuEta[i][j];
			jacobianDy[i][j] = jacobianPrzezDet[i][2] * tabFunkcjiKsztaltuKsi[i][j] + jacobianPrzezDet[i][3] * tabFunkcjiKsztaltuEta[i][j];
			jacobianDy[i][j] = jacobianPrzezDet[i][2] * tabFunkcjiKsztaltuKsi[i][j] + jacobianPrzezDet[i][3] * tabFunkcjiKsztaltuEta[i][j];
			jacobianDy[i][j] = jacobianPrzezDet[i][2] * tabFunkcjiKsztaltuKsi[i][j] + jacobianPrzezDet[i][3] * tabFunkcjiKsztaltuEta[i][j];
		}
	}
	return jacobianDy;
}

//macierz transponujemy {dN/dx}{dn/dx) i {dn/dy}{dn/dy} po 4 punktach calkowania kazdy

double*** LocalExecute::mnozenieMacierzy(int size, double** tab)
{
	macierzTransponowana = new double**[size];
	for (int i = 0; i < size; i++)
	{
		macierzTransponowana[i] = new double*[size];
	}
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			macierzTransponowana[i][j] = new double[size];
		}
	}

	for (int k = 0; k < size; k++)
	{
		for (int i = 0; i < size; i++)
		{
			macierzTransponowana[k][i][0] = tab[k][i] * tab[k][0];
			macierzTransponowana[k][i][1] = tab[k][i] * tab[k][1];
			macierzTransponowana[k][i][2] = tab[k][i] * tab[k][2];
			macierzTransponowana[k][i][3] = tab[k][i] * tab[k][3];
		}
	}
	return macierzTransponowana;
}

//mnozymy kazda z macierzy razy swoj wyznacznik Jacobianu

double*** LocalExecute::macierzRazyDet(int size, double*** tab)
{
	macierzDet = new double**[size];
	for (int i = 0; i < size; i++)
	{
		macierzDet[i] = new double*[size];
	}
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			macierzDet[i][j] = new double[size];
		}
	}

	for (int k = 0; k < size; k++)
	{
		for (int i = 0; i < size; i++)
		{
			macierzDet[k][i][0] = tab[k][i][0] * tabDetJ[i];
			macierzDet[k][i][1] = tab[k][i][1] * tabDetJ[i];
			macierzDet[k][i][2] = tab[k][i][2] * tabDetJ[i];
			macierzDet[k][i][3] = tab[k][i][3] * tabDetJ[i];
		}
	}

	return macierzDet;
}

//mnozymy wartosc razy alfe kondukcji

double*** LocalExecute::macierzIkondukcja(int size, double conduction, double*** tab, double*** tab2)
{
	macierzConduction = new double**[size];
	for (int i = 0; i < size; i++)
	{
		macierzConduction[i] = new double*[size];
	}
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			macierzConduction[i][j] = new double[size];
		}
	}
	for (int k = 0; k < size; k++)
	{
		for (int i = 0; i < size; i++)
		{
			macierzConduction[k][i][0] = conduction * (tab[k][i][0] + tab2[k][i][0]);
			macierzConduction[k][i][1] = conduction *  (tab[k][i][1] + tab2[k][i][1]);
			macierzConduction[k][i][2] = conduction * (tab[k][i][2] + tab2[k][i][2]);
			macierzConduction[k][i][3] = conduction *  (tab[k][i][3] + tab2[k][i][3]);
		}
	}
	return macierzConduction;
}

//otrzymujemy macierz H jeszcze bez warunkow brzegowych

double** LocalExecute::macierzH(int size)
{
	tabMacierzH = new double*[size];
	for (int i = 0; i < size; i++)
	{
		tabMacierzH[i] = new double[size];
	}
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			tabMacierzH[i][j] = macierzConduction[0][i][j] + macierzConduction[1][i][j] + macierzConduction[2][i][j] + macierzConduction[3][i][j];
		}
	}
	return tabMacierzH;
}

//majac wszystkie dane liczymy calke po objetosci mnozac funkcje ksztaltu przez siebie, wyznacznik jacobianu, gestosc i cieplo wlasciwe

double*** LocalExecute::macierzCalkowania(int size, double ro, double c)
{
	tabMacierzyCalkowania = new double**[size];
	for (int i = 0; i < size; i++)
	{
		tabMacierzyCalkowania[i] = new double*[size];
	}
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			tabMacierzyCalkowania[i][j] = new double[size];
		}
	}
	for (int k = 0; k < size; k++)
	{
		for (int i = 0; i < size; i++)
		{
			tabMacierzyCalkowania[k][i][0] = tabN[k][0] * tabN[k][i] * tabDetJ[k] * ro*c;
			tabMacierzyCalkowania[k][i][1] = tabN[k][1] * tabN[k][i] * tabDetJ[k] * ro*c;
			tabMacierzyCalkowania[k][i][2] = tabN[k][2] * tabN[k][i] * tabDetJ[k] * ro*c;
			tabMacierzyCalkowania[k][i][3] = tabN[k][3] * tabN[k][i] * tabDetJ[k] * ro*c;
		}
	}
	return tabMacierzyCalkowania;
}

double** LocalExecute::macierzC(int size)
{
	tabMacierzC = new double*[size];
	for (int i = 0; i < size; i++)
	{
		tabMacierzC[i] = new double[size];
	}
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			tabMacierzC[i][j] = tabMacierzyCalkowania[0][i][j] + tabMacierzyCalkowania[1][i][j] + tabMacierzyCalkowania[2][i][j] + tabMacierzyCalkowania[3][i][j];
		}
	}

	cout << "Macierz C:" << endl;
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			cout << tabMacierzC[i][j] << "  ";
		}
		cout << endl;
	}
	return tabMacierzC;
}

//pobieramy dlugosci bokow na podstawie danych do siatki (rozmiaru) wczytanego z pliku

double* LocalExecute::dlugoscBokow(int size)
{
	tabDlugosciBokow = new double[size];
	for (int i = 0; i < 4; i++)
	{
		tabDlugosciBokow[i] = sqrt(pow(element->nodeTab[i].x - element->nodeTab[(i + 1) % 4].x, 2) +
			pow(element->nodeTab[i].y - element->nodeTab[(i + 1) % 4].y, 2));
	}

	return tabDlugosciBokow;
}

//jacobian 1d

double* LocalExecute::detJBoku(int size)
{
	tabDetJBoku = new double[size];
	for (int i = 0; i < size; i++)
	{
		tabDetJBoku[i] = tabDlugosciBokow[i] / 2;
	}
	return tabDetJBoku;
}

//FUNKCJE KSZTALTU DLA ELEMENTU UNIWERSALNEGO

double*** LocalExecute::funkcjeKsztaltuBC(int size)
{
	tabFKBC = new double**[size];
	for (int i = 0; i < size; i++)
	{
		tabFKBC[i] = new double*[2];
	}
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			tabFKBC[i][j] = new double[size];
		}
	}

	//pow1
	{
		double eta = -1;
		double ksi = -(1 / sqrt(3));
		tabFKBC[0][0][0] = N1(ksi, eta);
		tabFKBC[0][0][1] = N2(ksi, eta);
		tabFKBC[0][0][2] = N3(ksi, eta);
		tabFKBC[0][0][3] = N4(ksi, eta);
	}
	{
		double eta = -1;
		double ksi = (1 / sqrt(3));
		tabFKBC[0][1][0] = N1(ksi, eta);
		tabFKBC[0][1][1] = N2(ksi, eta);
		tabFKBC[0][1][2] = N3(ksi, eta);
		tabFKBC[0][1][3] = N4(ksi, eta);
	}
	//pow2
	{
		double eta = -(1 / sqrt(3));
		double ksi = 1;
		tabFKBC[1][0][0] = N1(ksi, eta);
		tabFKBC[1][0][1] = N2(ksi, eta);
		tabFKBC[1][0][2] = N3(ksi, eta);
		tabFKBC[1][0][3] = N4(ksi, eta);
	}
	{
		double eta = (1 / sqrt(3));
		double ksi = 1;
		tabFKBC[1][1][0] = N1(ksi, eta);
		tabFKBC[1][1][1] = N2(ksi, eta);
		tabFKBC[1][1][2] = N3(ksi, eta);
		tabFKBC[1][1][3] = N4(ksi, eta);
	}
	//pow3
	{
		double eta = 1;
		double ksi = (1 / sqrt(3));
		tabFKBC[2][0][0] = N1(ksi, eta);
		tabFKBC[2][0][1] = N2(ksi, eta);
		tabFKBC[2][0][2] = N3(ksi, eta);
		tabFKBC[2][0][3] = N4(ksi, eta);
	}
	{
		double eta = 1;
		double ksi = -(1 / sqrt(3));
		tabFKBC[2][1][0] = N1(ksi, eta);
		tabFKBC[2][1][1] = N2(ksi, eta);
		tabFKBC[2][1][2] = N3(ksi, eta);
		tabFKBC[2][1][3] = N4(ksi, eta);
	}
	//pow4
	{
		double eta = (1 / sqrt(3));
		double ksi = -1;
		tabFKBC[3][0][0] = N1(ksi, eta);
		tabFKBC[3][0][1] = N2(ksi, eta);
		tabFKBC[3][0][2] = N3(ksi, eta);
		tabFKBC[3][0][3] = N4(ksi, eta);
	}
	{
		double eta = -(1 / sqrt(3));
		double ksi = -1;
		tabFKBC[3][1][0] = N1(ksi, eta);
		tabFKBC[3][1][1] = N2(ksi, eta);
		tabFKBC[3][1][2] = N3(ksi, eta);
		tabFKBC[3][1][3] = N4(ksi, eta);
	}
	return tabFKBC;
}

//wszystkie te tablice musimy pomnozyc razy konwekcje

double*** LocalExecute::mnozenieMacierzyAlfa(int size, double ***tab, double alfa)
{
	macierzAlfa = new double**[8];
	for (int i = 0; i < 8; i++)
	{
		macierzAlfa[i] = new double*[8];
	}
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < size; j++)
		{
			macierzAlfa[i][j] = new double[size];
		}
	}
	for (int l = 0; l < size; l++)
	{
		int k = 0;
		int k1 = 0;
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				macierzAlfa[k1][l][0] = tab[k][j][l] * tab[k][j][0] * alfa;
				macierzAlfa[k1][l][1] = tab[k][j][l] * tab[k][j][1] * alfa;
				macierzAlfa[k1][l][2] = tab[k][j][l] * tab[k][j][2] * alfa;
				macierzAlfa[k1][l][3] = tab[k][j][l] * tab[k][j][3] * alfa;
				k1++;
			}
			k++;
		}
	}

	return macierzAlfa;
}

//macierz mnozymy razy wyznacznik jacobianu 1d po powierzchni

double*** LocalExecute::mnozenieMacierzyAlfaDet(int size, double ***tab)
{
	macierzAlfaDet = new double**[size];
	for (int i = 0; i < size; i++)
	{
		macierzAlfaDet[i] = new double*[size];
	}
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			macierzAlfaDet[i][j] = new double[size];
		}
	}
	int k = 0;
	for (int l = 0; l < size; l++)
	{
		for (int i = 0; i < size; i++)
		{
			macierzAlfaDet[l][i][0] = (tab[k][i][0] + tab[k + 1][i][0]) * tabDetJBoku[l];
			macierzAlfaDet[l][i][1] = (tab[k][i][1] + tab[k + 1][i][1]) * tabDetJBoku[l];
			macierzAlfaDet[l][i][2] = (tab[k][i][2] + tab[k + 1][i][2]) * tabDetJBoku[l];
			macierzAlfaDet[l][i][3] = (tab[k][i][3] + tab[k + 1][i][3]) * tabDetJBoku[l];
		}
		k = k + 2;
	}
	return macierzAlfaDet;
}

//wyznaczamy macierz H z warunkami brzegowymi, warunki brzegowe ustalam w gridzie

double** LocalExecute::macierzHBC(int size)
{
	tabMacierzHBC = new double*[size];
	for (int i = 0; i < size; i++)
	{
		tabMacierzHBC[i] = new double[size];
	}
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			tabMacierzHBC[i][j] = (element->nodeTab[0].BC & element->nodeTab[1].BC) * macierzAlfaDet[0][i][j] +
				(element->nodeTab[1].BC & element->nodeTab[2].BC) * macierzAlfaDet[1][i][j] +
				(element->nodeTab[2].BC & element->nodeTab[3].BC) * macierzAlfaDet[2][i][j] +
				(element->nodeTab[3].BC & element->nodeTab[0].BC) * macierzAlfaDet[3][i][j];
		}
	}
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			tabMacierzHBC[i][j] += tabMacierzH[i][j];
		}
	}

	cout << "Macierz HBC:" << endl;
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			cout << tabMacierzH[i][j] << "  ";
		}
		cout << endl;
	}

	return tabMacierzHBC;
}

//wektor P jest to wektor obciazen, musimy ustalic odpowiednie warunki brzegowe, dlugosci bokow i na tej podstawie wyznaczniki jacobianu 1d
//wspolczynnik wymiany alfa, temperature otoczenia i funkcje ksztaltu dla elementu uniwersalnego

double* LocalExecute::wektorP(int size, double alfa, double t)
{
	double sum = 0;
	tabWektorP = new double[size];
	tabMacierzWektorP = new double*[size];

	for (int i = 0; i < size; i++)
	{
		tabMacierzWektorP[i] = new double[size];
	}

	for (int n = 0; n < size; n++)
	{
		tabMacierzWektorP[n][0] = (alfa * tabDetJBoku[n] * t * tabFKBC[n][0][0]) + (alfa * tabDetJBoku[n] * t * tabFKBC[n][1][0]);
		tabMacierzWektorP[n][1] = (alfa * tabDetJBoku[n] * t * tabFKBC[n][0][1]) + (alfa * tabDetJBoku[n] * t * tabFKBC[n][1][1]);
		tabMacierzWektorP[n][2] = (alfa * tabDetJBoku[n] * t * tabFKBC[n][0][2]) + (alfa * tabDetJBoku[n] * t * tabFKBC[n][1][2]);
		tabMacierzWektorP[n][3] = (alfa * tabDetJBoku[n] * t * tabFKBC[n][0][3]) + (alfa * tabDetJBoku[n] * t * tabFKBC[n][1][3]);
	}

	for (int i = 0; i < 4; ++i) {
		sum += tabMacierzWektorP[0][i] * (element->nodeTab[0].BC & element->nodeTab[1].BC) +
			tabMacierzWektorP[1][i] * (element->nodeTab[1].BC & element->nodeTab[2].BC) +
			tabMacierzWektorP[2][i] * (element->nodeTab[2].BC & element->nodeTab[3].BC) +
			tabMacierzWektorP[3][i] * (element->nodeTab[3].BC & element->nodeTab[0].BC);


		tabWektorP[i] = sum;


		sum = 0;

	}

	cout << "Wektor P:" << endl;
	for (int i = 0; i < size; i++)
	{
		cout << tabWektorP[i] << "  ";
	}
	cout << endl;
	cout << endl;
	return tabWektorP;
}

//OSOBNE FUNCKJE WYKONUJACE POPRZEDNIE FUNCKJE

void LocalExecute::funkcjeKsztaltu()
{
	File file;
	int size = File::getPointCount();
	double integrationPoint = File::getIntegrationPoint();

	setTabKsi(size, integrationPoint);
	setTabEta(size, integrationPoint);
	double **tabN = N(size);
	pochodnaFunkcjiKsztaltuEta(size);
	pochodnaFunkcjiKsztaltuKsi(size);
}

void LocalExecute::matrixJ()
{
	int size = file.getPointCount();

	obliczJacobian(size);
	detJacobian(size);
	obliczJacobianPrzezDet(size);

}

void LocalExecute::matrixH()
{
	int size = File::getPointCount();
	double conductivity = File::getConductivity();
	double convection = File::getConvection();

	jacobianDx = pochodnaJacobianuPoDx(size);
	jacobianDy = pochodnaJacobianuPoDy(size);

	pochodnaFunkcjiKsztaltuKsi(size);
	pochodnaFunkcjiKsztaltuEta(size);
	double*** mnozenieMatrixDx = mnozenieMacierzy(size, jacobianDx);
	double*** mnozenieMatrixDy = mnozenieMacierzy(size, jacobianDy);
	double*** matrixRazyDetDx = macierzRazyDet(size, mnozenieMatrixDx);
	double*** matrixRazyDetDy = macierzRazyDet(size, mnozenieMatrixDy);
	macierzIkondukcja(size, conductivity, matrixRazyDetDx, matrixRazyDetDy);
	macierzH(size);

	dlugoscBokow(size);
	detJBoku(size);
	tabFKBC = funkcjeKsztaltuBC(size);
	macierzAlfa = mnozenieMacierzyAlfa(size, tabFKBC, convection);
	mnozenieMacierzyAlfaDet(size, macierzAlfa);
	macierzHBC(size);
}

void LocalExecute::matrixC()
{
	int size = File::getPointCount();
	double ro = File::getDensity();
	double c = File::getSpecificHeat();
	macierzCalkowania(size, ro, c);
	macierzC(size);

}

void LocalExecute::vectorP()
{
	int size = File::getPointCount();
	double alfa = File::getConvection();
	double tempOt = File::getAmbientTemperature();
	wektorP(size, alfa, tempOt);
}

//int LocalExecute::iterator=0;