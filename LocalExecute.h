#pragma once
#include"File.h"
#include"Element.h"

class LocalExecute
{
	File file;

	double *tabEta;
	double *tabKsi;
	double **tabN;
	double *tabXp;
	double *tabYp;
	double **tabFunkcjiKsztaltuKsi;
	double **tabFunkcjiKsztaltuEta;
	double **jacobian;
	double *tabDetJ;
	double **jacobianPrzezDet;
	double **jacobianDx;
	double **jacobianDy;
	double ***macierzTransponowana;
	double ***macierzDet;
	double ***macierzConduction;
	double **tabMacierzH;
	double **tabMacierzC;
	double ***tabMacierzyCalkowania;
	double *tabDlugosciBokow;
	double *tabDetJBoku;
	double ***tabFKBC;
	double ***macierzAlfa;
	double ***macierzAlfaDet;

	double **tabMacierzHBC;
	double **tabMacierzWektorP;
	double *tabWektorP;

	static int iterator;

public:
	Element *element;
	double **tabElement;
	void showElement();
	LocalExecute();
	LocalExecute(Element*);
	~LocalExecute();
	void setElement(int);
	void setTabEta(int, double);
	void setTabKsi(int, double);
	double* getTabEta();
	double* getTabKsi();
	double N1(double, double);
	double N2(double, double);
	double N3(double, double);
	double N4(double, double);
	double** N(int);
	double FKszt1Eta(double);
	double FKszt2Eta(double);
	double FKszt3Eta(double);
	double FKszt4Eta(double);
	double FKszt1Ksi(double);
	double FKszt2Ksi(double);
	double FKszt3Ksi(double);
	double FKszt4Ksi(double);
	double calculateJacobianX(double*);
	double calculateJacobianY(double*);
	double detJ(double*);
	double *detJacobian(int);
	double **obliczJacobian(int);

	double** obliczJacobianPrzezDet(int);
	double** pochodnaFunkcjiKsztaltuKsi(int);
	double** pochodnaFunkcjiKsztaltuEta(int);
	double** pochodnaJacobianuPoDx(int);
	double** pochodnaJacobianuPoDy(int);
	double*** mnozenieMacierzy(int, double**);
	double*** macierzRazyDet(int, double***);
	double*** macierzIkondukcja(int, double, double***, double***);
	double** macierzH(int);

	double*** macierzCalkowania(int, double, double);
	double** macierzC(int);

	double* dlugoscBokow(int);
	double* detJBoku(int);
	double*** funkcjeKsztaltuBC(int);
	double*** mnozenieMacierzyAlfa(int, double***, double);
	double*** mnozenieMacierzyAlfaDet(int, double***);
	double** macierzHBC(int);

	double* wektorP(int, double, double);

	void matrixJ();
	void matrixH();
	void matrixC();
	void vectorP();
	void funkcjeKsztaltu();
};