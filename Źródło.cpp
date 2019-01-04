#include "File.h"
#include "Element.h"
#include "Grid.h"
#include "Node.h"
#include "LocalExecute.h"
#include <iomanip>

int main()
{
	Grid grid;

	//grid.showNodesCordinates();
	//grid.showElementID();
	grid.setValueofMatrixHAndVectorP();
	grid.showGlobalMatrixC();
	cout << endl;
	grid.showGlobalMatrixH();
	cout << endl;
	grid.showGlobalVectorP();
	
	grid.calculateTemperatures();

	system("pause");
	return 0;
}