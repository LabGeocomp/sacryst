#include "SATest.h"
#include <algorithm>

using namespace std;

double calculateEnergy(cParameters &params);

cSimulatedAnnealing SA;

int main(int argc, char **argv)
{
	char filename[500];
	if (argc > 2)
		strcpy_s(filename,500, argv[2]);
	if (argc < 2) {
		printf("It is necessary to write the name of the file!!!!!!\n");
		int c = 0;
		while (true) {
			if (_kbhit()) {
				char key = _getch();
				if (key == ';')
					break;

				filename[c] = key;
				c++;
				printf("%c", key);
			}
		}
		printf("\n");
		filename[c] = 0;
	}
	else
		strcpy_s(filename, argv[1]);

	// --> Configure filename and number of parameters
	if (argc < 3)
		SA.init(filename, "results.txt", "statistics.txt");
	else
		SA.init(filename, argv[2], "statistics.txt");

	SA.setCalculateEnergyFunction(calculateEnergy);
	cout << "RUN" << endl;
	SA.run();

	if (argc < 3)
		SA.write("results.txt");
	else
		SA.write(argv[2]);

	// --> used for saving the best cost and number of iteration for benchmak.
	//SA.write_best("Cost_benchmark.txt");

	//system("pause");
}


// --> Definition of the cost function. 
double calculateEnergy(cParameters &params)
{
	vector<double>::iterator itt = params.fbegin();
	double val = 0;
	while (itt != params.fend()){
		val += (*itt)*(*itt);
		itt++;
	}
	return val;

	//vector<double>::iterator itt = params.fbegin();
	//vector<double>::iterator ittn = params.fbegin() + 1;
	//double val = 0;
	//while (itt != params.fend() - 1) {
	//	val += 100 * (*ittn - (*itt) * (*itt)) * (*ittn - (*itt) * (*itt)) + (1 - (*itt) * (*itt)) * (1 - (*itt) * (*itt));
	//	itt++;
	//	ittn++;
	//}
	//return val;

}


void cSimulatedAnnealing::checkKeyBoard(void)
{
	if (_kbhit()) {
		char key = _getch();
		if (key == 'q') {
			SA.write(resFilename);
			SA.write_best("Cost_benchmark.txt");
			exit(1);
		}
	}
}
