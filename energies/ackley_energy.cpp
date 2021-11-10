#include "cParameters.h"
#include <vector>

#define M_PI       3.14159265358979323846

double calculateEnergy(cParameters &params)
{
	//*****************************Ackley*****************************
	std::vector<double>::iterator itt = params.fbegin();
	double d = params.getNumberFloats();
	double resp = 0;
	double val = 0;
	double val1 = 0;
	int aux = 2;
	//cout << sqrt(double(aux)) << endl;
	while (itt != params.fend()) {
		val += cos(2 * M_PI * (*itt));
		val1 += (*itt) * (*itt);
		itt++;
	}
	resp = 20 - exp((1 / d) * val) + exp(1) - 20 * exp(-0.2 * sqrt((1 / d) * val1));
	return resp;
}

