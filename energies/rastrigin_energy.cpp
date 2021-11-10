#include "cParameters.h"
#include <vector>

#define M_PI       3.14159265358979323846

double calculateEnergy(cParameters &params)
{
	//*****************************Rastrigin*****************************
	std::vector<double>::iterator itt = params.fbegin();
	double d = params.getNumberFloats();
	double val = 0;
	while (itt != params.fend()) {
		val += (*itt) * (*itt) - 10 * cos(2 * M_PI * (*itt));
		itt++;
	}
	val += 10 * d;
	return val;
}

