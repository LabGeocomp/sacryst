#include "cParameters.h"
#include <vector>

double calculateEnergy(cParameters &params)
{
    //*****************************Zakharov*****************************
	std::vector<double>::iterator itt = params.fbegin();
	double d = params.getNumberFloats();
	double resp = 0;
	double val = 0;
	double val1 = 0;
	int aux = 1;

	while (itt != params.fend()) {
		val += (*itt) * (*itt);
		val1 += (double(aux) / 2) * (*itt);
		itt++;
		aux++;
	}

	resp = val + val1 * val1 + val1 * val1 * val1 * val1;
	return resp;
}