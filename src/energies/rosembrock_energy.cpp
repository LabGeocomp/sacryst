#include "cParameters.h"
#include <vector>

double calculateEnergy(cParameters &params)
{
	//*****************************Rosembrock*****************************
	std::vector<double>::iterator itt = params.fbegin();
	std::vector<double>::iterator ittn = params.fbegin() + 1;
	double val = 0;
	while (itt != params.fend() - 1) {
		val += 100 * (*ittn - (*itt) * (*itt)) * (*ittn - (*itt) * (*itt)) + (1 - (*itt) * (*itt)) * (1 - (*itt) * (*itt));
		itt++;
		ittn++;
	}
	return val;
}

