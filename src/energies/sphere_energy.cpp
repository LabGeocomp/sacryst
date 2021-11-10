#include "cParameters.h"
#include <vector>

double calculateEnergy(cParameters &params)
{
	//*****************************Sphere*****************************
	std::vector<double>::iterator itt = params.fbegin();
	double val = 0;
	while (itt != params.fend()) {
		val += (*itt) * (*itt);
		itt++;
	}
	return val;
}