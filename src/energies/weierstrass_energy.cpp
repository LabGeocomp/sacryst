#include "cParameters.h"
#include <vector>

#define M_PI       3.14159265358979323846

double calculateEnergy(cParameters &params)
{
	//*****************************Weierstrass*****************************
	std::vector<double>::iterator itt = params.fbegin();
	double d = params.getNumberFloats();
	double resp = 0;
	double val = 0;
	double val1 = 0;

	while (itt != params.fend()) {
		for (int i = 1; i <= 20; i++) {
			val += pow(0.5, double(i)) * cos(2 * M_PI * pow(3, double(i)) * ((*itt) + 0.5));
		}
		itt++;
	}
	for (int j = 1; j <= 20; j++) {
		val1 += pow(0.5, double(j)) * cos(M_PI * pow(3, double(j)));
	}
	resp = val - d * val1;
	return resp;
}