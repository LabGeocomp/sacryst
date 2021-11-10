#include "cParameters.h"
#include <vector>

double calculateEnergy(cParameters &params)
{
	//*****************************Griewangk*****************************
	std::vector<double>::iterator itt = params.fbegin();
	double d = params.getNumberFloats();
	double val = 0;
	double val1 = 1;
	int aux = 1;
	//cout << sqrt(double(aux)) << endl;
	while (itt != params.fend()) {
		val += (*itt) * (*itt) / 4000;
		val1 *= cos((*itt) / sqrt(double(aux)));
		//cout <<"HERE!!!"<< *itt << "  " <<SA.getValue(0) <<"   "<< val << endl;
		itt++;
		aux++;
	}
	val -= val1;
	val += 1;
	return val;
}

