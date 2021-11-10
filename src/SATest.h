#ifndef _SATEST_
#define _SATEST_

#include <vector>
#include <iostream>
#include <random>
#include <math.h>  
#include <chrono>
#include <fstream>

#define C_LIMIT 20
#define BOLTZ 1

auto seed_norm = std::chrono::system_clock::now().time_since_epoch().count();
std::mt19937 generator((unsigned int)seed_norm);

double genrand_real1(void) {
	std::uniform_real_distribution<double> distribution(0.0, 1.0);
	return distribution(generator);
}

// --> This class represents the parameters
class cParameters {
public:
	cParameters() {};
	~cParameters() {};

	typedef std::vector<double>::iterator fIterator;
	typedef std::vector<double>::const_iterator  const_fIterator;
	typedef std::vector<int>::iterator iIterator;
	typedef std::vector<int>::const_iterator  const_iIterator;

	// --> Iterators operation
	fIterator fbegin(void) { return vFloats.begin(); };
	const_fIterator fbegin(void) const { return vFloats.begin(); };
	fIterator fend(void) { return vFloats.end(); };
	const_fIterator fend(void) const { return vFloats.end(); };
	iIterator ibegin(void) { return vInts.begin(); };
	const_iIterator ibegin(void) const { return vInts.begin(); };
	iIterator iend(void) { return vInts.end(); };
	const_iIterator iend(void) const { return vInts.end(); };


	// --> Retrive the num-th integer parameter
	int getIntParameter(int num) {
		if (num < 0 || num >(int)(vInts.size() - 1)) return 0;
		int cont = 0; iIterator it = ibegin();
		while (cont != num && it != iend()) {
			cont++; it++;
		}
		return *it;
	};

	// --> Set the num-th integer parameter
	void setIntParameter(int num, int param) {
		if (num < 0 || num >(int)(vInts.size() - 1)) return;
		int cont = 0; iIterator it = ibegin();
		while (cont != num && it != iend()) {
			cont++; it++;
		}
		*it = param;
	};

	// --> Retrive the num-th flaot parameter
	double getFloatParameter(int num) {
		if (num < 0 || num >(int)(vFloats.size() - 1)) return 0.0;
		int cont = 0; fIterator it = fbegin();
		while (cont != num && it != fend()) {
			cont++; it++;
		}
		return *it;
	};

	// --> Set the num-th Float parameter
	void setFloatParameter(int num, double param) {
		if (num < 0 || num >(int)(vFloats.size() - 1)) return;
		int cont = 0; fIterator it = fbegin();
		while (cont != num && it != fend()) {
			cont++; it++;
		}
		*it = param;
	};

	int getNumberFloats(void) const { return (int)vFloats.size(); };
	int getNumberInts(void) const { return (int)vInts.size(); };

	void clear(void) {
		vFloats.clear();
		vInts.clear();
		crystallization.clear();
		accepted.clear();
		rejected.clear();
	};

	void addFloat(double val) { vFloats.push_back(val); crystallization.push_back(0); accepted.push_back(0); rejected.push_back(0); };
	void addInt(int n) { vInts.push_back(n); };

	void setMax(std::vector<double> max, int imax) { this->max = max; this->imax = imax; };
	std::vector<double> getMax(void) const { return this->max; };
	void setMin(std::vector<double> min, int imin) { this->min = min; this->imin = imin; };
	std::vector<double> getMin(void) const { return this->min; };

	void pushMax(double val){
		max.push_back(val);
	};

	void pushMin(double val){
		min.push_back(val);
	};

	double getMin(int index) const {
		return min[0];
	};

	double getMax(int index) const {
		return max[0];
	};

	double getStep(int index) const {
		return max[0] - min[0];
	};

	int getLast(void) const { return this->last; };
	void setLast(int last) { this->last = last; };
	int getFirst(void) const { return this->first; };
	void setFirst(int first) { this->first = first; };

	void incrementCrystallization(int index) {
		if (index + 1 > 0)
			if (crystallization[index] < 5000000)
			crystallization[index] ++;

		rejected[index] ++;

	};
	void resetCrystallization(int index) {
		if (index + 1 > 0) {
			if (crystallization[index] > 0) {
				crystallization[index] = 0 ;
			}
			if (crystallization[index] < 0)
				crystallization[index] = 0;
		}
		accepted[index] ++;
	};

	void decreaseCrystallization(int index, int pass) {
		if (index + 1 > 0) {
			if (crystallization[index] > 0) {
				crystallization[index] -=pass;
			}
			if (crystallization[index] < 0)
				crystallization[index] = 0;
		}
		accepted[index] ++;
	};

	int getCrystallization(int index) { return crystallization[index]; };

	int getaccepted(int index) { return accepted[index]; };

	int getrejected(int index) { return rejected[index]; };

	void verifycrystalization() {
		int i = 1000;
		while (i < vFloats.size()) {
			if (accepted[i] > rejected[i]) {
				crystallization[i] /=2;
			}
			accepted[i] = 0;
			rejected[i] = 0;
			i++;
		}
	}

	double shuffle(int index) {
		double val = getMax(index);
		int count = 0;
		if (crystallization[index] < C_LIMIT) {
			while (val + vFloats[index] < getMin(index) || val + vFloats[index] > getMax(index) || count == 0) {
				val = 0;
				for (int k = 0; k < crystallization[index] + 1; k++)
					val += genrand_real1();
				val /= (double)(crystallization[index] + 1);
				val -= 0.5;
				val *= getStep(index);
				count = 1;
			}
		}
		else{
			std::normal_distribution<double> distribution(0.0, 1 / (double)(1 * exp(crystallization[index] + 2  - C_LIMIT)));
			while (val + vFloats[index] < getMin(index) || val + vFloats[index] > getMax(index) || count == 0) {
				val = distribution(generator)*getStep(index);
				count = 1;
			}
		}

		vFloats[index] += val;
		return val;
	};



private:

	// --> Vector of Floats
	std::vector<double> vFloats;

	// --> Vector of Integers
	std::vector<int> vInts;

	// --> Parameter range
	std::vector<double> min, max;

	// --> Pointer to first and last points
	int first, last, imin, imax;

	// --> Cristallization factor for this parameter
	std::vector<int> crystallization;

	// --> accepted and rejected by parameter
	std::vector<int> accepted;
	std::vector<int> rejected;
};

class cInterpretParameters {
public:
	cInterpretParameters() {
		parameters = (cParameters *)0;
	};
	~cInterpretParameters() {};

	// --> Configure the parameters
	void setParameters(cParameters *parameters) { this->parameters = parameters; };
	cParameters *getParameters(void) { return this->parameters; };

	
private:
	// --> Pointer to parameters
	cParameters *parameters;
};

class cEditParameters {
public:
	cEditParameters() {
		parameters = (cParameters *)0;
		pEditFloat = -1;
	};
	~cEditParameters() {};

	// --> Configure the parameters
	void setParameters(cParameters *parameters) { this->parameters = parameters; };
	cParameters *getParameters(void) { return this->parameters; };


	// --> Edit next float parameter
	void goNextFloatParameter(void) {
		if (parameters == (cParameters *)0) {
			std::cout << "There is no associated parameters." << std::endl;
			return;
		}

		if (pEditFloat < parameters->getNumberFloats() - 1)
			pEditFloat++;
	};

	// --> Edit previous float parameter
	void goPreviousFloatParameter(void) { if (pEditFloat > 0) pEditFloat--; };

	// --> Modify current float parameter
	void modifyFloatParameter(double step) {
		if (parameters == (cParameters *)0) {
			std::cout << "There is no associated parameters." << std::endl;
			return;
		}

		if (pEditFloat < 0) {
			std::cout << "There is no current float parameter set." << std::endl;
			return;
		}

		double param = parameters->getFloatParameter(pEditFloat);
		parameters->setFloatParameter(pEditFloat, param + step);
	};

	int getCurrentFloatParameter(void) const { return this->pEditFloat; };

private:
	// --> Pointer to parameters
	cParameters *parameters;

	// --> Pointer to the floating parameter that is edited
	int pEditFloat;
};

class cStatistics {
public:
	cStatistics() { first = true; };
	~cStatistics() {};

	void init(void) {
		first = true;
		#ifdef WRITE_STATISTICS
		std::ofstream fp;
		fp.open(filename, std::ios_base::out);
		if (!fp.is_open()) {
        	std::cout << "Failed to open " << filename << " for writing statistics." << std::endl;
			return;
    	}
		fp << "pTa, tMinE, tMaxE, tMedia, tVariance, tSpecifHeat, tVarCost, tAvgSquare, tAvgCost, sum, pAlfa, pAccepted, pNotAccepted, N_iter" << std::endl;
		fp.close();
		#endif
	};

	// --> Retrieve new value for alfa
	double getAlfa(void) const { return this->pAlfa; };

	// --> Set the storage for statistics
	void setFilename(std::string &filename) { this->filename = std::string(filename); };

	// --> Add a new energy value
	void pushEnergy(double value) { energy.push_back(value); };

	// --> clear energy vactor
	void clearEnergy() { energy.clear(); };

	// --> Calculate the statistics
	void print(double pTa, int pAccepted, int pNotAccepted, int N_iter) {

		std::vector<double>::iterator itf = energy.begin();
		double sum = 0, tMinE = *itf, tMaxE = *itf;
		for (itf = energy.begin(); itf != energy.end(); itf++) {
			if (tMinE > *itf) tMinE = *itf;
			if (tMaxE < *itf) tMaxE = *itf;
			sum += exp(-(*itf) /( BOLTZ*(pTa)));
		}

		std::vector<double> PiT;
		for (itf = energy.begin(); itf != energy.end(); itf++)
			PiT.push_back(1.0/energy.size());

		double tAvgCost = 0, tAvgSquare = 0;
		std::vector<double>::iterator itg = energy.begin();
		for (itf = PiT.begin(); itf != PiT.end(); itf++, itg++) {
			tAvgCost += (*itf) * (*itg);
			tAvgSquare += (*itf) * (*itg) * (*itg);
		}

		double tMedia = 0.0; int tCount = 0;
		for (itf = PiT.begin(); itf != PiT.end(); itf++, tCount++) 
			tMedia += *itf;
		tMedia /= (double)tCount;

		double std_amost = 0.0;
		double avg_amost = 0.0;
		for (itg = energy.begin(); itg != energy.end(); itg++)
			avg_amost += *itg;
		tMedia /= (double)tCount;
		for (itg = energy.begin(); itg != energy.end(); itg++)
			std_amost += (*itg- avg_amost)*(*itg - avg_amost);
		std_amost /= (double)(tCount);
		std_amost = sqrt(std_amost);


		double tVariance = 0.0;
		for (itf = PiT.begin(); itf != PiT.end(); itf++)
			tVariance += (*itf - tMedia) * (*itf - tMedia);
		tVariance /= (double)tCount;

		tVarCost = tAvgSquare - tAvgCost * tAvgCost;
		if (first) {
			pLastVarCost = tVarCost;
			pLastVariance = tVariance;
			pLastT = pTa;
			first = false;
		}

		double tVarC = (1.0 - 0.95) * tVarCost + 0.95 * pLastVarCost * pTa / pLastT;
		double tVar = (1.0 - 0.95) * tVariance + 0.95 * pLastVariance * pTa / pLastT;
		double tSpecifHeat = tVarCost / (pTa * pTa);

		pLastT = pTa;
		pLastVariance = tVariance;
		pLastVarCost = tVarCost;
		pAlfa = exp(-(0.05 * pTa) / sqrt(tVarCost));
		if (isnan(pAlfa)) pAlfa = 0.99;
		if (pAlfa < 0.8) pAlfa = 0.8;

		#ifdef WRITE_STATISTICS
		std::ofstream fp;
		fp.open(filename, std::ios_base::out|std::ios_base::app);
		if (!fp.is_open()) {
        	std::cout << "Failed to open " << filename << " for writing appended statistics." << std::endl;;
			return;
    	}
		fp << pTa << "; " << tMinE << "; " << tMaxE << "; " << tMedia << "; " << tVariance << "; " << tSpecifHeat << "; " << tVarCost << "; " << tAvgSquare << "; " << tAvgCost << "; " << sum << "; " << pAlfa << "; " << pAccepted << "; " << pNotAccepted << "; " << N_iter << "; " << std_amost << std::endl;
		fp.close();
		#endif
	};

	double get_variance(void) { return this->tVarCost; };

private:
	std::vector<double> energy;

	// --> File to store statistics
	std::string filename;

	// --> Last variational cost
	double pLastVarCost;

	// --> Last variance
	double pLastVariance;

	// --> Last temperature
	double pLastT;

	// --> New value for alfa
	double pAlfa;

	// --> Flag indicating first time
	bool first;

	// --> energy variance
	double tVarCost;
};

class cSimulatedAnnealing{
public:
	cSimulatedAnnealing() {};
	~cSimulatedAnnealing() {};

	void init(std::vector<double> &lowerLimit, std::vector<double> &upperLimit, int numbvar, int iteractionNumber, int acceptedEquilibrium, double initT, double Final_Temp, int Max_iter, std::string &res, std::string &wfilename) {
		this->Final_Temp = Final_Temp;
		this->Max_iter = Max_iter;

		// --> Reserve space for parameters
		pListParameters.clear();

		// --> Set storage for statistics
		statistics.setFilename(wfilename);

		// --> Configure parameters
		interpreter.setParameters(&pListParameters);

		// --> Read parameter file
		pListParameters.clear();
		for (const auto val: lowerLimit) { pListParameters.pushMin(val); } 
		for (const auto val: upperLimit) { pListParameters.pushMax(val); } 
		for (int aux = 0; aux < numbvar; aux++) pListParameters.addFloat(pListParameters.getMin(aux) + (pListParameters.getMax(aux) - pListParameters.getMin(aux))*genrand_real1());

		pTa = initT;
		pAccepted = pNotAccepted = 0;
		pAlfa = 0.95;

		firstShow = true;

		resFilename = (char*)res.c_str();

		this->pIteractionNumber = iteractionNumber;
		this->pAcceptedEquilibrium = acceptedEquilibrium;

	}

	void setValue(int index, double value) { pListParameters.setFloatParameter(index, value); };

	double getValue(int index) { return pListParameters.getFloatParameter(index); };
	double getEnergyCandidate(void) const { return this->pEnergyCandidate; };


	void setCalculateEnergyFunction(double(*calculateEnergy)(cParameters &)) {
		this->calculateEnergy = calculateEnergy;
	};

	void acceptedCandidate(long index, double cEnergyCandidate) {
		pListParameters = pListParametersCandidates;
		interpreter.setParameters(&pListParameters);
		pEnergyCandidate = cEnergyCandidate;
		pAccepted++;
		statistics.pushEnergy(cEnergyCandidate);
		if (std_temp > 100 ) pListParameters.resetCrystallization(index);
		else {
			pListParameters.decreaseCrystallization(index, 2);
		}
	};

	void run(void) {
		int count = 0;
		double step = 0;
		bestCont = 0;
		N_iter = 0;

		// --> Determine the initial Energy
		interpreter.setParameters(&pListParameters);
		pEnergyCandidate = pEnergyCandidateBest = (*calculateEnergy)(pListParameters);
		pListParametersBest = pListParameters;

		// --> Initiate the statistics calculation
		statistics.init();

		// --> While the stop criteria is not reached
		while (stopCriteria()) {

			// --> Initialize counters
			pAccepted = pNotAccepted = 0;

			// --> While the equilibrium is not reached
			while (equilibriumNotReached()) {
		
				// --> Copy the actual candidate
				pListParametersCandidates = pListParameters;

				// --> Calculate the new shuffling
				long index = (int)(genrand_real1() * ((double)pListParameters.getNumberFloats()));
				step = pListParametersCandidates.shuffle(index);
				interpreter.setParameters(&pListParametersCandidates);

				// --> Calculate the new Energy
				double cEnergyCandidate = (*calculateEnergy)(pListParametersCandidates);

				// --> New energy is lower?
				if (cEnergyCandidate < pEnergyCandidate){  //delta"e" = cEnergyCandidate-pEnergyCandidate (nova energia - a anterior)
					// --> Then accept the candidate
					acceptedCandidate(index, cEnergyCandidate);
				}
				else {
					// --> The probability is lower than the energy probability?
					if (exp((pEnergyCandidate - cEnergyCandidate) / (BOLTZ*(pTa))) > genrand_real1()){
						// --> Then accept the candidate
						acceptedCandidate(index, cEnergyCandidate);
					}
					else {
						// --> The candidate was not accepted
						pListParameters.incrementCrystallization(index);
						pNotAccepted++;
					}
				}

				// --> Check for best ever candidate
				if ((pEnergyCandidateBest > cEnergyCandidate) || firstShow) {
					if (pEnergyCandidateBest > cEnergyCandidate) {
						pListParametersBest = pListParameters;
						pEnergyCandidateBest = cEnergyCandidate;
						bestCont++;
					}
					pListParameters = pListParametersBest;
					interpreter.setParameters(&pListParameters);
					pEnergyCandidate = (*calculateEnergy)(pListParameters);
					firstShow = false;
				}
				N_iter++;
			}
			#ifdef PRINT_INTERMEDIATE_RESULTS
			std::cout << "*** " << pEnergyCandidateBest << ", " << pTa << ", " << pAccepted << ", " << (pAccepted + pNotAccepted) << ", " << N_iter << std:: endl;
			#endif
			statistics.print(pTa, pAccepted, pNotAccepted, N_iter);

			// --> Retrieve new value for alfa
			pAlfa = statistics.getAlfa();

			// --> Go to next temperature
			if (pTa > pow(10,-1000)){
				pTa = nextTemperature(); };

			// --> purge energy vector
			statistics.clearEnergy();

			std_temp  = sqrt(statistics.get_variance());
		}
		#ifdef PRINT_INTERMEDIATE_RESULTS
		std::cout << "End = " << pEnergyCandidateBest << ", " << pTa << ", " << pAccepted << ", " << (pAccepted + pNotAccepted) << ", " << N_iter << std::endl;
		#endif
	};

	void write(const std::string filename) {
		interpreter.setParameters(&pListParametersBest);
		std::ofstream fw;
		fw.open(filename, std::ios_base::out);
		if (!fw.is_open()) {
        	std::cout << "Failed to open " << filename << " for writing results." << std::endl;;
			return;
    	}
		fw << "*** Exit File " << std::endl;
		fw << std::endl << std::endl;
		fw << "***  Variable Values Crystallization NUmber of Acceptance and Rejections " << std::endl;
		std::vector<double>::iterator itt = pListParametersBest.fbegin();
		int i = 0;
		for(auto itt = pListParametersBest.fbegin(); itt != pListParametersBest.fend(); itt++) {
			fw  << *itt << " " << pListParametersBest.getCrystallization(i) << " " << pListParametersBest.getaccepted(i) << " " << pListParametersBest.getrejected(i) << std::endl;
			i++;
		}
		fw << "---" << std::endl;
		fw << std::endl << std::endl;
		fw << "*** Cost" << std::endl;
		fw << pEnergyCandidateBest << std::endl;
		fw << std::endl << std::endl;
		fw << "*** Number of iterations" << std::endl;
		fw << N_iter << std::endl;
		fw.close();
	};

	// --> Criteria for stop
	bool stopCriteria(void) {
		return (pTa > Final_Temp || firstShow) && N_iter < Max_iter;
	};

	// --> Determine the next temperature
	double nextTemperature(void) { return pTa * pAlfa;};

	// --> Criterium for equilibrium
	bool equilibriumNotReached(void) {
		return pAccepted + pNotAccepted < pIteractionNumber && pAccepted < pAcceptedEquilibrium && N_iter < Max_iter;
	};

	// --> Total number of accepted solutions
	long pAccepted;

	// --> Total number of not accepted solutions
	long pNotAccepted;

	// --> Number of accepted to reach equilibrium
	long pAcceptedEquilibrium;

	// --> Total number of iteractions to reach equilibrium
	long pIteractionNumber;

	// --> Total number of iteractions
	long N_iter;

	// --> Cooling factor
	double pAlfa;

	// --> Actual temperature
	double pTa;

	// --> Final temperature and iterations
	double Final_Temp;
	int Max_iter;

	// --> Energy for the actual candidate
	double pEnergyCandidate;

	// --> List of parameters
	cParameters pListParameters;
	cParameters pListParametersCandidates;
	cParameters pListParametersBest;


	// --> List of parameters for best energy
	double pEnergyCandidateBest;

	// --> Store the statistics
	cStatistics statistics;

	// --> Pointer for energy calculation
	double(*calculateEnergy)(cParameters &);

	cInterpretParameters interpreter;

	bool firstShow;

	char *resFilename;

	long bestCont;

	double std_temp;

	int std_mag;

};

#endif	/* _SATEST_ */
