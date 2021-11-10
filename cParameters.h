#ifndef _CPARAMETERS_
#define _CPARAMETERS_

#include <vector>

// --> This class represents the parameters
class cParameters {
public:
	cParameters();
	~cParameters();

	
	int getIntParameter(int num); // --> Retrive the num-th integer parameter
	void setIntParameter(int num, int param); // --> Set the num-th integer parameter
	double getFloatParameter(int num); // --> Retrive the num-th flaot parameter
	void setFloatParameter(int num, double param); // --> Set the num-th Float parameter
	void clear(void);
	void addFloat(double val);
	void addInt(int n);

	void incrementCrystallization(int index);
	void resetCrystallization(int index);
	void decreaseCrystallization(int index, int pass);
	int getCrystallization(int index) { return crystallization[index]; };
	int getaccepted(int index) { return accepted[index]; };
	int getrejected(int index) { return rejected[index]; };
	void verifycrystalization();
	double shuffle(int index);

	int getNumberFloats(void) const { return (int)vFloats.size(); };
	int getNumberInts(void) const { return (int)vInts.size(); };

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

	void setMax(std::vector<double> max) { this->max = max; };
	std::vector<double> getMax(void) const { return this->max; };
	void setMin(std::vector<double> min) { this->min = min; };
	std::vector<double> getMin(void) const { return this->min; };

	void pushMax(double val){ max.push_back(val); };
	void pushMin(double val){ min.push_back(val); };
	// TODO: use individual limits
	double getMin(int index) const { return min[0]; };
	double getMax(int index) const { return max[0]; };
	double getStep(int index) const { return max[0] - min[0]; };

private:

	// --> Vector of Floats
	std::vector<double> vFloats;

	// --> Vector of Integers
	std::vector<int> vInts;

	// --> Parameter range
	std::vector<double> min, max;

	// --> Cristallization factor for this parameter
	std::vector<int> crystallization;

	// --> accepted and rejected by parameter
	std::vector<int> accepted;
	std::vector<int> rejected;
};

#endif /* _CPARAMETERS_ */