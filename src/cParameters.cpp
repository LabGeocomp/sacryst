#include "cParameters.h"
#include <vector>
#include <random>

#define C_LIMIT 20

extern std::mt19937 generator;
double genrand_real1(void);

cParameters::cParameters() {}
cParameters::~cParameters() {}

// --> Retrive the num-th integer parameter
int cParameters::getIntParameter(int num) {
    if (num < 0 || num >(int)(vInts.size() - 1)) return 0;
    int cont = 0; iIterator it = ibegin();
    while (cont != num && it != iend()) {
        cont++; it++;
    }
    return *it;
}

// --> Set the num-th integer parameter
void cParameters::setIntParameter(int num, int param) {
    if (num < 0 || num >(int)(vInts.size() - 1)) return;
    int cont = 0; iIterator it = ibegin();
    while (cont != num && it != iend()) {
        cont++; it++;
    }
    *it = param;
}

// --> Retrive the num-th flaot parameter
double cParameters::getFloatParameter(int num) {
    if (num < 0 || num >(int)(vFloats.size() - 1)) return 0.0;
    int cont = 0; fIterator it = fbegin();
    while (cont != num && it != fend()) {
        cont++; it++;
    }
    return *it;
}

// --> Set the num-th Float parameter
void cParameters::setFloatParameter(int num, double param) {
    if (num < 0 || num >(int)(vFloats.size() - 1)) return;
    int cont = 0; fIterator it = fbegin();
    while (cont != num && it != fend()) {
        cont++; it++;
    }
    *it = param;
}

void cParameters::clear(void) {
    vFloats.clear();
    vInts.clear();
    crystallization.clear();
    accepted.clear();
    rejected.clear();
}

void cParameters::addFloat(double val) { 
    vFloats.push_back(val); 
    crystallization.push_back(0); 
    accepted.push_back(0); 
    rejected.push_back(0); 
}

void cParameters::addInt(int n) {
    vInts.push_back(n);
}

void cParameters::incrementCrystallization(int index) {
    if (index + 1 > 0)
        if (crystallization[index] < 5000000)
        crystallization[index] ++;

    rejected[index] ++;

}

void cParameters::resetCrystallization(int index) {
    if (index + 1 > 0) {
        if (crystallization[index] > 0) {
            crystallization[index] = 0 ;
        }
        if (crystallization[index] < 0)
            crystallization[index] = 0;
    }
    accepted[index] ++;
}

void cParameters::decreaseCrystallization(int index, int pass) {
    if (index + 1 > 0) {
        if (crystallization[index] > 0) {
            crystallization[index] -=pass;
        }
        if (crystallization[index] < 0)
            crystallization[index] = 0;
    }
    accepted[index] ++;
}

void cParameters::verifycrystalization() {
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

double cParameters::shuffle(int index) {
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
}
