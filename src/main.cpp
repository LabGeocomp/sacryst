#include <algorithm>
#include "include/args.hxx"
#include "cParameters.h"
#include "SATest.h"

double calculateEnergy(cParameters &params);

cSimulatedAnnealing SA;

int main(int argc, char **argv)
{
	// Parse coomand line arguments
	args::ArgumentParser parser("Benchmark for the simulated annealing with crystallization algorithm.", "No comments.");
	args::HelpFlag help(parser, "help", "Display this help menu", { 'h', "help" });
	args::CompletionFlag completion(parser, { "complete" });
	args::Positional<std::string> outfile(parser, "filename", "Output file name");
	args::ValueFlagList<double> lower(parser, "number", "Lower limit", { 'l' });
	args::ValueFlagList<double> upper(parser, "number", "Upper limit", { 'u' });
	args::ValueFlag<int> numvars(parser, "number", "Number of variables", { "numvars" });
	args::ValueFlag<int> maxits(parser, "number", "Maximum iteration per temperature", { "maxits" });
	args::ValueFlag<int> maxaccept(parser, "number", "Maximum accepted solutions per temperature", { "maxaccept" });
	args::ValueFlag<double> inittemp(parser, "number", "Initial temperature", { "inittemp" });
	args::ValueFlag<double> finaltemp(parser, "number", "Final temperature", { "finaltemp" });
	args::ValueFlag<int> totalits(parser, "number", "Maximum number of iterations", { "totalits" });
	try
    {
        parser.ParseCLI(argc, argv);
    }
    catch (args::Help)
    {
        std::cout << parser;
        return 0;
    }
    catch (args::ParseError e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    catch (args::ValidationError e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
	if (outfile) { std::cout << "outfile: " << args::get(outfile) << std::endl; }
	else { std::cerr << "Missing output results file name"; return 1;}
    if (lower) { for (const auto l: args::get(lower)) { std::cout << "l: " << l << std::endl; } }
	else { std::cerr << "Missing lower limits"; return 1;}
    if (upper) { for (const auto u: args::get(upper)) { std::cout << "u: " << u << std::endl; } }
	else { std::cerr << "Missing upper limits"; return 1;}
    if (numvars) { std::cout << "numvars: " << args::get(numvars) << std::endl; }
	else { std::cerr << "Missing number of variables (numvars)"; return 1;}
    if (maxits) { std::cout << "maxits: " << args::get(maxits) << std::endl; }
	else { std::cerr << "Missing maximum iteration per temperature (maxits)"; return 1;}
    if (maxaccept) { std::cout << "maxaccept: " << args::get(maxaccept) << std::endl; }
	else { std::cerr << "Missing maximum accepted solutions per temperature (maxaccept)"; return 1;}
    if (inittemp) { std::cout << "inittemp: " << args::get(inittemp) << std::endl; }
	else { std::cerr << "Missing initial temperature (inittemp)"; return 1;}
    if (finaltemp) { std::cout << "finaltemp: " << args::get(finaltemp) << std::endl; }
	else { std::cerr << "Missing final temperature (finaltemp)"; return 1;}
    if (totalits) { std::cout << "totalits: " << args::get(totalits) << std::endl; }
	else { std::cerr << "Missing maximum number of iterations (totalits)"; return 1;}

	std::vector<double> lowerVec = args::get(lower);
	std::vector<double> upperVec = args::get(upper);
	std::string outfileStr = args::get(outfile);
	std::string statoutfile("statistics.txt");
	SA.init(lowerVec, upperVec, args::get(numvars), args::get(maxits), args::get(maxaccept), args::get(inittemp), args::get(finaltemp), args::get(totalits), outfileStr, statoutfile);

	SA.setCalculateEnergyFunction(calculateEnergy);
	std::cout << "RUN" << std::endl;
	SA.run();

	SA.write(outfileStr.c_str());
	
    return 0;
}


#ifdef CUSTOM_ENERGY
// --> Definition of the cost function. 
double calculateEnergy(cParameters &params)
{
	std::vector<double>::iterator itt = params.fbegin();
	double val = 0;
	while (itt != params.fend()){
		val += (*itt)*(*itt);
		itt++;
	}
	return val;
}
#endif