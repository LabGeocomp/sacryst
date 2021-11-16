# sacryst
Simulated Annealing With Crystallization Heuristic

Use CMake to create the project.

There are three flags in the CMakeLists:

* BENCHMARK_TEST - Creates the execultables of the benchmark functions.
* STATISTICS_EXIT - Outputs a txt file with the full statistics.
* SCREEN_EXIT - Prints on the screen the intermediates results, such as the cost, the temperature and the number of accepted and rejected solutions.


## Usage

An example to execute the sphere test.

```bash
sacryst_sphere.exe results.txt -l -100 -u 100 --numvars 10 --maxits 100 --maxaccept 50 --inittemp 1e1 --finaltemp 1e-100 --totalits 100000
```

You have to input the following arguments:

* results.txt - the name of the output file;
* -l -100 - the lower limit value. In this example, it is -100;
* -u 100 - the upper limit value. In this example, it is 100;
* --numvars 10 - the number of variables of the objective function. In this example, it is 10;
* --maxits 100 - the maximum number of iterations for a certain temperature. In this example, it is 100;
* --maxaccept 50 - the maximum number of accepted solutions for a certain temperature. In this example, it is 50;
* --inittemp 1e1 - the initial temperature. In this example, it is 10;
* --finaltemp 1e-100 - the final temperature. In this example, it is 1e-100;
* --totalits 100000 - the total number of iterations. In this example, it is 100000;
