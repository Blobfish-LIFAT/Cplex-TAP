#include "instance.h"
#include "solver.h"
#include "SolverMPLS.h"
#include "combo.h"
#include <iostream>
#include <string>
#include <math.h>


// combo will not compile with sdl check sue to several C4703 warnings when initializing stuff in switches
// item : profit/wheight/keep
bool* exclude_nodes(cplex_tap::Instance tap, int lb_interest) {
	std::cout << "Filtering Nodes" << endl;
	int maxDist = 0;
	for (auto i = 0; i < tap.size(); ++i) {
		for (auto j = i + 1; j < tap.size(); ++j) {
			int d = tap.dist(i, j);
			if (d > maxDist)
				maxDist = d;
		}
	}
	maxDist++;

	// Debug for combo
	item items[3];
	items[0] = { 2, 2, false };
	items[1] = { 5, 3, false };
	items[2] = { 1, 5, false };

	item* start = items;
	item* end = (item*)(&items + 1) - 1;
	// run combo (pointer to first item, pointer to last item, max whieght, starting lb, starting ub, set x in items, no relaxation)
	stype s = combo(start, end, 6, 0, 99, true, false);

	
	bool out[3]; // out[tap.size()]
	for (auto i = 0; i < 3; ++i) { //i < tap.size()
		out[i] = items[i].x;
		std::cout << out[i];
	}
	std::cout << endl;
	return out;
}

int run_debug(bool progressive, double temps, double dist, std::string path) {
	using namespace cplex_tap;
	const auto tap = Instance(path);

	const auto solver = Solver(tap);

	// Easy
    int budget = lround(temps * tap.size() * 27.5f);
    int dist_bound = lround( dist * tap.size() * 4.5);

	double time = solver.solve_and_print(dist_bound, budget, progressive, false, false);
    std::cout << endl << "TIME TO SOLVE " << time << endl;


	return 0;
}

int run_debug_mpls(double temps, double dist, std::string path) {
    using namespace cplex_tap;
    const auto tap = Instance(path);

    const auto solver = SolverMPLS(tap);

    int budget = lround(temps * tap.size() * 27.5f);
    int dist_bound = lround( dist * tap.size() * 4.5);

    double time = solver.solve_and_print(dist_bound, budget, false, false, false);
    std::cout << endl << "TIME TO SOLVE " << time << endl;


    return 0;
}

int run_epsilon_test(char* argv[]) {
    int series_id = stoi(argv[1]);
	using namespace cplex_tap;

	ofstream res;
	std::stringstream outname;
	outname << "/users/21500078t/res_cplex" << series_id << ".csv";
	res.open(outname.str());
	res << "series_id;size;epsilon_time;epsilon_distance;time_solve" << endl;

	float eptimes[] = {0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90};
	float epdists[] = {0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90};
	int sizes[] = {100, 300, 500, 700};

    for(const int &size : sizes){
	    for (const float &epdist : epdists){
            for (const float &eptime : eptimes){
                std::cout << "Loading TAP instance " << size << endl;
                std::stringstream fname;
                fname << "/users/21500078t/cplex_test/instances/tap_" << series_id << "_" << size << ".dat";
                const auto tap = Instance(fname.str());

                const auto solver = Solver(tap);

                int budget = lround(eptime * size * 27.5f);
                int dist_bound = lround( epdist * size * 4.5);

                double time = solver.solve_and_print(dist_bound, budget, false, false, false);
                std::cout << endl << "TIME TO SOLVE " << time << endl;
                res << series_id << ";" << size << ";" << eptime << ";" << epdist << ";" << time << endl;
                res.flush();
            }
        }
    }

	res.flush();
	res.close();
	return 0;
}

int production(char* argv[]) {
    using namespace cplex_tap;
    const auto tap = Instance(argv[1]);
    const auto solver = Solver(tap);
    int budget = stoi(argv[2]);
    int dist_bound = stoi(argv[3]);

    double time = solver.solve_and_print(dist_bound, budget, false, false, true);

    return 0;
}

int main(int argc, char* argv[]) {
    return production(argv);
    //return run_epsilon_test(argv);
	//return run_debug(false, 0.15, 0.20,"/users/21500078t/cplex_test/instances/tap_8_500.dat");
	//return run_debug_mpls(0.15, 0.20,"/users/21500078t/cplex_test/instances/tap_12_300.dat");
}

