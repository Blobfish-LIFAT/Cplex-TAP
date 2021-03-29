#include "instance.h"
#include "solver.h"
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

int run_debug() {
	using namespace cplex_tap;

	const auto tap = Instance("/users/21500078t/cplex_test/instances/tap_1_100.dat");

	const auto solver = Solver(tap);

	double time = solver.solve_and_print(200, 50, false);
    std::cout << endl << "TIME TO SOLVE " << time << endl;


	return 0;
}

int run_scale_test(int series_id) {
	using namespace cplex_tap;

	ofstream res;
	std::stringstream outname;
	outname << "/users/21500078t/res_cplex" << series_id << ".csv";
	res.open(outname.str());
	res << "series_id;size;epsilon_time;epsilon_distance;time_solve" << endl;

	float eptimes[] = {0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95};
	float epdists[] = {0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95};
	int sizes[] = { 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000 };

	for (const float &eptime : eptimes){
	    for (const float &epdist : epdists){
            for(const int &size : sizes){
                std::cout << "Loading TAP instance " << size << endl;
                std::stringstream fname;
                fname << "/users/21500078t/cplex_test/instances/tap_" << series_id << "_" << size << ".dat";
                const auto tap = Instance(fname.str());

                const auto solver = Solver(tap);

                int budget = lround(eptime * size * 27.5f);
                int dist_bound = lround( epdist * size * 4.5);

                double time = solver.solve_and_print(dist_bound, budget, false);
                std::cout << endl << "TIME TO SOLVE " << time << endl;
                res << series_id << ";" << eptime << ";" << epdist << ";" << size << ";" << time << endl;
                res.flush();
            }
        }
    }

	res.flush();
	res.close();
	return 0;
}

int main(int argc, char* argv[]) {
    //int tap_series = stoi(argv[1]);
	//return run_scale_test(tap_series);

	return run_debug();
}

