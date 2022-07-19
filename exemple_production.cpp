#include "instance.h"
#include "solver.h"
#include "SolverVPLSHammingSX.h"
#include <iostream>
#include <string>
#include <math.h>

int run_exact_test(char* argv[]) {

	using namespace cplex_tap;

	ofstream res;
	std::stringstream outname;
	outname << "/users/21500078t/res_colgen/" << argv[2];
	res.open(outname.str());
	res << "series_id;size;epsilon_time;epsilon_distance;time_solve;z;nb_nodes;solution" << endl;
    res.precision(17);

	float eptimes[] = {0.6};
	float epdists[] = {0.3};
	int sizes[] = {40, 60, 80, 100, 200, 300};

    for(const int &size : sizes){
	    for (const float &epdist : epdists){
            for (const float &eptime : eptimes){
                    for (int series_id = 0; series_id < 30; ++series_id) {
                        std::cout << "Loading TAP instance " << size << endl;
                        std::stringstream fname;
                        fname << argv[1] << "/tap_" << series_id << "_" << size << ".dat";
                        const auto tap = Instance(fname.str());

                        const auto solver = Solver(tap);

                        //int budget = lround(eptime * size * 27.5f);
                        //int budget = lround(eptime * size * 6.f);
                        int budget = lround(eptime * size * 27.5f); //f1
                        //int dist_bound = lround(epdist * size * 4.5);
                        //int dist_bound = lround(epdist * size * 7.f);
                        int dist_bound = lround(epdist * size * 5.5f); //f1

                        Solution sol = solver.solve(dist_bound, budget, false, false, false, "");
                        std::cout << endl << "TIME TO SOLVE " << sol.time << endl;
                        res << series_id << ";" << size << ";" << eptime << ";" << epdist << ";" << sol.time << ";"
                            << sol.z << ";" << sol.nodes << ";";
                        for (int i = 0; i < sol.sequence.size(); ++i) {
                            res << sol.sequence.at(i);
                            if (i != sol.sequence.size() - 1)
                                res << ",";
                        }
                        res << endl;
                        res.flush();

                    }
            }
        }
    }

	res.flush();
	res.close();
	return 0;
}

int production(char* argv[]) {
    using namespace cplex_tap;
    const auto tap = Instance(argv[3]);
    const auto solver = Solver(tap);

    //int budget = lround(stod(argv[1]) * tap.size() * 27.5f);
    int budget = lround(stod(argv[1]));
    //int dist_bound = lround( stod(argv[2]) * tap.size() * 4.5);
    int dist_bound = lround( stod(argv[2]));

    Solution sol = solver.solve(dist_bound, budget, false, false, false, "");

    return 0;
}

int run_debug() {
    std::string fname = "/users/21500078t/instances/tap_7_60.dat";

    using namespace cplex_tap;
    const auto tap = Instance(fname);

    const auto solver = Solver(tap);

    // Easy
    int budget = lround( 0.25 * tap.size() * 27.5f);
    int dist_bound = lround( 0.35 * tap.size() * 4.5);

    //std::cout << budget << "/" << dist_bound << std::endl;

    Solution sol = solver.solve(dist_bound, budget, false, false, false, "/home/alex/instances/tap_1_300.warm");
    std::cout << endl << "TIME TO SOLVE " << sol.time << endl;
    std::cout << sol.time << ";" << sol.z << ";"  << endl;


    return 0;
}

int main(int argc, char* argv[]) {
    std::cout.precision(17);
    //run_debug();
    //return experiments_vpls(argv);
    //return production(argv);
    run_exact_test(argv);
}

