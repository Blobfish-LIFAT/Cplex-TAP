#include "instance.h"
#include "solver.h"
#include "SolverVPLS.h"
#include "SolverVPLSDet.h"
#include "SolverVPLSHamming.h"
#include "SolverVPLSHammingSX.h"
#include <iostream>
#include <string>
#include <math.h>

int run_exact_test(char* argv[]) {

	using namespace cplex_tap;

	ofstream res;
	std::stringstream outname;
	outname << "/users/21500078t/res2022/" << argv[2];
	res.open(outname.str());
	res << "series_id;size;epsilon_time;epsilon_distance;time_solve;z;nb_nodes;ss;solution" << endl;
    res.precision(17);

	float eptimes[] = {0.6};
	float epdists[] = {0.3};
	int sizes[] = {400, 500, 600, 700};//40, 60, 80, 100, 200, 300

    for(const int &size : sizes){
	    for (const float &epdist : epdists){
            for (const float &eptime : eptimes){
                    for (int series_id = 0; series_id < 30; ++series_id) {
                        int ss = -1;//lround(0.8*size);

                        std::cout << "Loading TAP instance " << size << endl;
                        std::stringstream fname;
                        fname << argv[1] << "/tap_" << series_id << "_" << size << ".dat";
                        const auto tap = Instance(fname.str());

                        const auto solver = Solver(tap, ss);

                        //int budget = lround(eptime * size * 27.5f);
                        //int budget = lround(eptime * size * 6.f);
                        int budget = lround(eptime * size * 27.5f); //f1
                        //int dist_bound = lround(epdist * size * 4.5);
                        //int dist_bound = lround(epdist * size * 7.f);
                        int dist_bound = lround(epdist * size * 5.5f); //f1

                        Solution sol = solver.solve_and_print(dist_bound, budget, false, false, false, false, "");
                        std::cout << endl << "TIME TO SOLVE " << sol.time << endl;
                        res << series_id << ";" << size << ";" << eptime << ";" << epdist << ";" << sol.time << ";"
                            << sol.z << ";" << sol.nodes << ";" << ss << ";";
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

int experiments_vpls(char* argv[]){
    using namespace cplex_tap;
    // ep
    double eptime = 0.6;
    double epdist = 0.3;

    // result file
    ofstream csv;
    std::stringstream outname;
    outname << "/users/21500078t/res2022/" << argv[1];
    csv.open(outname.str());
    csv << "instance;size;type;edist;etime;solve_time;z;solution" << std::endl;
    csv.precision(17); // avoid systematic rounding

    int sizes[] = {400, 500, 600, 700};//40, 60, 80, 100,  200,  300,
    for (int size : sizes){
        for (int i = 0; i < 30; ++i) {
            if (size == 40 && (i == 9 || i == 13 || i == 19 || i == 21)) continue;
            std::string fname = "/users/21500078t/instances/filtered_instances/tap_" + std::to_string(i) + "_" + std::to_string(size) + ".dat";
            std::string warmname = "/users/21500078t/instances/filtered_instances/tap_" + std::to_string(i) + "_" + std::to_string(size) + ".warm";
            const auto tap = Instance(fname);
            int budget = lround( eptime * size * 27.5f);
            int dist_bound = lround( epdist * size * 4.5);
            const auto solverRandom = SolverVPLS(tap, 7, 15, 20, 90);
            const auto solverDet = SolverVPLSDet(tap, 5, 20, 20, 120);
            const auto solverHammingS = SolverVPLSHamming(tap, 7, 15, 20, 90);// 10 60
            const auto solverRandomSX = SolverVPLSHammingSX(tap, 5, 50, 20, 120);// 30 90
            Solution random = solverRandom.solve_and_print(dist_bound, budget, false, false, false, true, warmname);
            Solution det = solverDet.solve_and_print(dist_bound, budget, false, false, false, true, warmname);
            Solution hammingS = solverHammingS.solve_and_print(dist_bound, budget, false, false, false, true, warmname);
            Solution hammingSX = solverRandomSX.solve_and_print(dist_bound, budget, false, false, false, true, warmname);
            csv << "tap_" << i << "_" << size <<".dat;" << size << ";" << "random;" << epdist << ";" << eptime << ";" << random.time << ";" << random.z << ";" << std::endl;
            csv << "tap_" << i << "_" << size <<".dat;" << size << ";" << "det;" << epdist << ";" << eptime << ";" << det.time << ";" << det.z << ";" << std::endl;
            csv << "tap_" << i << "_" << size <<".dat;" << size << ";" << "s;" << epdist << ";" << eptime << ";" << hammingS.time << ";" << hammingS.z << ";" << std::endl;
            csv << "tap_" << i << "_" << size <<".dat;" << size << ";" << "sx;" << epdist << ";" << eptime << ";" << hammingSX.time << ";" << hammingSX.z << ";" << std::endl;
            csv.flush();
        }
    }
    csv.close();
}

int production(char* argv[]) {
    using namespace cplex_tap;
    const auto tap = Instance(argv[3]);
    const auto solver = Solver(tap);

    //int budget = lround(stod(argv[1]) * tap.size() * 27.5f);
    int budget = lround(stod(argv[1]));
    //int dist_bound = lround( stod(argv[2]) * tap.size() * 4.5);
    int dist_bound = lround( stod(argv[2]));

    Solution sol = solver.solve_and_print(dist_bound, budget, false, false, false, false, "");

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

    Solution sol = solver.solve_and_print(dist_bound, budget, false, false, false, false, "/home/alex/instances/tap_1_300.warm");
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

