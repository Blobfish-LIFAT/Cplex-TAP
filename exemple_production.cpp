#include "instance.h"
#include "solver.h"
#include "SolverVPLS.h"
#include "SolverVPLSDet.h"
#include "SolverVPLSHamming.h"
#include "SolverVPLSHammingSX.h"
#include <iostream>
#include <string>
#include <math.h>

int run_exact_test() {

	using namespace cplex_tap;

	ofstream res;
	std::stringstream outname;
	outname << "/users/21500078t/res_cplex_exact_vtkpaper_600.csv";
	res.open(outname.str());
	res << "series_id;size;epsilon_time;epsilon_distance;time_solve;z;solution" << endl;
    res.precision(17);

	float eptimes[] = {0.25};
	float epdists[] = {0.35};
	int sizes[] = {20, 40, 60, 80, 100, 150, 200, 250, 300, 350, 400, 450, 500};

    for(const int &size : sizes){
	    for (const float &epdist : epdists){
            for (const float &eptime : eptimes){
                    for (int series_id = 0; series_id < 30; ++series_id) {
                        std::cout << "Loading TAP instance " << size << endl;
                        std::stringstream fname;
                        fname << "./instances/tap_" << series_id << "_" << size << ".dat";
                        const auto tap = Instance(fname.str());

                        const auto solver = Solver(tap);

                        int budget = lround(eptime * size * 27.5f);
                        int dist_bound = lround(epdist * size * 4.5);

                        Solution sol = solver.solve_and_print(dist_bound, budget, false, false, false, false, "");
                        std::cout << endl << "TIME TO SOLVE " << sol.time << endl;
                        res << series_id << ";" << size << ";" << eptime << ";" << epdist << ";" << sol.time << ";" << sol.z << ";" ;
                        for (int i = 0; i < sol.sequence.size(); ++i) {
                            res << sol.sequence.at(i);
                            if (i != sol.sequence.size()-1)
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
    ofstream csv;
    csv.open (argv[2]);
    if (stoi(argv[1]) == 0) csv << "instance;size;type;edist;etime;solve_time;z;solution" << std::endl;
    csv.precision(17);
    int sizes[] = {20, 40, 60, 80, 100, 150, 200, 250, 300, 350, 400, 450, 500};
    for (int size : sizes){
        if (stoi(argv[1]) > size)
            continue;
        for (int i = 0; i < 30; ++i) {
            if (size == 40 && (i == 9 || i == 13 || i == 19 || i == 21)) continue;
            std::string fname = "./instances/tap_" + std::to_string(i) + "_" + std::to_string(size) + ".dat";
            std::string warmname = "./instances/tap_" + std::to_string(i) + "_" + std::to_string(size) + ".warm";
            const auto tap = Instance(fname);
            int budget = lround( 0.25 * tap.size() * 27.5f);
            int dist_bound = lround( 0.35 * tap.size() * 4.5);
            const auto solverRandom = SolverVPLS(tap, 7, 15, 20, 90);
            const auto solverDet = SolverVPLSDet(tap, 5, 20, 20, 120);
            const auto solverHammingS = SolverVPLSHamming(tap, 7, 15, 20, 90);// 10 60
            const auto solverRandomSX = SolverVPLSHammingSX(tap, 5, 50, 20, 120);// 30 90
            Solution random = solverRandom.solve_and_print(dist_bound, budget, false, false, false, true, warmname);
            Solution det = solverDet.solve_and_print(dist_bound, budget, false, false, false, true, warmname);
            Solution hammingS = solverHammingS.solve_and_print(dist_bound, budget, false, false, false, true, warmname);
            Solution hammingSX = solverRandomSX.solve_and_print(dist_bound, budget, false, false, false, true, warmname);
            csv << "tap_" << i << "_" << size <<".dat;" << size << ";" << "random;" << 0.35 << ";" << 0.25 << ";" << random.time << ";" << random.z << ";" << std::endl;
            csv << "tap_" << i << "_" << size <<".dat;" << size << ";" << "det;" << 0.35 << ";" << 0.25 << ";" << det.time << ";" << det.z << ";" << std::endl;
            csv << "tap_" << i << "_" << size <<".dat;" << size << ";" << "s;" << 0.35 << ";" << 0.25 << ";" << hammingS.time << ";" << hammingS.z << ";" << std::endl;
            csv << "tap_" << i << "_" << size <<".dat;" << size << ";" << "sx;" << 0.35 << ";" << 0.25 << ";" << hammingSX.time << ";" << hammingSX.z << ";" << std::endl;
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
    std::string fname = "/home/alex/instances/tap_1_300.dat";

    using namespace cplex_tap;
    const auto tap = Instance(fname);

    const auto solver = SolverVPLSHammingSX(tap, 5, 50, 20, 120);

    // Easy
    int budget = lround( 0.25 * tap.size() * 27.5f);
    int dist_bound = lround( 0.35 * tap.size() * 4.5);

    Solution sol = solver.solve_and_print(dist_bound, budget, false, false, false, true, "/home/alex/instances/tap_1_300.warm");
    std::cout << endl << "TIME TO SOLVE " << sol.time << endl;
    std::cout << sol.time << ";" << sol.z << ";"  << endl;


    return 0;
}

int main(int argc, char* argv[]) {
    std::cout.precision(17);
    //run_debug();
    return experiments_vpls(argv);
    //return production(argv);
    //run_exact_test();
}

