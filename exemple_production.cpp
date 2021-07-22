#include "instance.h"
#include "solver.h"
#include "SolverVPLS.h"
#include "SolverVPLSDet.h"
#include "SolverVPLSHamming.h"
#include "SolverVPLSHammingSX.h"
#include <iostream>
#include <string>
#include <math.h>

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

                Solution sol = solver.solve_and_print(dist_bound, budget, false, false, false, false, "");
                std::cout << endl << "TIME TO SOLVE " << sol.time << endl;
                res << series_id << ";" << size << ";" << eptime << ";" << epdist << ";" << sol.time << endl;
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
    const auto tap = Instance(argv[3]);
    const auto solver = Solver(tap);

    int budget = lround(stod(argv[1]) * tap.size() * 27.5f);
    int dist_bound = lround( stod(argv[2]) * tap.size() * 4.5);

    Solution sol = solver.solve_and_print(dist_bound, budget, false, false, false, false, "");

    return 0;
}

int run_debug(bool progressive, double temps, double dist, std::string path) {
    using namespace cplex_tap;
    const auto tap = Instance(path);

    const auto solver = Solver(tap);

    // Easy
    int budget = lround(temps * tap.size() * 27.5f);
    int dist_bound = lround( dist * tap.size() * 4.5);

    Solution sol = solver.solve_and_print(dist_bound, budget, progressive, false, false, false, "");
    std::cout << endl << "TIME TO SOLVE " << sol.time << endl;


    return 0;
}

//bin tbound dbound h initTime epochTime instanceFile warmFile
int run_debug_vpls(char* argv[]) {
    using namespace cplex_tap;
    const auto tap = Instance(argv[6]);

    //const auto solver = SolverVPLSHammingSX(tap, 25, stoi(argv[3]), stoi(argv[4]), stoi(argv[5]));
    const auto solver = SolverVPLSDet(tap, 25, stoi(argv[3]), stoi(argv[4]), stoi(argv[5]));

    int budget = lround( stod(argv[1]) * tap.size() * 27.5f);
    int dist_bound = lround( stod(argv[2]) * tap.size() * 4.5);

    bool seed = true;
    if (string(argv[7]) == "none"){
        seed = false;
    }
    Solution sol = solver.solve_and_print(dist_bound, budget, false, false, false, seed, argv[7]);
    std::cout << endl << "TIME TO SOLVE " << sol.time << endl;


    return 0;
}

int experiments_vpls(char* argv[]){
    using namespace cplex_tap;
    ofstream csv;
    csv.open (argv[2]);
    std::cout << "instance;size;type;edist;etime;solve_time;z" << std::endl;
    int sizes[] = {20, 40, 60, 80, 100, 150, 200, 250, 300, 350, 400, 450, 500};
    for (int size : sizes){
        if (stoi(argv[1]) > size)
            continue;
        for (int i = 0; i < 30; ++i) {
            std::string fname = "./instances/tap_" + std::to_string(i) + "_" + std::to_string(size) + ".dat";
            std::string warmname = "./instances/tap_" + std::to_string(i) + "_" + std::to_string(size) + ".warm";
            const auto tap = Instance(fname);
            int budget = lround( 0.25 * tap.size() * 27.5f);
            int dist_bound = lround( 0.35 * tap.size() * 4.5);
            const auto solverRandom = SolverVPLS(tap, 25, 15, 20, 90);
            const auto solverDet = SolverVPLSDet(tap, 25, 20, 20, 120);
            const auto solverHammingS = SolverVPLSHamming(tap, 25, 10, 20, 60);
            const auto solverRandomSX = SolverVPLSHammingSX(tap, 25, 30, 20, 90);
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

int main(int argc, char* argv[]) {
    // uncomment for test campaign binary
    return experiments_vpls(argv);

	// uncomment for VPLS binary
    if (argc < 7){
        cout << "USE: ./binary tbound dbound h initTime epochTime instanceFile warmFile" << endl;
        exit(2);
    } else {
        std::cout << "PARAMS: " << argv[1] << " " << argv[2] << " " << argv[3] << " " << argv[4] << " " << argv[5] << endl;
    }
	return run_debug_vpls(argv);
}

