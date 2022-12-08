#include "instance.h"
#include "solver.h"
#include "SolverVPLSHammingSX.h"
#include <iostream>
#include <string>
#include <math.h>
#include "CGTAPInstance.h"
#include "pricingSolver.h"
#include "initSolver.h"

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

                        int budget = lround(eptime * size * 27.5f);
                        int dist_bound = lround(epdist * size * 4.5f);

                        Solution sol = solver.solve(dist_bound, budget, false, "");
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

    Solution sol = solver.solve(dist_bound, budget, false, "");

    return 0;
}

int run_debug(char* argv[]) {
    using namespace cplex_tap;

    const auto tap = Instance("/home/alex/instances/tap_1_200.dat");
    const auto solver = Solver(tap);

    int budget = lround( 0.6 * tap.size() * 27.5f);
    int dist_bound = lround( 0.3 * tap.size() * 4.5);


    Solution sol = solver.solve(dist_bound, budget, false, "");
    std::cout << endl << "TIME TO SOLVE " << sol.time << endl;
    std::cout << sol.time << ";" << sol.z << ";"  << endl;

    return 0;
}

void dump_instance(cplex_tap::CGTAPInstance ist, std::string path){
    // Make the queries
    std::vector<cplex_tap::Query> queries;
    for (int i = 0; i < ist.getNbDims(); ++i) {
        for (int j = 0; j < ist.getNbDims(); ++j) {
            if (i != j){
                for (int k = 0; k < ist.getAdSize(j); ++k) {
                    for (int l = k+1; l < ist.getAdSize(j); ++l) {
                        std::vector<std::pair<string, int> > lPred = {{ist.getDimName(j), l}};
                        std::vector<std::pair<string, int> > rPred = {{ist.getDimName(j), k}};
                        queries.emplace_back(cplex_tap::Query(ist.getTableName(), "sum", ist.getDimName(i), ist.getMeasureName(0), ist.getMeasureName(0), lPred, rPred));
                    }
                }
            }
        }
    }

    //Dump
    std::fstream out(path, std::ios::out);
    out << queries.size() << endl;
    auto interest = JVMAdapter::getInterest(queries, ist);
    for (int i = 0; i < queries.size(); ++i) {
        out << interest[i];
        if (i != interest.size() - 1)
            out << " ";
    }
    out << endl;
    auto time = JVMAdapter::getTime(queries, ist);
    for (int i = 0; i < queries.size(); ++i) {
        out << time[i];
        if (i != time.size() - 1)
            out << " ";
    }
    out << endl;
    for (int i = 0; i < queries.size(); ++i) {
        auto q = queries[i];
        out << q.getAgg() << " ";
        out << ist.getMeasureId(q.getMeasureLeft()) << " " << ist.getMeasureId(q.getMeasureRight()) << " ";
        out << ist.getDimId(q.getGbAttribute()) << " ";
        for (auto pair : q.getLeftPredicate()) {
            out << pair.first << "=" << pair.second << "&";
        }
        out << " ";
        for (auto pair : q.getRightPredicate()) {
            out << pair.first << "=" << pair.second << "&";
        }

        out << endl;
    }

    out.close();
}

int main(int argc, char* argv[]) {
    std::cout.precision(17);

    std::string demo = argv[1]; //

    auto cgIST = cplex_tap::CGTAPInstance(demo);
    //dump_instance(cgIST, "/home/alex/CLionProjects/Cplex-TAP/ist_dump.dat");

    //cplex_tap::initSolver init(cgIST, 50);
    //auto startingSet = init.build_starting_set();

    std::cout<< "done.." << std::endl;
    time_t start, end;
    start = clock();
    cplex_tap::pricingSolver solver = cplex_tap::pricingSolver(cgIST, 125, 1000, stoi(argv[2]));
    cplex_tap::Solution s = solver.solve();
    end = clock();
    double time_to_sol = (double)(end - start) / (double)CLOCKS_PER_SEC;
    cout << "[TIME] TOTAL " << time_to_sol << endl;



    //run_debug(argv);
    //return production(argv);
    //run_exact_test(argv);
}

