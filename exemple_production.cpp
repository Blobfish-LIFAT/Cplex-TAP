#include "instance.h"
#include "solver.h"
#include "SolverVPLSHammingSX.h"
#include <iostream>
#include <string>
#include <math.h>
#include "CGTAPInstance.h"
#include "pricingSolver.h"
#include "KnapsackSolver.h"
#include "princingCPSolver.h"

#include "InitTargets.h"

inline constexpr auto hash_djb2a(const std::string_view sv) {
    unsigned long hash{ 5381 };
    for (unsigned char c : sv) {
        hash = ((hash << 5) + hash) ^ c;
    }
    return hash;
}

inline constexpr auto operator"" _sh(const char *str, size_t len) {
    return hash_djb2a(std::string_view{ str, len });
}


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
    cout << "OPT ? " << sol.optimal << endl << "Z=" << sol.z << endl;
    for (int i = 0; i < sol.sequence.size(); ++i) {
        cout << sol.sequence[i];
        if (i != sol.sequence.size() - 1)
            cout << " ";
    }
    cout << endl;
    return 0;
}

vector<cplex_tap::Query> getAllQueries(const cplex_tap::CGTAPInstance *ist);

void dump_instance(vector<cplex_tap::Query> queries, cplex_tap::CGTAPInstance ist, std::string path){
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
        out << q.getAgg() << ";";
        out << ist.getMeasureId(q.getMeasureLeft()) << ";" << ist.getMeasureId(q.getMeasureRight()) << ";";
        out << ist.getDimId(q.getGbAttribute()) << ";";
        for (auto pair : q.getLeftPredicate()) {
            out << pair.first << "=" << pair.second << "&";
        }
        out << ";";
        for (auto pair : q.getRightPredicate()) {
            out << pair.first << "=" << pair.second << "&";
        }

        out << endl;
    }

    out.close();

    std::cout << "Dump complete" << std::endl;
}

vector<cplex_tap::Query> getAllQueries(const cplex_tap::CGTAPInstance *ist) {
    vector<cplex_tap::Query> queries;
    for (int i = 0; i < ist->getNbDims(); ++i) {
        for (int j = 0; j < ist->getNbDims(); ++j) {
            if (i != j){
                for (int k = 0; k < ist->getAdSize(j); ++k) {
                    for (int l = k+1; l < ist->getAdSize(j); ++l) {
                        vector<pair<string, int> > lPred = {{ist->getDimName(j), l}};
                        vector<pair<string, int> > rPred = {{ist->getDimName(j), k}};
                        queries.emplace_back(cplex_tap::Query(ist->getTableName(), "sum", ist->getDimName(i), ist->getMeasureName(0), ist->getMeasureName(0), lPred, rPred));
                    }
                }
            }
        }
    }
    return queries;
}




int run_debug(char* argv[]) {
    using namespace cplex_tap;

    auto cgIST = cplex_tap::CGTAPInstance("/home/alex/tap_instances/demo_cg_6");

    std::cout<< "--- LOADING DONE ---" << std::endl;
    time_t start, end;
    start = clock();

    auto ini = KMppInit(cgIST, 500, 200);
    auto starting = ini.build(50);

    double time_to_init = (double)(clock() - start) / (double)CLOCKS_PER_SEC;
    std::cout << "Generated " << starting.size() << " queries" << std::endl;
    std::cout<< "--- INIT COMPLETE ["<< time_to_init <<"]---" << std::endl;

    start = ::clock();
    princingCPSolver solver = princingCPSolver(cgIST, 500, 200, starting);
    Solution s = solver.solve();
    std::cout<< "--- Solved ["<< (double)(clock() - start) / (double)CLOCKS_PER_SEC <<"]---" << std::endl;

    end = clock();
    double time_to_sol = (double)(end - start) / (double)CLOCKS_PER_SEC;
    cout << "[TIME] TOTAL " << time_to_sol << endl;

    return 0;
}

int main(int argc, char* argv[]) {
    std::cout.precision(17);
    //return run_debug(argv);

    //(ep_t,ep_d)
    int ep_t = stoi(argv[3]);
    int ep_d = stoi(argv[4]);
    std::string ist_path = argv[1];
    std::string init_profile = argv[2];

    int cplex_sym = stoi(argv[5]);
    long rng_seed = stol(argv[6]);

    auto cgIST = cplex_tap::CGTAPInstance(ist_path);
    vector<cplex_tap::Query> starting_queries;

    // Dump instance for debug
    dump_instance(getAllQueries(&cgIST), cgIST, "/home/alex/tap_ist_dump");

    time_t start, end;
    start = clock();

    switch (hash_djb2a(init_profile)) {
        case "rd_300"_sh:
            starting_queries = rd_300(cgIST, rng_seed);
            break;
        case "rd_50"_sh:
            starting_queries = rd_50(cgIST, rng_seed);
            break;
        case "rd_100"_sh:
            starting_queries = rd_100(cgIST, rng_seed);
            break;
        case "rd_150"_sh:
            starting_queries = rd_150(cgIST, rng_seed);
            break;
        case "rd_200"_sh:
            starting_queries = rd_200(cgIST, rng_seed);
            break;
        case "rd_div_50"_sh:
            starting_queries = rd_div_50(cgIST, rng_seed);
            break;
        case "rd_div_100"_sh:
            starting_queries = rd_div_100(cgIST, rng_seed);
            break;
        case "rd_div_150"_sh:
            starting_queries = rd_div_150(cgIST, rng_seed);
            break;
        case "rd_div_25_25"_sh:
            starting_queries = rd_div_25_25(cgIST, rng_seed);
            break;
        case "rd_div_50_50"_sh:
            starting_queries = rd_div_50_50(cgIST, rng_seed);
            break;
        case "rd_div_int_2_24_24"_sh:
            starting_queries = rd_div_int_2_24_24(cgIST, rng_seed);
            break;
        case "rd_div_int_2_49_49"_sh:
            starting_queries = rd_div_int_2_49_49(cgIST, rng_seed);
            break;
        case "rd_div_int_2_74_74"_sh:
            starting_queries = rd_div_int_2_74_74(cgIST, rng_seed);
            break;
        case "rd_div_int_2_99_99"_sh:
            starting_queries = rd_div_int_2_99_99(cgIST, rng_seed);
            break;
        case "rd_div_int_2_149_149"_sh:
            starting_queries = rd_div_int_2_149_149(cgIST, rng_seed);
            break;
        case "rd_div_int_17_17_17"_sh:
            starting_queries = rd_div_int_17_17_17(cgIST, rng_seed);
            break;
        case "rd_div_int_33_33_33"_sh:
            starting_queries = rd_div_int_33_33_33(cgIST, rng_seed);
            break;
        case "rd_div_int_50_50_50"_sh:
            starting_queries = rd_div_int_50_50_50(cgIST, rng_seed);
            break;
        case "rd_div_int_68_66_66"_sh:
            starting_queries = rd_div_int_68_66_66(cgIST, rng_seed);
            break;
        case "rd_div_int_100_100_100"_sh:
            starting_queries = rd_div_int_100_100_100(cgIST, rng_seed);
            break;
        case "rd_int_div_17_17_17"_sh:
            starting_queries = rd_int_div_17_17_17(cgIST, rng_seed);
            break;
        case "rd_int_div_33_33_33"_sh:
            starting_queries = rd_int_div_33_33_33(cgIST, rng_seed);
            break;
        case "rd_int_div_50_50_50"_sh:
            starting_queries = rd_int_div_50_50_50(cgIST, rng_seed);
            break;
        case "rd_int_div_68_66_66"_sh:
            starting_queries = rd_int_div_68_66_66(cgIST, rng_seed);
            break;
        case "rd_int_div_100_100_100"_sh:
            starting_queries = rd_int_div_100_100_100(cgIST, rng_seed);
            break;
        case "rd_int_div_2_24_24"_sh:
            starting_queries = rd_int_div_2_24_24(cgIST, rng_seed);
            break;
        case "rd_int_div_2_49_49"_sh:
            starting_queries = rd_int_div_2_49_49(cgIST, rng_seed);
            break;
        case "rd_int_div_2_74_74"_sh:
            starting_queries = rd_int_div_2_74_74(cgIST, rng_seed);
            break;
        case "rd_int_div_2_99_99"_sh:
            starting_queries = rd_int_div_2_99_99(cgIST, rng_seed);
            break;
        case "rd_int_div_2_149_149"_sh:
            starting_queries = rd_int_div_2_149_149(cgIST, rng_seed);
            break;
        case "rdsrt_50"_sh:
            starting_queries = rdsrt_50(cgIST, rng_seed);
            break;
        case "rdsrt_100"_sh:
            starting_queries = rdsrt_100(cgIST, rng_seed);
            break;
        case "rdsrt_150"_sh:
            starting_queries = rdsrt_150(cgIST, rng_seed);
            break;
        case "rdsrt_200"_sh:
            starting_queries = rdsrt_200(cgIST, rng_seed);
            break;
        case "rdsrt_300"_sh:
            starting_queries = rdsrt_300(cgIST, rng_seed);
            break;
        case "rdsk_50"_sh:
            starting_queries = rdsk_50(cgIST, rng_seed);
            break;
        case "rdsk_100"_sh:
            starting_queries = rdsk_100(cgIST, rng_seed);
            break;
        case "rdsk_150"_sh:
            starting_queries = rdsk_150(cgIST, rng_seed);
            break;
        case "rdsk_200"_sh:
            starting_queries = rdsk_200(cgIST, rng_seed);
            break;
        case "rdsk_300"_sh:
            starting_queries = rdsk_300(cgIST, rng_seed);
            break;
        case "rdks_50"_sh:
            starting_queries = rdks_50(cgIST, ep_t, ep_d, rng_seed);
            break;
        case "rdks_100"_sh:
            starting_queries = rdks_100(cgIST, ep_t, ep_d, rng_seed);
            break;
        case "rdks_200"_sh:
            starting_queries = rdks_200(cgIST, ep_t, ep_d, rng_seed);
            break;
        case "rdks_150"_sh:
            starting_queries = rdks_150(cgIST, ep_t, ep_d, rng_seed);
            break;
        case "rdks_300"_sh:
            starting_queries = rdks_300(cgIST, ep_t, ep_d, rng_seed);
            break;
        case "kmeanspp_50"_sh:
            starting_queries = kmeanspp_50(cgIST, ep_t, ep_d, rng_seed);
            break;
        case "kmeanspp_100"_sh:
            starting_queries = kmeanspp_100(cgIST, ep_t, ep_d, rng_seed);
            break;
        case "kmeanspp_150"_sh:
            starting_queries = kmeanspp_150(cgIST, ep_t, ep_d, rng_seed);
            break;
        case "kmeanspp_200"_sh:
            starting_queries = kmeanspp_200(cgIST, ep_t, ep_d, rng_seed);
            break;
        case "kmeanspp_300"_sh:
            starting_queries = kmeanspp_300(cgIST, ep_t, ep_d, rng_seed);
            break;
        case "ks"_sh:
            std::cout<< "--- INIT COMPLETE ["<< 0 <<"]---" << std::endl;
            start = clock();
            auto queries = getAllQueries(&cgIST);
            cplex_tap::KnapsackSolver solver = cplex_tap::KnapsackSolver(cgIST);
            cplex_tap::Solution s = solver.solve(getAllQueries(&cgIST), ep_t, ep_d);
            std::cout << "[STEP][END] - z*=" << s.z << std::endl;
            cout << "[TIME][ITER][s] " << s.time << endl;
            cout << "[TIME] TOTAL " << (double)(::clock() - start) / (double)CLOCKS_PER_SEC << endl;
            return 0;
    }

    end = clock();
    double time_to_init = (double)(end - start) / (double)CLOCKS_PER_SEC;
    std::cout << "Generated " << starting_queries.size() << " queries" << std::endl;
    std::cout<< "--- INIT COMPLETE ["<< time_to_init <<"]---" << std::endl;

    start = clock();

    cplex_tap::pricingSolver solver = cplex_tap::pricingSolver(cgIST, ep_d, ep_t, starting_queries);

    solver.setCplexSym(cplex_sym);

    cplex_tap::Solution s = solver.solve();

    end = clock();
    double time_to_sol = (double)(end - start) / (double)CLOCKS_PER_SEC;
    cout << "[TIME] TOTAL " << time_to_sol << endl;


}

