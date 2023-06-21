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
#include <sstream>
#include <regex>
#include "solver.h"
#include "InitTargets.h"

static cplex_tap::Instance buildRMPInstance(vector<cplex_tap::Query>& queries, cplex_tap::CGTAPInstance pricingIST) {
    vector<double> interest = JVMAdapter::getInterest(queries, pricingIST);
    vector<double> time = JVMAdapter::getTime(queries, pricingIST);
    vector<vector<int>> distMatrix;
    for (int i = 0; i < queries.size(); ++i) {
        vector<int> line;
        for (int j = 0; j < queries.size(); ++j) {
            line.emplace_back(queries[i].dist(queries[j]));
        }
        distMatrix.emplace_back(line);
    }
    return {static_cast<int>(queries.size()), interest, time, distMatrix};
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
    int n = 0;
    for (int i = 0; i < ist->getNbDims(); ++i) {
        for (int j = 0; j < ist->getNbDims(); ++j) {
            if (i != j){
                for (int k = 0; k < ist->getAdSize(j); ++k) {
                    for (int l = k+1; l < ist->getAdSize(j); ++l) {
                        vector<pair<string, int> > lPred = {{ist->getDimName(j), l}};
                        vector<pair<string, int> > rPred = {{ist->getDimName(j), k}};
                        n++;
                        //queries.emplace_back(cplex_tap::Query(ist->getTableName(), "sum", ist->getDimName(i), ist->getMeasureName(0), ist->getMeasureName(0), lPred, rPred));
                    }
                }
            }
        }
    }
    cout << "[INFO] Generated all queries " << n++ << endl;
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

    string solver_conf = argv[5];
    long rng_seed = stol(argv[6]);

    auto cgIST = cplex_tap::CGTAPInstance(ist_path);
    vector<cplex_tap::Query> starting_queries;

    // Dump instance for debug
    //dump_instance(getAllQueries(&cgIST), cgIST, "/home/alex/tap_ist_dump");

    //getAllQueries(&cgIST);
    //exit(0);

    time_t start, end;
    start = clock();

    // Parse and select init method
    std::stringstream ini(init_profile);
    std::string segment;
    std::vector<std::string> seglist;
    while(std::getline(ini, segment, '_'))
        seglist.push_back(segment);

    // Simple regular expression matching
    const std::string fnames[] = {"foo.txt", "bar.txt", "baz.dat", "zoidberg"};
    const std::regex rd_regex("rd_[0-9]+");
    const std::regex rd_div_regex("rd_div_[0-9]+");
    const std::regex rd_int_regex("rd_int_[0-9]+");
    const std::regex rd_div_int_regex("rd_div_int_[0-9]+");
    const std::regex rd_int_div_regex("rd_int_div_[0-9]+");
    const std::regex rdsrt_regex("rdsrt_[0-9]+");
    const std::regex rdks_regex("rdks_[0-9]+");
    const std::regex kmeanspp_regex("kmeanspp_[0-9]+");

    if (std::regex_match(init_profile, rd_regex)){
        int xx = stoi(seglist[1]);
        starting_queries = rd_xx(xx, cgIST, rng_seed);
    }
    else if (std::regex_match(init_profile, rd_div_regex)){
        int xx = stoi(seglist[2]);
        starting_queries = rd_div_xx(xx, cgIST, rng_seed);
    }
    else if (std::regex_match(init_profile, rd_int_regex)){
        int xx = stoi(seglist[2]);
        starting_queries = rd_int_xx(xx, cgIST, rng_seed);
    }
    else if (std::regex_match(init_profile, rd_div_int_regex)){
        int uval = stoi(seglist[3]);
        int split = uval / 2;
        starting_queries = rd_div_int_xx_yy_zz(1, split -1, split, cgIST, rng_seed);
    }
    else if (std::regex_match(init_profile, rd_int_div_regex)){
        int uval = stoi(seglist[3]);
        int split = uval / 2;
        starting_queries = rd_int_div_xx_yy_zz(1, split -1, split, cgIST, rng_seed);
    }
    else if (std::regex_match(init_profile, rdsrt_regex)){
        int xx = stoi(seglist[1]);
        starting_queries = rdstr_xx(xx, cgIST, rng_seed);
    }
    else if (std::regex_match(init_profile, rdks_regex)){
        int xx = stoi(seglist[1]);
        starting_queries = rdks_xx(xx, cgIST, ep_t, ep_d, rng_seed);
    }
    else if (std::regex_match(init_profile, kmeanspp_regex)){
        int xx = stoi(seglist[1]);
        starting_queries = kmeanspp_xx(xx, cgIST, ep_t, ep_d, rng_seed);
    }
    else {
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

    // LP
    start = clock();

    cplex_tap::pricingSolver solver = cplex_tap::pricingSolver(cgIST, ep_d, ep_t, starting_queries);
    solver.setCplexSym(0);
    solver.setPricingItTimeout(12);
    solver.setGlobalTimeout(1600);
    solver.setSelectedConf("best");
    cplex_tap::Solution s = solver.solve();

    end = clock();
    double time_to_sol = (double)(end - start) / (double)CLOCKS_PER_SEC;
    cout << "[TIME] TOTAL " << time_to_sol << endl;

    //MIP
    start = clock();

    cplex_tap::pricingSolver solver_mip  = cplex_tap::pricingSolver(cgIST, ep_d, ep_t, starting_queries);
    solver_mip.setCplexSym(0);
    solver_mip.setUseFloat(false);
    solver_mip.setPricingItTimeout(600);
    solver_mip.setGlobalTimeout(1600);
    solver_mip.setSelectedConf("default");
    s = solver_mip.solve();

    end = clock();
    time_to_sol = (double)(end - start) / (double)CLOCKS_PER_SEC;
    cout << "[TIME] TOTAL " << time_to_sol << endl;

    /*
     * Test staring pools
     *
    auto solver = cplex_tap::KnapsackSolver(cgIST);
    cplex_tap::Solution s = solver.solve(starting_queries, ep_t, ep_d);
    std::cout << "[POOL][ks] - z*=" << s.z << std::endl;

    const cplex_tap::Instance &rmpInstance = buildRMPInstance(starting_queries, cgIST);
    if (starting_queries.size() < 1001) {
        auto solver_math = cplex_tap::SolverVPLSHammingSX(rmpInstance, 15, 15, 30, 20);
        s = solver_math.solve(ep_d, ep_t, false, "");
        std::cout << "[POOL][math] - z*=" << s.z << std::endl;
    }

    if (starting_queries.size() < 501){
        auto solver_exact = cplex_tap::Solver(rmpInstance);
        s = solver_exact.solve(ep_d, ep_t, false, "");
        std::cout << "[POOL][cplex] - z*=" << s.z << std::endl;
    }*/

}



