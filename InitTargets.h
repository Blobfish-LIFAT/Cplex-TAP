//
// Created by alex on 06/02/23.
//

#ifndef CPLEX_TEST_INITTARGETS_H
#define CPLEX_TEST_INITTARGETS_H
/*
 *    --- Init modes ---
 */
#include "RandomInit.h"
#include "IntensificationInit.h"
#include "DiversificationInit.h"
#include "randThenSortInit.h"
//#include "RandomInitSkew.h"
#include "KSInit.h"
#include "KMppInit.h"

std::vector<cplex_tap::Query> rd_xx(int xx, cplex_tap::CGTAPInstance ist, long seed){
    cplex_tap::RandomInit rdi = cplex_tap::RandomInit(ist, seed);
    return rdi.build(xx);
}

auto rd_1 = std::bind(rd_xx, 1, std::placeholders::_1, std::placeholders::_2);

std::vector<cplex_tap::Query> rd_div_xx(int xx, cplex_tap::CGTAPInstance ist, long seed){
    cplex_tap::DiversificationInit dvi = cplex_tap::DiversificationInit(ist, rd_1(ist, seed));
    return dvi.build(xx);
}


std::vector<cplex_tap::Query> rd_int_xx(int xx, cplex_tap::CGTAPInstance ist, long seed){
    cplex_tap::IntensificationInit dvi = cplex_tap::IntensificationInit(ist, rd_1(ist, seed));
    return dvi.build(xx);
}


std::vector<cplex_tap::Query> rd_div_int_xx_yy_zz(int xx, int yy, int zz, cplex_tap::CGTAPInstance ist, long seed){
    cplex_tap::DiversificationInit dvi = cplex_tap::DiversificationInit(ist, rd_xx(xx, ist, seed));
    auto tmp = dvi.build(xx + yy);
    cplex_tap::IntensificationInit ini = cplex_tap::IntensificationInit(ist, tmp);
    return ini.build(xx + yy + zz);
}

std::vector<cplex_tap::Query> rd_int_div_xx_yy_zz(int xx, int yy, int zz, cplex_tap::CGTAPInstance ist, long seed){
    cplex_tap::IntensificationInit dvi = cplex_tap::IntensificationInit(ist, rd_xx(xx, ist, seed));
    auto tmp = dvi.build(xx + yy);
    cplex_tap::DiversificationInit ini = cplex_tap::DiversificationInit(ist, tmp);
    return ini.build(xx + yy + zz);
}


std::vector<cplex_tap::Query> rdstr_xx(int xx, cplex_tap::CGTAPInstance ist, long seed){
    cplex_tap::randThenSortInit rdi = cplex_tap::randThenSortInit(ist, seed);
    return rdi.build(xx);
}


std::vector<cplex_tap::Query> rdks_xx(int xx, cplex_tap::CGTAPInstance ist, int ept, int epd, long seed){
    cplex_tap::KSInit rdi = cplex_tap::KSInit(ist, epd, ept, seed);
    return rdi.build(xx);
}

std::vector<cplex_tap::Query> kmeanspp_xx(int xx, cplex_tap::CGTAPInstance ist, int ept, int epd, long seed){
    cplex_tap::KMppInit rdi = cplex_tap::KMppInit(ist, epd, ept, seed);
    return rdi.build(xx);
}



#endif //CPLEX_TEST_INITTARGETS_H
