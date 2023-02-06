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
#include "RandomInitSkew.h"
#include "KSInit.h"
#include "KMppInit.h"

std::vector<cplex_tap::Query> rd_xx(int xx, cplex_tap::CGTAPInstance ist, long seed){
    cplex_tap::RandomInit rdi = cplex_tap::RandomInit(ist, seed);
    return rdi.build(xx);
}

auto rd_1 = std::bind(rd_xx, 1, std::placeholders::_1, std::placeholders::_2);
auto rd_10 = std::bind(rd_xx, 10, std::placeholders::_1, std::placeholders::_2);
auto rd_50 = std::bind(rd_xx, 50, std::placeholders::_1, std::placeholders::_2);
auto rd_100 = std::bind(rd_xx, 100, std::placeholders::_1, std::placeholders::_2);
auto rd_150 = std::bind(rd_xx, 150, std::placeholders::_1, std::placeholders::_2);

std::vector<cplex_tap::Query> rd_div_xx(int xx, cplex_tap::CGTAPInstance ist, long seed){
    cplex_tap::DiversificationInit dvi = cplex_tap::DiversificationInit(ist, rd_1(ist, seed));
    return dvi.build(xx);
}

auto rd_div_50 = std::bind(rd_div_xx, 50, std::placeholders::_1, std::placeholders::_2);
auto rd_div_100 = std::bind(rd_div_xx, 100, std::placeholders::_1, std::placeholders::_2);
auto rd_div_150 = std::bind(rd_div_xx, 150, std::placeholders::_1, std::placeholders::_2);

std::vector<cplex_tap::Query> rd_div_xx_yy(int xx, int yy, cplex_tap::CGTAPInstance ist, long seed){
    cplex_tap::DiversificationInit dvi = cplex_tap::DiversificationInit(ist, rd_xx(xx, ist, seed));
    return dvi.build(xx + yy);
}

auto rd_div_25_25 = std::bind(rd_div_xx_yy, 25, 25, std::placeholders::_1, std::placeholders::_2);
auto rd_div_50_50 = std::bind(rd_div_xx_yy, 50, 50, std::placeholders::_1, std::placeholders::_2);

std::vector<cplex_tap::Query> rd_div_int_xx_yy_zz(int xx, int yy, int zz, cplex_tap::CGTAPInstance ist, long seed){
    cplex_tap::DiversificationInit dvi = cplex_tap::DiversificationInit(ist, rd_xx(xx, ist, seed));
    auto tmp = dvi.build(xx + yy);
    cplex_tap::IntensificationInit ini = cplex_tap::IntensificationInit(ist, tmp);
    return ini.build(xx + yy + zz);
}

auto rd_div_int_2_24_24 = std::bind(rd_div_int_xx_yy_zz, 2, 24, 24, std::placeholders::_1, std::placeholders::_2);
auto rd_div_int_2_49_49 = std::bind(rd_div_int_xx_yy_zz, 2, 49, 49, std::placeholders::_1, std::placeholders::_2);
auto rd_div_int_2_74_74 = std::bind(rd_div_int_xx_yy_zz, 2, 74, 74, std::placeholders::_1, std::placeholders::_2);
auto rd_div_int_17_17_17 = std::bind(rd_div_int_xx_yy_zz, 17, 17, 17, std::placeholders::_1, std::placeholders::_2);
auto rd_div_int_33_33_33 = std::bind(rd_div_int_xx_yy_zz, 33, 33, 33, std::placeholders::_1, std::placeholders::_2);
auto rd_div_int_50_50_50 = std::bind(rd_div_int_xx_yy_zz, 50, 50, 50, std::placeholders::_1, std::placeholders::_2);

std::vector<cplex_tap::Query> rd_int_div_xx_yy_zz(int xx, int yy, int zz, cplex_tap::CGTAPInstance ist, long seed){
    cplex_tap::IntensificationInit dvi = cplex_tap::IntensificationInit(ist, rd_xx(xx, ist, seed));
    auto tmp = dvi.build(xx + yy);
    cplex_tap::DiversificationInit ini = cplex_tap::DiversificationInit(ist, tmp);
    return ini.build(xx + yy + zz);
}

auto rd_int_div_17_17_17 = std::bind(rd_int_div_xx_yy_zz, 17, 17, 17, std::placeholders::_1, std::placeholders::_2);
auto rd_int_div_33_33_33 = std::bind(rd_int_div_xx_yy_zz, 33, 33, 33, std::placeholders::_1, std::placeholders::_2);
auto rd_int_div_50_50_50 = std::bind(rd_int_div_xx_yy_zz, 50, 50, 50, std::placeholders::_1, std::placeholders::_2);
auto rd_int_div_2_24_24 = std::bind(rd_int_div_xx_yy_zz, 2, 24, 24, std::placeholders::_1, std::placeholders::_2);
auto rd_int_div_2_49_49 = std::bind(rd_int_div_xx_yy_zz, 2, 49, 49, std::placeholders::_1, std::placeholders::_2);
auto rd_int_div_2_74_74 = std::bind(rd_int_div_xx_yy_zz, 2, 74, 74, std::placeholders::_1, std::placeholders::_2);

std::vector<cplex_tap::Query> rdstr_xx(int xx, cplex_tap::CGTAPInstance ist, long seed){
    cplex_tap::randThenSortInit rdi = cplex_tap::randThenSortInit(ist, seed);
    return rdi.build(xx);
}

auto rdsrt_50 = std::bind(rdstr_xx, 50, std::placeholders::_1, std::placeholders::_2);
auto rdsrt_100 = std::bind(rdstr_xx, 100, std::placeholders::_1, std::placeholders::_2);
auto rdsrt_150 = std::bind(rdstr_xx, 150, std::placeholders::_1, std::placeholders::_2);

std::vector<cplex_tap::Query> rdsk_xx(int xx, cplex_tap::CGTAPInstance ist, long seed){
    cplex_tap::RandomInitSkew rdi = cplex_tap::RandomInitSkew(ist, seed);
    return rdi.build(xx);
}

auto rdsk_50 = std::bind(rdsk_xx, 50, std::placeholders::_1, std::placeholders::_2);
auto rdsk_100 = std::bind(rdsk_xx, 100, std::placeholders::_1, std::placeholders::_2);
auto rdsk_150 = std::bind(rdsk_xx, 150, std::placeholders::_1, std::placeholders::_2);

std::vector<cplex_tap::Query> rdks_xx(int xx, cplex_tap::CGTAPInstance ist, int ept, int epd, long seed){
    cplex_tap::KSInit rdi = cplex_tap::KSInit(ist, epd, ept, seed);
    return rdi.build(xx);
}

auto rdks_50 = std::bind(rdks_xx, 50, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
auto rdks_100 = std::bind(rdks_xx, 100, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
auto rdks_150 = std::bind(rdks_xx, 150, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);

std::vector<cplex_tap::Query> kmeanspp_xx(int xx, cplex_tap::CGTAPInstance ist, int ept, int epd, long seed){
    cplex_tap::KMppInit rdi = cplex_tap::KMppInit(ist, epd, ept, seed);
    return rdi.build(xx);
}

auto kmeanspp_50 = std::bind(kmeanspp_xx, 50, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
auto kmeanspp_100 = std::bind(kmeanspp_xx, 100, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
auto kmeanspp_150 = std::bind(kmeanspp_xx, 150, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);

#endif //CPLEX_TEST_INITTARGETS_H
