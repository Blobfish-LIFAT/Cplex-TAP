#include "instance.h"
#include "solver.h"
#include <iostream>
#include <string>

int main() {
	using namespace cplex_tap;

	std::cout << "Loading TAP instance\n";
	const auto tap = Instance("C:\\Users\\achan\\source\\repos\\cplex_test\\small_test_instance.txt");
	const auto solver = Solver(tap);

	solver.solve_and_print();

	return 0;
}