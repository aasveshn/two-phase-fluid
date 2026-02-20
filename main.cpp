#include <iostream>
#include <vector>
#include <algorithm>
#include <memory>
#include "Mesh/Mesh.h"
#include "Solver/Solver.h"




int main()
{
	Mesh mesh("input/mesh.txt");
	auto solver = std::make_unique<Solver>(mesh);
	mesh.SetInitialCondidions(solver->getPhases());
	solver -> Solve();
		
	return 0;
}
