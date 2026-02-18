#include <iostream>
#include <vector>
#include <algorithm>
#include <memory>
#include "Mesh/Mesh.h"
#include "Solver/Solver.h"




int main()
{
	Mesh mesh("input/mesh.txt");
	mesh.SetInitialCondidions();
	//Solver solver(mesh);
	//solver.Solve();
	auto solver = std::make_unique<Solver>(mesh);
	solver -> Solve();
		
	return 0;
}
