#include "thermalSimulation.hpp"

#include <iostream>
#include <fstream>

/*template <int dim>
void transThermal<dim>::transMesh(float time){
	typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.endc();
}*/

int main(){
	std::cout<< "Additive manufacturing start" << std::endl;
	const int dim = 3;
	Triangulation<dim> triangulation;
	transThermal<dim> thermal(triangulation);
	thermal.readFullMesh();
	thermal.defaultSetup();
	thermal.update_active_fe_indices();
	thermal.setup_system();
	thermal.assemble_system();
}
