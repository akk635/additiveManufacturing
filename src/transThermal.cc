#include "thermalSimulation.hpp"

#include <iostream>
#include <fstream>
#include <typeinfo>

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
	thermal.assignLayerIterators();
	thermal.defaultSetup(); // to set the intial fe for the first layer
	int layer = 3;

	for (int i = 0; i < layer; i++){
		thermal.setup_system();
		if (i < layer-1)
		thermal.pretimeStepping();
		//thermal.postCoolingStep();
		const std::string filename = "solution-"+ Utilities::int_to_string(i, 3)+ ".vtk";
		thermal.outputResults(filename);
		if(i<layer-1)
			thermal.update_active_fe_indices();
	}
	std::cout << "looping done" << std::endl;
}
