#include "thermalSimulation.hpp"

#include <iostream>
#include <fstream>
#include <typeinfo>

int main(){
	std::cout<< "Additive manufacturing start" << std::endl;
	const int dim = 3;
	Triangulation<dim> triangulation;
	transThermal<dim> thermal(triangulation);
	thermal.readFullMesh();
	thermal.assignLayerIterators();
	thermal.defaultSetup();
	// int layers = thermal.nZ;
	int layers = 3;
	for (int i = 1; i < layers; i++){
		thermal.update_active_fe_indices();
		thermal.setup_system();

		thermal.pretimeStepping();
//		thermal.postCoolingStep();
		thermal.solveTag = true;
		const std::string filename = "solution-"+ Utilities::int_to_string(i, 3)+ ".vtk";
		thermal.outputResults(filename);
	}
	std::cout << "looping done" << std::endl;
}
