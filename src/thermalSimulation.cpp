/*
 * thermalSimulation.cpp
 *
 *  Created on: Mar 18, 2015
 *      Author: karthik
 */
#include "thermalSimulation.hpp"

template <int dim>
transThermal<dim>::transThermal(Triangulation<dim> & tria):triangulation(tria), dof_handler(tria),activeFE(1){
	fe_collection.push_back(zeroFE);
	fe_collection.push_back(activeFE);
}

// To read the external mesh from salome platform
template <int dim>
void transThermal<dim>::readFullMesh(){
	GridIn<dim> grid_in;
	grid_in.attach_triangulation (triangulation);
	std::ifstream input_file("../MeshFiles/60Structure.unv");
	grid_in.read_unv (input_file);

	std::cout << "   Number of active cells: "
			<< triangulation.n_active_cells()
			<< std::endl
			<< "   Total number of cells: "
			<< triangulation.n_cells()
			<< std::endl;
}

template <int dim>
void transThermal<dim>::defaultSetup(){
	// Initial additive manufacturing layers
	currentLayerCell = dof_handler.begin_active();
	topLayerCell = currentLayerCell->neighbor(dim-1);

	std::cout << "Height layer 1" << "\t" << currentLayerCell->center()[dim-1] << std::endl;
	std::cout << "Height layer 2" << "\t" << topLayerCell->center()[dim-1] << std::endl;


}

template <int dim>
void transThermal<dim>::update_active_fe_indices(){
	int i = 0;
	for (typename hp::DoFHandler<dim>::active_cell_iterator
		       cell = dof_handler.begin();
		       cell != dof_handler.end(); ++cell){
		if(cell->center()[dim-1] <= currentLayerCell->center()[dim-1]){
			cell->set_active_fe_index (1);
			i++;
		}else{
			cell->set_active_fe_index (0);
		}
	}
	std::cout << "#layerCells" << "\t" << i << std::endl;
}

template <int dim>
void transThermal<dim>::incrementLayer(){
	// at the end of each time step
	if (currentLayerCell->face(dim-1)->at_boundary())
		return;

	currentLayerCell = topLayerCell;
	topLayerCell = currentLayerCell->neighbor(dim-1);

}

template <int dim>
void transThermal<dim>::setup_system(){
	dof_handler.distribute_dofs (fe_collection);
	CompressedSparsityPattern compressed_sparsity_pattern(dof_handler.n_dofs(), dof_handler.n_dofs());
	DoFTools::make_sparsity_pattern (dof_handler, compressed_sparsity_pattern);
	sparsity_pattern.copy_from (compressed_sparsity_pattern);

	system_matrix.reinit (sparsity_pattern);
	solution.reinit (dof_handler.n_dofs());
	system_rhs.reinit (dof_handler.n_dofs());
	std::ofstream out ("sparsity_pattern.1");
	sparsity_pattern.print_gnuplot(out);
}

template <int dim>
void transThermal<dim>::assemble_system(){

}
