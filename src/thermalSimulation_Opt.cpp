/*
 * thermalSimulation.cpp
 *
 *  Created on: Mar 1, 2015
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

	std::vector<double> height;
	height.push_back(currentLayerCell->center()[dim-1]);

	while(1){
		nZ++;
		if (currentLayerCell->face(dim-1)->at_boundary())
			break;
		currentLayerCell = currentLayerCell->neighbor(dim-1);
		height.push_back(currentLayerCell->center()[dim-1]);
	}

	layerIterator =  new typename hp::DoFHandler<dim>::active_cell_iterator[nZ];
	// Setting the initial fe index of cells
	for (typename hp::DoFHandler<dim>::active_cell_iterator
	       cell = dof_handler.begin_active();
	       cell != dof_handler.end(); ++cell){
		cell->set_active_fe_index (0); // assigning zero FE for all the cells
		*(layerIterator + std::lower_bound(height.begin(),height.end(),cell->center()[dim-1])).insert(cell);
	}
}

template <int dim>
void transThermal<dim>::update_active_fe_indices(){
	int i = 0;
	// Updating the active cell fe indexes
	for (typename hp::DoFHandler<dim>::active_cell_iterator
		       cell = layerIterator[0].begin_active();
		       cell != layerIterator[0].end(); ++cell){
		i++;
	}
	std::cout<<"Initial layer cells" <<"\t"<< i << std::endl;
}

template <int dim>
void transThermal<dim>::setup_system(){
	dof_handler.distribute_dofs(fe_collection);
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
void transThermal<dim>::incrementLayer(){
	if(time == 0.0)
		topLayerCell = dof_handler.begin_active();

	if (topLayerCell->face(dim-1)->at_boundary())
		return;

	topLayerCell = topLayerCell->neighbor(dim-1);
}

template <int dim>
void transThermal<dim>::assemble_system(){
	/*	QGauss<dim>  quadrature_formula(2); // 2 qpoints in each direction
	FEValues<dim> fe_values (finite_element, quadrature_formula,
			update_values | update_gradients | update_JxW_values);
	const unsigned int   dofs_per_cell = finite_element.dofs_per_cell;
	const unsigned int   n_q_points    = quadrature_formula.size();

	FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
	Vector<double>       cell_rhs (dofs_per_cell);

	std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
	DoFHandler<dim>::active_cell_iterator
	cell = dof_handler.begin_active(),
	endc = dof_handler.end();
	for (; cell!=endc; ++cell){
		fe_values.reinit(cell);
		cell_matrix = 0;
		cell_rhs = 0;

	}*/
}
