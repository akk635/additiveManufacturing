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
	cellsIterator firstLayerCell = currentLayerCell;

	std::vector<double> height;
	height.push_back((double)currentLayerCell->center()[dim-1] * 10e4);



	while(1){
		nZ++;
		if (currentLayerCell->face(dim-1)->at_boundary())
			break;
		currentLayerCell = currentLayerCell->neighbor(dim-1);
		height.push_back((double)currentLayerCell->center()[dim-1]* 10e4);
	}

	int layerIndex = 0;
	layerIterator =  new std::vector<cellsIterator> [nZ];
	std::vector<double>::iterator layer;
	int i = 0;
	// Setting the initial fe index of cells
	for (typename hp::DoFHandler<dim>::active_cell_iterator
			cell = dof_handler.begin_active();
			cell != dof_handler.end(); cell++){
		cell->set_active_fe_index(0); // assigning zero FE for all the cells
		layer = std::lower_bound(height.begin(),height.end(),(double) cell->center()[dim-1]* 10e4);
		layerIndex = (int)(layer - height.begin());
		std::cout<< typeid(layer-height.begin()).name()<< std::endl;
		std::cout<< "layer" << layerIndex;
		std::cout<< " height"<<cell->center()[dim-1]* 10e4 << std::endl;
		layerIterator[layerIndex].push_back(cell);
	}
	std::cout<<'0'<<std::endl;
	// TODO: boundary indices
	for (auto cell = layerIterator[0].begin();
			cell != layerIterator[0].end(); ++cell){
		(*cell)->face(dim-1)->set_boundary_indicator(1);
		std::cout<<'1'<<std::endl;
	}
}

template <int dim>
void transThermal<dim>::update_active_fe_indices(){
	int i = 0;
	// Updating the active cell fe indexes
	for (auto cell = layerIterator[0].begin();
			cell != layerIterator[0].end(); ++cell){
		i++;
		(*cell)->set_active_fe_index(0);
		std::cout<<i<<std::endl;
		// assigning the bnd flags
		(*cell)->face(dim-1)->set_boundary_indicator(1);
	}
	std::cout<<"Initial layer cells" <<"\t"<< i << std::endl;
}

template <int dim>
void transThermal<dim>::setup_system(){
	dof_handler.distribute_dofs(fe_collection);
	CompressedSparsityPattern compressed_sparsity_pattern(dof_handler.n_dofs());
	DoFTools::make_sparsity_pattern (dof_handler, compressed_sparsity_pattern);
	sparsity_pattern.copy_from (compressed_sparsity_pattern);
	// matrix definitions valid for serial computations
	mass_matrix.reinit(sparsity_pattern);
	laplace_matrix.reinit(sparsity_pattern);
	system_matrix.reinit (sparsity_pattern);
/*	MatrixCreator::create_mass_matrix(dof_handler,
			QGauss<dim>(activeFE.degree+1),
			mass_matrix); // function pointer allows to specify the coefficient of the matrix entry
	MatrixCreator::create_laplace_matrix(dof_handler,
			QGauss<dim>(activeFE.degree+1),
			laplace_matrix);*/
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
	FEValues<dim> fe_values (activeFE, quadrature_formula,
			update_values | update_gradients | update_JxW_values); // encompassing objects of FE, Quadrature and mappings
	const unsigned int   dofs_per_cell = activeFE.dofs_per_cell;
	const unsigned int   n_q_points    = quadrature_formula.size();

	FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
	Vector<double>       cell_rhs (dofs_per_cell);

	std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

	for(auto cell : triangulation.active_cell_iterators()){
		fe_values.reinit(); // recomputes FE_Values on each cell
		// resetting the local contributions
		//TODO: to change when it can be computed parallely
		cell_matrix = 0;
		cell_rhs = 0;
		// integration over the quadrature points
		for (unsigned int q_index=0; q_index<n_q_points; ++q_index){
			// local laplace matrix computation
			for (unsigned int i=0; i<dofs_per_cell; ++i)
				for (unsigned int j=0; j<dofs_per_cell; ++j)
					cell_matrix(i,j) += (fe_values.shape_grad (i, q_index) *
							fe_values.shape_grad (j, q_index) *
							fe_values.JxW (q_index));
			// local rhs computation
			for (unsigned int i=0; i<dofs_per_cell; ++i)
				// TODO: to include the effects of the dirichlet bnd cdns
				cell_rhs(i) += 0;
		}
		// to get the corresponding global matrix indices
		cell->get_dof_indices (local_dof_indices);
		for (unsigned int i=0; i<dofs_per_cell; ++i)
			for (unsigned int j=0; j<dofs_per_cell; ++j)
				system_matrix.add (local_dof_indices[i],
						local_dof_indices[j],
						cell_matrix(i,j));
		for (unsigned int i=0; i<dofs_per_cell; ++i)
			system_rhs(local_dof_indices[i]) += cell_rhs(i);
	}
	// Interpolating the boundary values and proper boundary variables
	std::map<types::global_dof_index,double> boundary_values;
	VectorTools::interpolate_boundary_values (dof_handler,
			0,
			ConstantFunction<dim>(80),
			boundary_values);

	// Modifying the system of equations
	MatrixTools::apply_boundary_values (boundary_values,
			system_matrix,
			solution,
			system_rhs);*/
}

template <int dim>
void transThermal<dim>::outputResults(){
	DataOut<dim> data_out;
}
