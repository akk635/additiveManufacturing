/*
 * thermalSimulation.cpp
 *
 *  Created on: Mar 1, 2015
 *      Author: karthik
 */
#include "thermalSimulation.hpp"

template <int dim>
transThermal<dim>::transThermal(Triangulation<dim> & tria):triangulation(tria), dof_handler(tria),activeFE(1), theta(0.5){
	fe_collection.push_back(zeroFE);
	fe_collection.push_back(activeFE);
}

// To read the external mesh from salome platform
// read the deal.II documentation to know more about the GridIn class
template <int dim>
void transThermal<dim>::readFullMesh(){
	GridIn<dim> grid_in;
	grid_in.attach_triangulation (triangulation);
	std::ifstream input_file("./MeshFiles/60Structure.unv");
	grid_in.read_unv (input_file);

	std::cout << "   Number of active cells: "
			<< triangulation.n_active_cells()
			<< std::endl
			<< "   Total number of cells: "
			<< triangulation.n_cells()
			<< std::endl;
}

template <int dim>
void transThermal<dim>::assignLayerIterators(){
	// Initial additive manufacturing layers
	currentLayerCell = dof_handler.begin_active();
	cellsIterator firstLayerCell = currentLayerCell;

	std::vector<double> height;
	height.push_back((double)currentLayerCell->center()[dim-1]); // In this case (dim-1) specifies Z-direction
	while(1){
		nZ++;
		if (currentLayerCell->face(dim-1)->at_boundary())
			break;
		currentLayerCell = currentLayerCell->neighbor(dim-1);
		height.push_back((double)currentLayerCell->center()[dim-1]);
	}
	layerIterator =  new std::vector<cellsIterator> [nZ];

	int layerIndex = 0;
	std::vector<double>::iterator layer;
	// Comparison function to compare the double values
	compareFunction<double> compDouble((int)1e4);
	// Setting the initial fe index of cells
	for (typename hp::DoFHandler<dim>::active_cell_iterator
			cell = dof_handler.begin_active();
			cell != dof_handler.end(); cell++){
		cell->set_active_fe_index(0); // assigning zero FE for all the cells
		layer = std::lower_bound(height.begin(),height.end(),(double) cell->center()[dim-1],compDouble);
		layerIndex = (int)(layer-height.begin());
		layerIterator[layerIndex].push_back(cell);
	}
}

template <int dim>
void transThermal<dim>::defaultSetup(){
	// Updating the active cell fe indexes of first layer
	currentLayer = 0;
	for (auto cell = layerIterator[currentLayer].begin();
			cell != layerIterator[currentLayer].end(); ++cell){
		(*cell)->set_active_fe_index(1); // 2 times de-referencing of the cell because the iterator is one hop away from the data
		(*cell)->face(dim)->set_boundary_indicator(bottom_Bnd); // face(int face_Index) follow a specific numbering pattern
	}
/*	currentLayer++;
	for (auto cell = layerIterator[currentLayer].begin();
			cell != layerIterator[currentLayer].end(); ++cell){
		(*cell)->set_active_fe_index(1); // 2 times de-referencing of the cell because the iterator is one hop away from the data
		(*cell)->face(dim-1)->set_boundary_indicator(currentLayer + layer_Bnd_offset);
	}*/
	std::cout << "defaultSetup done" << std::endl;
}

template <int dim>
void transThermal<dim>::update_active_fe_indices(){
	currentLayer++;
	if (currentLayer == nZ) return;
	for (auto cell = layerIterator[currentLayer].begin();
			cell != layerIterator[currentLayer].end(); ++cell){
		(*cell)->set_active_fe_index(1);
		(*cell)->face(dim-1)->set_boundary_indicator(currentLayer + layer_Bnd_offset);
	}
}

template <int dim>
void transThermal<dim>::setup_system(){
	setup_DOFs();
	setup_matrix();
	setup_solution();
	system_rhs.reinit (dof_handler.n_dofs());
	std::cout << "Setup done" << std::endl;
}

template <int dim>
void transThermal<dim>::setup_DOFs(){
	dof_handler.distribute_dofs(fe_collection);
	CompressedSparsityPattern compressed_sparsity_pattern(dof_handler.n_dofs());
	DoFTools::make_sparsity_pattern (dof_handler, compressed_sparsity_pattern);
	sparsity_pattern.copy_from (compressed_sparsity_pattern);
}

template <int dim>
void transThermal<dim>::setup_matrix(){
	mass_matrix.reinit(sparsity_pattern);
	laplace_matrix.reinit(sparsity_pattern);
	system_matrix.reinit (sparsity_pattern);
	// create the mass and laplace matrices
	QGauss<dim> feQuadrature(activeFE.degree+1);
	QGauss<dim> zeroQuadrature(1);
	hp::QCollection<dim>  q_collection;
	q_collection.push_back(zeroQuadrature);
	q_collection.push_back(feQuadrature);

	matrixCoefficient<dim> temp = matrixCoefficient<dim>(rho*cP);
	MatrixCreator::create_mass_matrix(dof_handler,
			q_collection,
			mass_matrix,
			(const Function<dim> *)& temp); // default coeff. of matrix entries
	matrixCoefficient<dim> temp1 = matrixCoefficient<dim>(kT);
	MatrixCreator::create_laplace_matrix(dof_handler,
			q_collection,
			laplace_matrix,
			(const Function<dim> *)& temp1);
}

template <int dim>
void transThermal<dim>::setup_solution(){
	// Based upon the sparsity pattern may be needed to be implemented
	solution.reinit(dof_handler.n_dofs());
	std::cout << "# DOF's" << dof_handler.n_dofs() << std::endl;
	// Reinitializing the previous solutions
	int cellIndex;
	if (old_solution.size() > 0){
		std::cout << "CurrentLayer" << currentLayer << std::endl;
		for (int i = 0; i < currentLayer; i++){
			int j = 0;
			for (auto cell = layerIterator[i].begin();
					cell != layerIterator[i].end(); ++cell){
				cellIndex = i * layerIterator[i].size() + j++;
				std::vector<types::global_dof_index> tempDOF((*cell)->get_fe().n_dofs_per_cell());
				(*cell)->get_dof_indices(tempDOF);
				assert(checkItrOrder[cellIndex] == (*cell)->center()[0]);
				for (int dof = 0; dof < tempDOF.size(); dof++){
					assert(tempDOF[dof] < dof_handler.n_dofs());
					solution[tempDOF[dof]] = old_solution[prev_DOFs[cellIndex][dof]]; // Assuming same cell traversal for iterators
				}
			}
		}
	}
	old_solution.reinit(dof_handler.n_dofs());
	old_solution = solution;
}

template <int dim>
void transThermal<dim>::store_PrevDOFs(){
	// old_solution has been properly set
	for (int i = 0; i <= currentLayer; i++){
		for (auto cell = layerIterator[i].begin();
				cell != layerIterator[i].end(); ++cell){
			std::vector<types::global_dof_index> tempDOF((*cell)->get_fe().n_dofs_per_cell());
			(*cell)->get_dof_indices(tempDOF);
			prev_DOFs.push_back(tempDOF);
			checkItrOrder.push_back((*cell)->center()[0]);
		}
	}
}

template <int dim>
void transThermal<dim>::assemble_system(double time_step){
	Vector<double> tmp;
	tmp.reinit(dof_handler.n_dofs());
	// \left( M + k_n \Theta A\right) U^n = M U^(n-1) - k_n (1-\Theta) A U^(n-1) + k_n \left[ (1- \Theta)F^(n-1)+ \Theta F^(n)\right]
	// The source terms F in for this particular problem is zero
	mass_matrix.vmult(system_rhs, old_solution);
	laplace_matrix.vmult(tmp, old_solution);
	system_rhs.add(-(1 - theta) * time_step, tmp);

	system_matrix.copy_from(mass_matrix);
	system_matrix.add(theta * time_step, laplace_matrix);
}

template <int dim>
void transThermal<dim>::outputResults(const std::string filename){
	DataOut<dim,hp::DoFHandler<dim> > data_out;
	data_out.attach_dof_handler(dof_handler);
	data_out.add_data_vector(solution, "T");
	data_out.build_patches();
	std::ofstream output(filename.c_str());
	data_out.write_vtk(output);
}

template <int dim>
void transThermal<dim>::solve_system(){
	SolverControl solver_control(1000, 1e-8 * system_rhs.l2_norm());
	SolverCG<> cg(solver_control);
	PreconditionSSOR<> preconditioner;
	preconditioner.initialize(system_matrix, 1.0);
	cg.solve(system_matrix, solution, system_rhs,
			preconditioner);

	/*std::cout << "     " << solver_control.last_step()
            								  << " CG iterations." << std::endl;*/
}
