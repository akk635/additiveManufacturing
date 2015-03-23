/*
 * thermalSimulation.hpp
 *
 *  Created on: Mar 1, 2015
 *      Author: karthik
 */

#ifndef THERMALSIMULATION_HPP_
#define THERMALSIMULATION_HPP_

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/numerics/data_out.h>

#include <iostream>
#include <fstream>

using namespace dealii;

template <int dim>
class transThermal{
private:
	enum{
		zeroFE_id,
		activeFE_id

	};
	Triangulation<dim> &triangulation;
	hp::DoFHandler<dim> dof_handler; // to enumerate the DOF's in a mesh
	hp::FECollection<dim> fe_collection; // to ensure different finite elements for different cells

	// Finite element basis function to represent the active and zero finite elements
	FE_Q<dim> activeFE;
	FE_Nothing<dim> zeroFE;
	SparsityPattern      sparsity_pattern;
	SparseMatrix<double> system_matrix;
	Vector<double>       solution;
	Vector<double>       system_rhs;
	int nZ = 0; // no. of layers in Z-direction

	double time = 0.0;
	// Helper functions for layer wise increment
	typename hp::DoFHandler<dim>::active_cell_iterator topLayerCell;
	typename hp::DoFHandler<dim>::active_cell_iterator currentLayerCell;
	typedef typename hp::DoFHandler<dim>::active_cell_iterator cellsIterator;

	std::set<cellsIterator>* layerIterator;

public:
	transThermal(Triangulation<dim> & tria);
	~transThermal(){};
	void readFullMesh();
	void defaultSetup(); // Setting the cell active fe index and also initialize the top and current layers
	void update_active_fe_indices();
	void incrementLayer();
	void setup_system();
	void assemble_system();
	void solve_system(){};
};

#include "thermalSimulation_Opt.cpp"

#endif /* THERMALSIMULATION_HPP_ */
