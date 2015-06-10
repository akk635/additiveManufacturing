/*
 * thermalSimulation.hpp
 *
 *  Created on: Mar 1, 2015
 *      Author: karthik
 */

#ifndef THERMALSIMULATION_HPP_
#define THERMALSIMULATION_HPP_

// classes handling the meshes and enumeration of DOF's
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
// accessing the cells information
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>
// classes for the finite element type creation
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_nothing.h>
// creates the sparsitypattern for the sparse matrices
#include <deal.II/dofs/dof_tools.h>
// for assembling the linear system
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
// for the treatment of the boundary values
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
// for the linear algebra operations
#include <deal.II/lac/vector.h>
#include <deal.II/base/aligned_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/numerics/data_out.h>
#include <iostream>
#include <fstream>
#include "helper.hpp"

using namespace dealii;

template <int dim>
class transThermal{
private:
	enum{
		zeroFE_id,
		activeFE_id
	};
	// enums for boundary index
	enum{
		default_Bnd,
		bottom_Bnd,
		layer_Bnd_offset
	};
	Triangulation<dim> &triangulation;
	hp::DoFHandler<dim> dof_handler; // to enumerate the DOF's in a mesh
	hp::FECollection<dim> fe_collection; // to ensure different finite elements for different cells

	// Finite element basis function to represent the active and zero finite elements
	FE_Q<dim> activeFE;
	FE_Nothing<dim> zeroFE;
	SparsityPattern      sparsity_pattern;
	SparseMatrix<double> mass_matrix;
	SparseMatrix<double> laplace_matrix;
	SparseMatrix<double> system_matrix;
	Vector<double>       solution;
	Vector<double>       old_solution;
	Vector<double>       system_rhs;
	int currentLayer; // Numbering started from zero

    double               time = 0.0;
    unsigned int         timestep_number;
    // For the time stepping scheme
    const double         theta;
    const double heatingTime = 0.6;
    double coolingTime = 0.6;

    // Material parameters in SI units
    // TODO: shd consider the non-linear parameters
    int rho = 8390;
    int cP = 370;
    int kT = 6;
	// Helper functions for layer wise increment
	typename hp::DoFHandler<dim>::active_cell_iterator currentLayerCell;
	typedef typename hp::DoFHandler<dim>::active_cell_iterator cellsIterator;
	std::vector<cellsIterator>* layerIterator;
	std::vector<std::vector<types::global_dof_index> > prev_DOFs; // To store the DOF indices of the previous solution
	std::vector<double> checkItrOrder;

public:
	transThermal(Triangulation<dim> & tria);
	~transThermal(){};
	int nZ = 0; // total no. of layers in Z-direction
    bool solveTag = true;
	void readFullMesh();
	void assignLayerIterators();
	void defaultSetup(); // Setting the cell active fe index and also initialize the top and current layers
	void update_active_fe_indices();
	void setup_system();
	void setup_DOFs();
	void setup_matrix(); // to accommodate the change in the boundary conditions during the heating and cooling conditions
	void setup_solution(); // If it is not a symmetric sparsity pattern then it needs to be implemented
	void store_PrevDOFs();
	void assemble_system(double time_step);
	void heatingBndCdn();
	void coolingBndCdn();
	void solve_system();
	void outputResults(const std::string filename);
	// Running steps
	void pretimeStepping(); // during the actual laser heating for 600 ms
	void postCoolingStep(); // Cooling for 14 secs after the laser heating
};



template<unsigned int dim, typename T=double>
class matrixCoefficient : public Function<dim,T>
{
private:
	 T coefficient;
public:
	 matrixCoefficient(T argValue):Function<dim, T>(){coefficient = argValue;}
	 ~matrixCoefficient(){}
  virtual T value (const Point<dim>  &p,
                        const unsigned int component = 0) const{
	  Assert(component == 0, ExcInternalError());
	  return coefficient;
  }
};

#include "thermalSimulation_Opt.cpp"
#include "transientCdns.cpp"

#endif /* THERMALSIMULATION_HPP_ */
