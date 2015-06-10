/*
 * transientCdns.cpp
 *
 *  Created on: Apr 10, 2015
 *      Author: karthik
 */
#include "thermalSimulation.hpp"

template <int dim>
void transThermal<dim>::heatingBndCdn(){
/*Dirichlet boundary cdn T = 80 at bottom surface
	and T = 1255 at top layer*/
	// Interpolating the boundary values
	std::map<types::global_dof_index,double> boundary_values;
	VectorTools::interpolate_boundary_values (dof_handler,
			bottom_Bnd,
			ConstantFunction<dim>(80),
			boundary_values);
	VectorTools::interpolate_boundary_values (dof_handler,
			currentLayer+layer_Bnd_offset,
			ConstantFunction<dim>(1255),
			boundary_values);

	// Modifying the system of equations
	MatrixTools::apply_boundary_values (boundary_values,
			system_matrix,
			solution,
			system_rhs);
}

template <int dim>
void transThermal<dim>::coolingBndCdn(){
	// Interpolating the boundary values
	std::map<types::global_dof_index,double> boundary_values;
	VectorTools::interpolate_boundary_values (dof_handler,
			bottom_Bnd,
			ConstantFunction<dim>(80),
			boundary_values);

	Vector<double> forcing_terms;
    forcing_terms.reinit(dof_handler.n_dofs());
	QGauss<dim-1> feQuadrature(activeFE.degree+1);
	QGauss<dim-1> zeroQuadrature(1);
	hp::QCollection<dim-1>  q_collection;
	q_collection.push_back(zeroQuadrature);
	q_collection.push_back(feQuadrature);
	int layers[] = {bottom_Bnd,currentLayer+layer_Bnd_offset};
	// Neumann boundary condition on the rhs for both the top and bottom layers need to implement a function for now a zeroFunction has been considered
	VectorTools::create_boundary_right_hand_side (dof_handler,
	                                              q_collection,
	                                              ZeroFunction<dim>(),
	                                              forcing_terms, std::set<types::boundary_id>(layers, layers+2));
	system_rhs += forcing_terms;

	// Modifying the system of equations
	MatrixTools::apply_boundary_values (boundary_values,
			system_matrix,
			solution,
			system_rhs);
}

template <int dim>
void transThermal<dim>::pretimeStepping(){
	assemble_system(heatingTime);
	heatingBndCdn();
	if (solveTag)
		solve_system();
	std::cout << "Solution vector size" << solution.size()<<std::endl;
	time += heatingTime;
	old_solution = solution;
	store_PrevDOFs();
}

template <int dim>
void transThermal<dim>::postCoolingStep(){
	// TODO: Need to create adaptive time-stepping scheme later
	setup_matrix();
	for(double localTime = 0.0; localTime < 14.0; localTime+=coolingTime){
		assemble_system(coolingTime);
		coolingBndCdn();
		solve_system();
		time += coolingTime;
		old_solution = solution;
	}
	store_PrevDOFs();
}
