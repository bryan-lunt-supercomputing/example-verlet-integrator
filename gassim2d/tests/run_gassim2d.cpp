#include <random>

#include "gassim2d.h"
#include "gassim2d_util.h"

#include <string>
#include <sstream>


#define SIGMA_AR_NM 0.34
#define EPSILON_AR_AGNM2PERNS2 1.65
#define KB_AGPERKELVIN_NM2PERS2 1.38064852E-2
#define MASS_AR_AG 6.6335209E-5

int main(int argc, char **argv){

	lattice_t fill_this = create_lattice(0.36,
											{.x=-1.0,.y=-1.0},
											{.x=2.0,.y=2.0}
										);

	int N_particles = fill_this.N_total;;
	int N_ghosts = 0;


	int output_at_iteration = 1000;
	double timestep = 0.0000001;//timestep is given in nanoseconds, this is 1/10th of an ns
	double N_steps = 100000;
	double T_total = timestep*N_steps;

	gas_simulation my_sim;

	phys_particle_t all_particles[N_particles];
	phys_particle_t all_ghosts[N_ghosts];
	phys_vector_t forces[N_particles];

	my_sim = (gas_simulation) { .mass = MASS_AR_AG, .epsilon = EPSILON_AR_AGNM2PERNS2, .sigma = SIGMA_AR_NM, .cutoff = 3.0, .v_max = 1000.0 };

	int j;
	for(j = 0;j<N_particles;j++){
		phys_particle_t p;
		lattice_position_t foobar = get_lattice_position(&fill_this,j);
		p.p = (phys_vector_t){ .x = foobar.x, .y = foobar.y };
		p.v = (phys_vector_t){ .x = 0.0, .y = 0.0 };
		all_particles[j] = p;
	}

	maxwell_boltzmann(&my_sim, 1.0*KB_AGPERKELVIN_NM2PERS2, all_particles, N_particles );

	/* STARTING VELOCITIES */
	if( true ){
	fprintf(stderr, "%d\n",N_particles);//Number of particles
	fprintf(stderr, "\n");//Comment line
	for(int j = 0;j<N_particles;j++)
		fprintf(stderr, "Ar %f %f\n",j, all_particles[j].v.x, all_particles[j].v.y);
	}


	printf("%d\n",N_particles);//Number of particles
	printf("\n");//Comment line
	for(int j = 0;j<N_particles;j++)
		printf("Ar %f %f\n",j, all_particles[j].p.x, all_particles[j].p.y);




	int i;
	double T_elapsed = 0.0;
	for(i=0;i<N_steps;i++){

		find_forces_environment(&my_sim, all_particles, N_particles, forces);
		find_forces_self(&my_sim, all_particles, N_particles, forces);
		find_forces_ghost(&my_sim, all_particles, N_particles, forces, all_ghosts, N_ghosts);
		apply_relativity(&my_sim, all_particles, N_particles);
		update_states(&my_sim, timestep, all_particles, N_particles, forces);

		T_elapsed += timestep;
		if(0 == i % output_at_iteration){
			fprintf(stderr, "f= %d \n", i);
			printf("%d\n",N_particles);//Number of particles
			printf("\n");//Comment line
			for(int j = 0;j<N_particles;j++)
				printf("Ar %f %f\n",j, all_particles[j].p.x, all_particles[j].p.y);
		}
	}

	return 0;
}
