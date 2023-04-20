
#include "gassim.h"

#include <stdlib.h>
//#include <stdio.h>
#include <math.h>

void clear_forces(const int num_particles, phys_vector_t *forces){
	int i;
	for(i=0;i<num_particles;i++){
		forces[i] = ZERO_VECTOR;
	}
}


void find_forces_environment(const gas_simulation *sim, phys_particle_t *gas, const int num_particles, phys_vector_t *forces){
	//Currently not implemented.
	return clear_forces(num_particles, forces);
}

void find_forces_self(const gas_simulation *sim, phys_particle_t *gas, const int num_particles, phys_vector_t *forces){
	int i, j;
	for(i = 0;i<num_particles;i++){
		for(j = i+1;j<num_particles;j++){
			phys_vector_t r = VEC_DIFF(gas[i].p,gas[j].p);
			sim_working_t r_mag = VEC_NORM(r);
			if( r_mag > sim->cutoff ) { continue; }
			sim_working_t force_ij = F_LJ(r_mag, *sim); //felt by j, caused by i
			//printf("FOUND A FORCE %e \n",force_ij);
			force_ij = pow(r_mag,-1.0)*force_ij;//So that we don't need to normalize r.
						/*For efficiency we would put that into F_LJ (because it already takes a power of r_mag),
						but right now we value correctness
						above all*/
			forces[i] = VEC_P_SCALtVEC( forces[i], -force_ij, r);
			forces[j] = VEC_P_SCALtVEC( forces[j], +force_ij, r);
			/*
			forces[i] = (phys_vector_t){ .x = forces[i].x - force_ij*r.x,
										.y = forces[i].y - force_ij*r.y,
										.z = forces[i].z - force_ij*r.z };

			forces[j] = (phys_vector_t){ .x = forces[j].x + force_ij*r.x,
										.y = forces[j].y + force_ij*r.y,
										.z = forces[j].z + force_ij*r.z };
			*/
		}
	}
}

void find_forces_ghost(const gas_simulation *sim, phys_particle_t *gas, const int num_particles, phys_vector_t *forces, phys_particle_t *ghost, const int num_ghost ){
	int i, j;
	for(i = 0;i<num_particles;i++){
		for(j = 0;j<num_ghost;j++){
			phys_vector_t r = VEC_DIFF(ghost[j].p,gas[i].p);
			sim_working_t r_mag = VEC_NORM(r);
			if( r_mag > sim->cutoff ) { continue; }
			sim_working_t force_ji = F_LJ(r_mag, *sim); //felt by i, caused by j
			force_ji = pow(r_mag,-1.0)*force_ji;//So that we don't need to normalize r.
						/*For efficiency we would put that into F_LJ (because it already takes a power of r_mag),
						but right now we value correctness
						above all*/

			forces[i] = VEC_P_SCALtVEC(forces[i],force_ji,r);
			/* forces[i] = (phys_vector_t){ .x = forces[i].x + force_ji*r.x,
										.y = forces[i].y + force_ji*r.y,
										.z = forces[i].z + force_ji*r.z };
			*/
		}
	}
}


void apply_relativity(const gas_simulation *sim, phys_particle_t *gas, const int num_particles){
	//Does not actually apply relativity.
	int i;
	for(i = 0;i< num_particles;i++){
		sim_working_t v_mag = VEC_NORM(gas[i].v);
		if( v_mag <= sim->v_max){ continue; }
		sim_working_t norm_factor = sim->v_max * pow(v_mag, -1.0);

		gas[i].v = SCALtVEC(norm_factor, gas[i].v);
		/*
		gas[i].v = (phys_vector_t){ .x = norm_factor*gas[i].v.x,
									.y = norm_factor*gas[i].v.y,
									.z = norm_factor*gas[i].v.z };
		*/
	}
}

void update_states(const gas_simulation *sim, const double timestep, phys_particle_t *gas, const int num_particles, phys_vector_t *forces) {

	/* When there are many particles, everything up here can be considered as
	"free". Do whatever calculation is possible here. Value clarity.
	*/
	const sim_working_t half_dt_squared = 0.5*pow(timestep,2.0);
	const sim_working_t time_over_mass = timestep*pow(sim->mass,-1.0); //Everything would be faster if we take the masses and radii to be 1.0 each...
	const sim_working_t local_half_timesqr_over_mass = 0.5*pow(timestep,2.0)*pow(sim->mass,-1.0);
	//const double half_t = 0.5*timestep;
	int i;
	for(i = 0;i< num_particles;i++){
		//X(t) = (1/2)a*dt^2 + V_0*dt + X_0
		//X(t) = (1/2)(F/m)*dt^2 + V_0*dt + X_0
		//local_half_timesqr_over_mass = (1/2)(dt^2)/(m)
		//So,
		//X(t) = lhtom*F + V_o*dt + X_0

		phys_vector_t newpos_partial = VEC_P_SCALtVEC(gas[i].p, timestep, gas[i].v); //displacement just because of velocity
		gas[i].p = VEC_P_SCALtVEC(newpos_partial, local_half_timesqr_over_mass, forces[i]); //aditional displacement from velocity gained through acceleration

		/*
		const double a_dx = forces[i].x * local_half_timesqr_over_mass;
		const double a_dy = forces[i].y * local_half_timesqr_over_mass;
		const double a_dz = forces[i].z * local_half_timesqr_over_mass;
		gas[i].p = (phys_vector_t){ .x = gas[i].p.x + timestep*gas[i].v.x + a_dx,
									.y = gas[i].p.y + timestep*gas[i].v.y + a_dy,
									.z = gas[i].p.z + timestep*gas[i].v.z + a_dz };
		*/

		gas[i].v = VEC_P_SCALtVEC(gas[i].v, time_over_mass, forces[i]);//velocity changes according to acceleration F=mA -> DF = mA*time
		/*
		const double dv_x = forces[i].x*time_over_mass;
		const double dv_y = forces[i].y*time_over_mass;
		const double dv_z = forces[i].z*time_over_mass;
		gas[i].v = (phys_vector_t){ .x = gas[i].v.x + dv_x,
									.y = gas[i].v.y + dv_y,
									.z = gas[i].v.z + dv_z };
		*/

	}
}

/**
*This uses the Marsaglia algorithm, which generates a pair of independent gaussian distributed doubles at a time.
*
*Code adapted (taken directly) from: https://rosettacode.org/wiki/Statistics/Normal_distribution#C
*/
void marsaglia_alg(double *destination, int even_numtogen){
	int i;
	for (i = 0; i < even_numtogen; i += 2 )
        {
            double x,y,rsq,f;
            do {
                x = 2.0 * rand() / (double)RAND_MAX - 1.0;
                y = 2.0 * rand() / (double)RAND_MAX - 1.0;
                rsq = x * x + y * y;
            }while( rsq >= 1. || rsq == 0. );
            f = sqrt( -2.0 * log(rsq) / rsq );
            destination[i]   = x * f;
            destination[i+1] = y * f;
    	}
}


/**
* Set the velocities of the particles to random values according to a Maxwell-Boltzmann distribution.
*
* You need to provide your own value for Boltzmann's constant. This is so that you can set the mass, temp, etc. in arbitrary units.
* Our program will use nanometers / nanosecond = m/s and Kelvins, but a user might want something else.
*/
void maxwell_boltzmann(const gas_simulation *sim, const double kT, phys_particle_t *gas, const int num_particles ){
	if( num_particles <= 0){ return; }
	//there is at least one particle.

	/* The marsaglia algorithm generates multiple doubles at a time.
		This loop prevents us from wasting them by keeping a buffer/queue of
		random doubles. We get doubles out of that and into the velocities.
	*/

	const double maxwell_boltzmann_stddev = pow(kT*pow(sim->mass,-1.0),0.5);//The standard deviation of velocities under the Maxwell-Boltzmann distribution.
	
	//Variables related to batching the output of the Marsaglia algorithm. It
	const int batch_size = 2;
	double tmp_dataspace[VECTOR_DIMENSIONALITY*batch_size];//we generate enough for two vectors at once. (Guaranteeing we generate an even number and don't waste.)
	int batch_offset = 0;
	
	int i;//loop iterator
	for(i=0;i<num_particles;i++){
		
		if( i%batch_size == 0){
			//generate a new batch of normal distribution doubles.
			//We will scale for the stddev later.
			batch_offset = 0;
			marsaglia_alg(tmp_dataspace, VECTOR_DIMENSIONALITY*batch_size );
		}

		phys_vector_t tmp;
		tmp.x = tmp_dataspace[batch_offset++];
		tmp.y = tmp_dataspace[batch_offset++];
		tmp.z = tmp_dataspace[batch_offset++];
		
		gas[i].v = SCALtVEC(maxwell_boltzmann_stddev, tmp); //scale for stddev.
	}
}
