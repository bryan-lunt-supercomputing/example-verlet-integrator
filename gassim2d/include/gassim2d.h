#ifndef GASSIM_H
#define GASSIM_H


#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef double sim_working_t;//A lot of float hardware is internally 80bit, use that when available.
typedef double sim_scalar_t;

typedef struct { sim_scalar_t x; sim_scalar_t y;} phys_vector_t;
typedef struct { phys_vector_t p; phys_vector_t v; } phys_particle_t;

typedef struct { phys_vector_t min; phys_vector_t max; } bounding_box_t;

typedef struct { int x; int y; } int_vector_t; //C bools and C++ bools incompatible.
typedef struct { sim_scalar_t mass; sim_scalar_t epsilon; sim_scalar_t sigma; sim_scalar_t cutoff; sim_scalar_t v_max; bounding_box_t bounds; int_vector_t enforce_bounds; } gas_simulation;


#define VECTOR_DIMENSIONALITY 2
#define ZERO_VECTOR (phys_vector_t){ .x = 0.0, .y = 0.0};

inline phys_vector_t vec_difference(const phys_vector_t a, const phys_vector_t b){
	return (phys_vector_t){.x = a.x - b.x , .y = a.y - b.y};
}
#define VEC_DIFF(A,B) (phys_vector_t){.x = A.x - B.x , .y = A.y - B.y}

inline sim_scalar_t vec_norm(const phys_vector_t a){
	return sqrt(pow(a.x,2.0) + pow(a.y,2.0));
}
#define VEC_NORM(A) sqrt(pow(A.x,2.0) + pow(A.y,2.0))

#define VEC_P_VEC(A,B) (phys_vector_t){ .x = A.x + B.x, .y = A.y + B.y}
#define VEC_P_SCALtVEC(A,S,B) (phys_vector_t){ .x = A.x + S*B.x, .y = A.y + S*B.y}
#define SCALtVEC(S,A) (phys_vector_t){ .x = S*A.x, .y = S*A.y}

inline sim_working_t potential_LJ(const sim_working_t r, const gas_simulation sim){
	sim_working_t sixth_power = pow(sim.sigma*pow(r,-1.0),6.0);
	sim_working_t twelth_power = pow(sixth_power,2.0);
	return 4*sim.epsilon*(twelth_power - sixth_power);
}

inline sim_working_t force_LJ(const sim_working_t r, const gas_simulation sim){
	sim_working_t sixth_power = pow(sim.sigma*pow(r,-1.0),6.0);
	sim_working_t twelth_power = pow(sixth_power,2.0);
		/*
		It might be better to just do the full pow, which could exploit hardware parallelism...
		*/
	return -48.0*sim.epsilon*pow(r,-1.0)*(twelth_power - 0.5*sixth_power);
}

#define U_LJ(r_mag, s) potential_LJ(r_mag, s)
#define F_LJ(r_mag, s) force_LJ(r_mag, s)


void clear_forces(const int num_particles, phys_vector_t *forces);
void find_forces_environment(const gas_simulation *sim, phys_particle_t *gas, const int num_particles, phys_vector_t *forces);
void find_forces_self(const gas_simulation *sim, phys_particle_t *gas, const int num_particles, phys_vector_t *forces);
void find_forces_ghost(const gas_simulation *sim, phys_particle_t *gas, const int num_particles, phys_vector_t *forces, phys_particle_t *ghost, const int num_ghost );
void apply_relativity(const gas_simulation *sim, phys_particle_t *gas, const int num_particles);
void update_states(const gas_simulation *sim, const sim_scalar_t timestep, phys_particle_t *gas, const int num_particles, phys_vector_t *forces);


void marsaglia_alg(double *destination, int even_numtogen);
void maxwell_boltzmann(const gas_simulation *sim, const sim_scalar_t kT, phys_particle_t *gas, const int num_particles );


#ifdef __cplusplus
} /*extern C*/
#endif

#ifdef __cplusplus
#include <ostream>
inline std::ostream& operator <<(std::ostream& os, const phys_vector_t &object){
	return os << "{.x="<<object.x << ",.y="<<object.y<<"}";
}

inline std::ostream& operator <<(std::ostream& os, const phys_particle_t &object){
	return os << "{.p=" << object.p << ",.v="<< object.v<<"}";
}

inline std::ostream& operator <<(std::ostream& os, const bounding_box_t &object){
	return os << "{.min=" << object.min << ",.max="<< object.max<<"}";
}
#endif

#endif /* GASSIM_H */
