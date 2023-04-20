#include <math.h>

#include "gassim2d_util.h"


static const double sine_sixty = 0.8660254038;//sin(M_PI/3.0);
static const double sine_thirty = 0.5;//sin(M_PI/6.0);
static const double cosine_thirty = 0.8660254038; //= sine_sixty;
static const double regular_tetrahedron_height = 0.81649658092; //sqrt(2/3)

lattice_t create_lattice(const double pitch, const lattice_position_t min, const lattice_position_t max ){
	lattice_t our_lattice;
	our_lattice.pitch = pitch;
	our_lattice.min = min;
	our_lattice.max = max;

	double width = max.x - min.x;
	double height = max.y - min.y;

	our_lattice.N_x = (int)floor(width / pitch);
	our_lattice.N_y = (int)floor(height / pitch);
	
	//Ensure there is at least one row and one column
	if( our_lattice.N_x <= 0 ){ our_lattice.N_x = 1; }
	if( our_lattice.N_y <= 0 ){ our_lattice.N_y = 1; }

	our_lattice.N_total = our_lattice.N_x * our_lattice.N_y;

	return our_lattice;
}


lattice_position_t get_lattice_position(const lattice_t *a_lattice, const int i ){
	double p = a_lattice->pitch;

	int which_row = i/(a_lattice->N_x);
	int which_col = i%(a_lattice->N_x);

	return (lattice_position_t){ .x = a_lattice->min.x + p*which_col, .y = a_lattice->min.y + p*which_row };
}
