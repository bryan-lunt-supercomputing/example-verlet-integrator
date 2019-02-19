#include <math.h>

#include "gassim_util.h"


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
	double depth = max.z - min.z;

	our_lattice.N_x = (int)floor(width / pitch);
	if( our_lattice.N_x <= 0 ) our_lattice.N_x = 1;
	our_lattice.N_y = (int)floor(height / pitch);
	if( our_lattice.N_y <= 0 ) our_lattice.N_y = 1;
	our_lattice.N_z = (int)floor(depth / pitch);
	if( our_lattice.N_z <= 0 ) our_lattice.N_z = 1;

	our_lattice.N_total = our_lattice.N_x * our_lattice.N_y * our_lattice.N_z;

	return our_lattice;
}


lattice_position_t get_lattice_position(const lattice_t *a_lattice, const int i ){
	double p = a_lattice->pitch;
	int plane_size = a_lattice->N_x * a_lattice->N_y;

	int which_plane = i/(plane_size);
	int index_in_plane = i%plane_size;

	int which_row = index_in_plane/(a_lattice->N_x);
	int which_col = index_in_plane%(a_lattice->N_x);

	return (lattice_position_t){ .x = a_lattice->min.x + p*which_col, .y = a_lattice->min.y + p*which_row, .z = a_lattice->min.z + p*which_plane };
}
