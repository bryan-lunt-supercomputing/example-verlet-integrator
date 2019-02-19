#ifndef GASSIM_UTIL_H
#define GASSIM_UTIL_H

#include <math.h>


#ifdef __cplusplus
extern "C" {
#endif

typedef struct { double x; double y; } lattice_position_t;
typedef struct { double pitch;
				lattice_position_t min; lattice_position_t max;
				int N_x; int N_y;
				int N_total;
} lattice_t;

lattice_t create_lattice(const double pitch, const lattice_position_t min, const lattice_position_t max );
lattice_position_t get_lattice_position(const lattice_t *a_lattice, const int i );


#ifdef __cplusplus
} /*extern C*/
#endif


#endif /*GASSIM_UTIL_H*/
