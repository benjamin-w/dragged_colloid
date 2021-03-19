#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>
#include <time.h>

#define long unsigned long
int system_size = 16;
int dim=3;

typedef struct bdry_stack
{
	long position;
	char bdry_type;
} bdry_stack;

// GSL RNG
const gsl_rng_type *T;
gsl_rng *r;
int seed;

void initialise(double*, double*, long);
void evolve(double*, double*, int, long, long );

int main(int argc, char *argv[])
{
	double* phi;
	double* y;
	bdry_stack* boundary;
	
	long n_sites = powl(system_size, dim);
	long n_timestep = 100;
	
	phi = malloc(n_sites * sizeof(double));
	y = malloc(dim * sizeof(double));
	boundary = malloc( (n_sites - powl((system_size - 1), dim)) * sizeof(bdry_stack));
	initialise(y, phi, n_sites);

	evolve(y, phi, system_size, n_sites, n_timestep);

	return 0;
}

void initialise(double* y_colloid,double* phi, bdry_stack* boundary, long n_sites)
{
	gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc (T);
        if(seed==-1) seed = ((int) (((int) clock() ) % 100000));
        gsl_rng_set(r, seed);
        printf("# RNG Seed %i\n",seed);
	
	long i;
	for(i = 0; i < n_sites; i++)
	{
		phi[i] = gsl_ran_gaussian_ziggurat(r,1.0);
	}

	for(i = 0; i < dim; i++)
	{
		y_colloid[i] = 0.0;
	}

	long N1 = system_size;
	long N2 = system_size*system_sizel
	long x,y,z;
	long bdry_counter = 0;
	// Three bits count all possible bdry points in three dimensions: 
	for(i = 0; i < n_sites; i++)
	{
		z = (i / N2);
		y = ((i - z*N2) / N1);
		x = (i - z*N2 - y*N1);
		if(x == 0 || x == (system_size - 1))
		{
			boundary[bdry_counter].position = i;
			boundary[bdry_counter].bdry_type += 1;
		}
		
		if(y == 0 || y == (system_size - 1))
		{
			boundary[bdry_counter].position = i;
			boundary[bdry_counter].bdry_type += 2;
		}

		if(z == 0 || z == (system_size - 1))
		{
			boundary[bdry_counter].position = i;
			boundary[bdry_counter].bdry_type += 4;
		}

	}	

}

void evolve(double* y_colloid, double* phi, int system_size, long n_sites, long n_timestep)
{
	// Model B
	long tstep;

	int q = pow(2, dim); // coordination number
	long N1 = n_sites;
	long N2 = n_sites*n_sites;
	
	double* new_phi = malloc(n_sites * sizeof(double));
	long* neighbours = malloc(q * sizeof(long));
	double buffer;
	long pos;
	
	int i;

	for(tstep = 0; tstep < n_timestep; tstep++)
	{
		// Evaluate Lattice Laplacian
		for(pos = 0 ; pos < n_sites; pos++)
		{
			for(i = 0; i < q; i++)
			{	
				neighbours[i] = 0;
			}
		}

	}
}

// TODO: neighbours
