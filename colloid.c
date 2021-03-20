#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>
#include <time.h>

#define ALLOC(p,n)  (p)=malloc( (n) * sizeof(*(p))); if( (p) == NULL){printf("Allocation of '%s' failed. Terminate. \n", #p); exit(2); } 
#define CALLOC(p,n)  (p)=calloc( (n) , sizeof(*(p))); if( (p) == NULL){printf("Allocation of '%s' failed. Terminate. \n", #p); exit(2); } 
#define long unsigned long // This is to make long hold twice as much
#define DD printf("# Debug: line %d \n",__LINE__);

int system_size = 4;
int dim=3;

double* phi; // field
double* y_colloid; // colloid position
long** neighbours; // table with neighbours 	
	


// GSL RNG
const gsl_rng_type *T;
gsl_rng *r;
int seed;

void initialise(long);
void evolve(int, long, long );

int main(int argc, char *argv[])
{
	int i;
	long n_sites = powl(system_size, dim);
	long n_timestep = 100;
	
	ALLOC(phi, n_sites);
	ALLOC(y_colloid, dim);
	ALLOC(neighbours, n_sites);
	for(i = 0; i < n_sites; i++)
	{
		CALLOC(neighbours[i], 2 * dim); // At each index there are 2*D neighbours
	}

	initialise(n_sites);

	evolve(system_size, n_sites, n_timestep);

	return 0;
}

void initialise(long n_sites)
{
	// All that needs to be done once
	// GSL random number generator setup
	gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc (T);
        if(seed==-1) seed = ((int) (((int) clock() ) % 100000)); // If no seed provided, draw a random one
        gsl_rng_set(r, seed);
        printf("# RNG Seed %i\n",seed);
	
	// Initialise physical variables: phi, y
	long i;
	for(i = 0; i < n_sites; i++)
	{
		 *(phi + i) = gsl_ran_gaussian_ziggurat(r,1.0); // We sample from an infinite temperature state. Maybe something else would be better.
	}

	for(i = 0; i < dim; i++)
	{
		*(y_colloid + i) = 0.0;
	}

	// What are each position's neighbours? (x,y,z) -> j via j = x + N y + N^2 z
	long pos = 0;
	int dir;

	long x, y, z;
	long Nsquare = system_size * system_size;
	long zshift, yshift;

	
	for(z = 0; z < system_size; z++)
	{
		zshift = (z * Nsquare); // Add this to get pos in this z-plane
		for(y = 0; y < system_size; y++)
		{
			yshift = (y * system_size); // Add this to get pos in this y-plane
			for( x = 0; x < system_size; x++)
			{
				//printf("%lu %lu %lu -> %lu \n", x, y, z, pos);
				for(dir = 0; dir < 2 * dim; dir++)
				{
					neighbours[pos][0] = (((x + 1) % system_size) + yshift + zshift);
					neighbours[pos][1] = (((x - 1) % system_size) + yshift + zshift);

					neighbours[pos][2] = (x + ((y+1) % system_size) * system_size + zshift);
					neighbours[pos][3] = (x + ((y-1) % system_size) * system_size + zshift);

					neighbours[pos][4] = (x + yshift + ((z+1) % system_size) * Nsquare);
					neighbours[pos][5] = (x + yshift + ((z-1) % system_size) * Nsquare);
				}
				pos++;	
			}
		}
	}
	

}

void evolve(int system_size, long n_sites, long n_timestep)
{
	// Model B
	long tstep;

	int q = pow(2, dim); // coordination number
	long N1 = n_sites;
	long N2 = n_sites*n_sites;
	
	double* new_phi;
	ALLOC(new_phi, n_sites);
	
	double* laplacian_phi;
	ALLOC(laplacian_phi, n_sites);

	double* laplacian_square_phi;
	ALLOC(laplacian_square_phi, n_sites);

	double* laplacian_phi_cubed;
	ALLOC(laplacian_phi_cubed, n_sites);

	double* noise_field;
	ALLOC(noise_field, dim * n_sites);

	double* noise_gradient;
	ALLOC(noise_gradient, n_sites);	

	double 
	double buffer;
	long pos;
	
	int i;

	for(tstep = 0; tstep < n_timestep; tstep++)
	{
		// Evaluate Lattice Laplacian
		laplacian( &laplacian_phi, phi, n_sites);			// Write Laplacian of phi into laplacian_phi
		laplacian( &laplacian_square_phi, laplacian_phi, n_sites);	// Write (D^2)^2 phi into laplacian_square_phi 
		laplacian_of_cube( &laplacian_phi_cubed, phi, n_sites);		// Write D2 [phi(x)^3] into laplacian_phi_cubed
		generate_noise_field(&noise_field, n_sites, dim);		// Fill \vec{Lambda} with randomness
		gradient_field(&noise_gradient, noise_field, n_sites);		// Compute gradient noise term
	
		// Add together to new step d/dt phi = -a * D2 phi - b D4 phi - u D2 (phi^3) + D * noise 
		phi_time_step(&new_phi, laplacian_phi, laplacian_square_phi, laplacian_phi_cubed, noise_gradient, parameters); 
		swap(&new_phi, &phi);
		measure(phi); // Evaluate all sorts of correlators etc. here
	}
}

void laplacian(double** output_ptr,  double* ptr, long n_sites)
{
	int buffer = 0.0;
	for(i = 0; i < 2*dim; i++)
	{
		buffer += ptr[neighbours[pos][i]]; 
	}
	buffer -= (2*dim*phi[pos]);
	return buffer;
}

double laplacian_square( double* ptr, long n_sites, long pos);
{
	// D4 f(x) = D2f(x+h) + D2f(x-h) - 2 D2f(x)
	
}
