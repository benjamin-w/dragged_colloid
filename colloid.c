#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>
#include <time.h>
#include <unistd.h>

#define ALLOC(p,n)  (p)=malloc( (n) * sizeof(*(p))); if( (p) == NULL){printf("Allocation of '%s' failed. Terminate. \n", #p); exit(2); } 
#define CALLOC(p,n)  (p)=calloc( (n) , sizeof(*(p))); if( (p) == NULL){printf("Allocation of '%s' failed. Terminate. \n", #p); exit(2); } 
#define long unsigned long // This is to make long hold twice as much
#define DD printf("# Debug: line %d \n",__LINE__);

#define DIM (3) // dimension 


typedef struct parameter {
	double mass;
	double lambda;
	double quartic_u;
	double temperature;
	double relativeD;
	double trap_strength;
	int rng_seed;
	int system_size;
	double delta_t;
} parameter;


// GSL RNG
const gsl_rng_type *T;
gsl_rng *r;

void default_parameters(parameter*);
void initialise(double**, long*, long***, long, parameter*);
void evolve(double**, long*, long**, long, long, parameter* );
void laplacian(double**, double*, long**, long);
void laplacian_of_cube(double**, double*, long**, long);
void generate_noise_field(double**, long, parameter*);
void gradient_field(double**, double*, long**, long);
void phi_evolve(double**, double*, double*, double*, double*, long, parameter*);
void measure(double**, long, parameter*);
void printhelp(void);
void print_source(void);

int main(int argc, char *argv[])
{
	setlinebuf(stdout);

	// input
	parameter params;

	default_parameters(&params);

	opterr = 0;
	int c = 0;
        while( (c = getopt (argc, argv, "L:r:l:u:T:D:k:S:N:hX") ) != -1)
	{                switch(c)
                        {
				case 'L':
					params.system_size = atoi(optarg);
					break;
				case 'r':
					params.mass = atof(optarg);
					break;
				case 'l':
					params.lambda = atof(optarg);
					break;
				case 'u':
					params.quartic_u = atof(optarg);
					break;
				case 'T':
					params.temperature = atof(optarg);
					break;
				case 'D':
					params.relativeD = atof(optarg);
					break;
				case 'k':
					params.trap_strength = atof(optarg);
					break;
				case 'S':
					params.rng_seed = atoi(optarg);
					break;
				case 'X':
					print_source();
					exit(2);
				case 'h':
					printhelp();
					exit(2);
				default:
					printhelp();
                                exit(EXIT_FAILURE);
                        }
	}// getopt ends
	
	int i;
	long n_sites = powl(params.system_size, DIM);
	long n_timestep = 1000;

	
	
	double* phi; // field on lattice NxNxN
	long y_colloid; // colloid position encoded via lattice position
	long** neighbours; // i x j - table with j neighbours of site i
	

	ALLOC(phi, n_sites);
	ALLOC(neighbours, n_sites);

	for(i = 0; i < n_sites; i++)
	{
		CALLOC(neighbours[i], 2 * DIM); // At each index there are 2*D neighbours
	}
	
	initialise(&phi, &y_colloid, &neighbours, n_sites, &params);
	evolve(&phi, &y_colloid, neighbours, n_sites, n_timestep, &params);

	return 0;
}

void default_parameters(parameter* params)
{
	params->mass = 0.0;
	params->lambda = 0.1;
	params->quartic_u = 0.25;
	params->temperature = 1.0;
	params->relativeD = 1.0;
	params->trap_strength = 1.0;
	params->rng_seed = -1; // if seed is -1 (not given by user), it will be picked randomly in 'initialise'
	params->system_size = 32;
	params->delta_t = 0.0001;
}

void initialise(double** phi, long* y_colloid, long*** neighbours, long n_sites, parameter* params)
{
	// All that needs to be done once
	// GSL random number generator setup
	gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc (T);
        if(params->rng_seed==-1) params->rng_seed = ((int) (((int) clock() ) % 100000)); // If no seed provided, draw a random one
        gsl_rng_set(r, params->rng_seed);
        printf("# RNG Seed %i\n",params->rng_seed);
	
	// Initialise physical variables: phi, y
	long i;
	for(i = 0; i < n_sites; i++)
	{
		 (*phi)[i] = gsl_ran_gaussian_ziggurat(r,1.0); // We sample from an infinite temperature state. Maybe something else would be better.
	}

	// y - colloid
	*y_colloid = (n_sites/2);

	// What are each position's neighbours? (x,y,z) -> j via j = x + N y + N^2 z
	long pos = 0;
	int dir;

	long x, y, z;
	long system_size = params->system_size;
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
				for(dir = 0; dir < 2 * DIM; dir++)
				{
					(*neighbours)[pos][0] = (((x + 1) % system_size) + yshift + zshift);
					(*neighbours)[pos][1] = (((x - 1) % system_size) + yshift + zshift);

					(*neighbours)[pos][2] = (x + ((y+1) % system_size) * system_size + zshift);
					(*neighbours)[pos][3] = (x + ((y-1) % system_size) * system_size + zshift);

					(*neighbours)[pos][4] = (x + yshift + ((z+1) % system_size) * Nsquare);
					(*neighbours)[pos][5] = (x + yshift + ((z-1) % system_size) * Nsquare);
				}
				pos++;	
			}
		}
	}
	
	
}

void evolve(double** phi, long* y_colloid, long** neighbours, long n_sites, long n_timestep, parameter* params)
{
	// Model B
	long tstep;

	double* laplacian_phi;
	ALLOC(laplacian_phi, n_sites);

	double* laplacian_square_phi;
	ALLOC(laplacian_square_phi, n_sites);

	double* laplacian_phi_cubed;
	ALLOC(laplacian_phi_cubed, n_sites);

	double* noise_field;
	ALLOC(noise_field, DIM * n_sites);

	double* noise_gradient;
	ALLOC(noise_gradient, n_sites);	


	for(tstep = 0; tstep < n_timestep; tstep++)
	{
		// I'm sure this can be optimised, but at least its explicit enough to see what's going on
		// Evaluate Lattice Laplacian
		laplacian( &laplacian_phi, *phi, neighbours, n_sites);			// Write Laplacian of phi into laplacian_phi
		laplacian( &laplacian_square_phi, laplacian_phi, neighbours, n_sites);	// Write (D^2)^2 phi into laplacian_square_phi 
		laplacian_of_cube(&laplacian_phi_cubed, *phi, neighbours, n_sites);		// Write D2 [phi(x)^3] into laplacian_phi_cubed
		generate_noise_field(&noise_field, n_sites, params);		// Fill \vec{Lambda} with randomness
		
		gradient_field(&noise_gradient, noise_field, neighbours,  n_sites);		// Compute gradient noise term
		
		// Add together to new step d/dt phi = -a * D2 phi - b D4 phi - u D2 (phi^3) + D * noise 
		phi_evolve(phi, laplacian_phi, laplacian_square_phi, laplacian_phi_cubed, noise_gradient, n_sites, params); // Only phi-coupling
		
		//y_evolve(
		
		measure(phi, tstep, params); // Evaluate all sorts of correlators etc. here
	}
	
}

void laplacian(double** laplacian,  double* field, long** neighbours,  long n_sites)
{
	// Returns Laplacian as calculated from six cubic neighbour cells
	int buffer;
	long pos, i;
	for(pos = 0; pos < n_sites; pos++)
	{
		buffer = 0;
		for(i = 0; i < 2*DIM; i++)
		{
			buffer += field[neighbours[pos][i]]; 
		}
		buffer -= (2*DIM*field[pos]);
		(*laplacian)[pos] = buffer;
	}
}

void laplacian_of_cube(double** laplacian,  double* field, long** neighbours,  long n_sites)
{
	// Returns the Laplacian of field^3
	int buffer;
	long pos, i;
	for(pos = 0; pos < n_sites; pos++)
	{
		buffer = 0;
		for(i = 0; i < 2*DIM; i++)
		{
			buffer += (field[neighbours[pos][i]])*(field[neighbours[pos][i]])*(field[neighbours[pos][i]]); 
		}
		buffer -= (2*DIM*field[pos]*field[pos]*field[pos]);
		(*laplacian)[pos] = buffer;
	}
}

void generate_noise_field(double** noise_field, long n_sites, parameter* params)
{
	// This function generates a completely uncorrelated 3-dim. random field on the lattice
	long i;
	double noise_intensity = sqrt(2.0 * params->temperature * params->delta_t);
	for(i = 0; i < (DIM * n_sites); i++)
	{
		(*noise_field)[i] = noise_intensity * gsl_ran_gaussian_ziggurat(r, 1.0);
	}
}

void gradient_field(double** grad_noise, double* noise, long** neighbours, long n_sites)
{
	// Computes the gradient of a (noisy) field
	long i;
	int j; // neighbour of i
	double buffer;
	for(i = 0; i < n_sites; i++)
	{
		buffer = 0;
		for(j = 0; j < DIM; j++)
		{
			buffer += noise[neighbours[i][2*j]];
			buffer -= noise[neighbours[i][2*j+1]];
		}
		(*grad_noise)[i] = 0.5*buffer;
	}
}


void phi_evolve(double** phi, double* laplacian_phi, double* laplacian_square_phi, double* laplacian_phi_cubed, double* noise_gradient, long n_sites, parameter* params)
{
	// Physics happens
	long i;
	double delta_t = params->delta_t;
	for(i = 0; i < n_sites; i++)
	{
		(*phi)[i] += (delta_t*(-laplacian_square_phi[i] - params->mass*laplacian_phi[i] - params->quartic_u * laplacian_phi_cubed[i]) + noise_gradient[i]);
	}

}
/*
void deep_swap();
*/
void measure(double** phi, long tstep, parameter* params)
{
	long i;
	printf("%g\t", tstep*params->delta_t);
	for(i = 0; i < params->system_size; i++)
	{
		printf("%g\t",((*phi)[0]) * (*phi)[i] );
	}
	printf("\n");
}

void printhelp(void)
{
	printf("# Colloid in Gaussian Field\n# '%s' built %s\n\
# Use with options flags\n\
# -L Length of three-dimensional lattice\n\
# -r Mass of Gaussian field\n\
# -l Voupling strength between colloid and field\n\
# -u Quartic coupling strength\n\
# -T Temperature of bath\n\
# -D Relative motility colloid/field\n\
# -k Strength of harmonic trap\n\
# -S Seed for RNG\n\
# -h To see this helpscreen\n\
# -X Output source code\n",__FILE__,__DATE__);
}

void print_source()
{
    printf("/* Source Code %s, created %s */\n",__FILE__,__DATE__);
    FILE *fp;
    int c;
   
    fp = fopen(__FILE__,"r");

    do {
         c = getc(fp); 
 	 putchar(c);    }
    while(c != EOF);     
    fclose(fp);
}
