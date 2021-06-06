// 02.05.2021 - Colloid in a fluctuating scalar field
// Davide Venturelli & Benjamin Walter
// Last update 01.06.2021
// Compile as $ gcc -o colloid colloid.c -lm -lgsl
/* COMMENTS:
- Stochastic Runge-Kutta II for colloid evolution, Euler-Maruyama for field evolution (can be enhanced, but the price is O(N) at least). We could even think of anisotropic resolution (better around the colloid).
- No boundary conditions on the colloid displacement; they only get enforced when locating the nearest site.
- Space is measured in units of the lattice spacing.
- Point-like colloid version.

TODO
- check SEED
- B: All I have done is MODEL B, I haven't worked on 'evolveA', hoping there'll be a way to merge evolveA and evolveB into a single function.

*/

// LIBRARIES, TYPES, DEFINITIONS

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>
#include <time.h>
#include <unistd.h>

#define ALLOC(p,n)  (p)=malloc( (n) * sizeof(*(p))); if( (p) == NULL){printf("Allocation of '%s' failed. Terminate. \n", #p); exit(2); } 
#define CALLOC(p,n)  (p)=calloc( (n) , sizeof(*(p))); if( (p) == NULL){printf("Allocation of '%s' failed. Terminate. \n", #p); exit(2); } 
//#define long unsigned long						 						// This is to make long hold twice as much
#define DD printf("# Debug: line %d \n",__LINE__);
#define DEBUG (1)
#define MOD 2															// 0 for Model A, 2 for Model B

// Data types
typedef struct{
	double mass;														// Mass of the field
	double lambda;														// Field-colloid coupling strength
	double quartic_u;													// Self-interaction coupling strength
	double temperature;													// Temperature of the common bath
	double relativeD;													// Ratio of colloid to field mobility
	double trap_strength;												// Stiffness of the harmonic trap
	int rng_seed;														// Seed for random number generator
	int system_size;													// Side of the DIM-dimensional lattice
	double delta_t;														// Time-discretization step
	long n_timestep;													// Number of timesteps
	int mc_runs;														// Number of Monte Carlo iterations
} parameter;

typedef struct{				// This structure bundles all observables 
	double** field_average;		// Saves the measured field-average <phi[i]> for certain subset of i's (eg along an axis) AND at all writing times
	double** field_correlation;	// Saves the measured field-correlator < phi[0] phi[i]> for certain subset of i's (eg along an axis) AND at all writing times
	double* colloid_msd;		// Saves the measured mean square displacement of the colloid at all writing times
	double write_time_delta;	// This is the gap between writing times (at the moment linear, maybe later exponential writing times?) 
	int write_count;
//	double* colloid_fpt_distribution;	// Saves the first-passage time distribution of the colloid
} observables;

// GSL RNG (it turns out that Ziggurat is faster than Box-Muller)
const gsl_rng_type *T;
gsl_rng *r;


// FUNCTION PROTOTYPES

void default_parameters(parameter*);
void initialise(double**, double**, long***, parameter*);
void initialise_observable(observables*, parameter*);
void wipe(double**, double**, parameter*);
void evolveB(double**, double**, long**, parameter*, observables*);
void evolveA(double**, double**, long**, parameter*, observables*);
void laplacian(double**, double*, long**, long);
void laplacian_of_cube(double**, double*, long**, long);
void generate_noise_field(double**, long, parameter*);
void gradient_field(double**, double*, long**, long);
void phi_evolveB(double**, double*, double*, double*, double*, long, parameter*, double*, long **);
void phi_evolveA(double**, double*, double*, long, parameter*, double*);
void measure(double**, double**, long, parameter*, observables*);
void print_observables(observables*, parameter*);
void print_params(parameter*);
void printhelp(void);
void print_source(void);
void neighborhood(long***, int);
int ind2j(int, int, int);
int vec2ind(int*, int);
int closest_site(double *, int);
unsigned modulo( int , unsigned );
unsigned int factorial(unsigned int);
int intpow(int, int);
double floatpow(double, int);

// Global variables (TODO?)
int DIM = 1;

// MAIN BODY

int main(int argc, char *argv[]){
	setlinebuf(stdout);
	
	// Input
	parameter params;
	default_parameters(&params);

	opterr = 0;
	int c = 0;
	// ./colloid -r 2 -L 23 (order doesn't count)
    while( (c = getopt (argc, argv, "L:r:l:u:T:d:k:S:t:N:M:D:hX") ) != -1){
    	switch(c){
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
			case 'd':
				params.relativeD = atof(optarg);
				break;
			case 'k':
				params.trap_strength = atof(optarg);
				break;
			case 'S':
				params.rng_seed = atoi(optarg);
				break;
			case 't':
				params.delta_t = atof(optarg);
				break;
			case 'N':
				params.n_timestep = atol(optarg);
				break;
			case 'M':
				params.mc_runs = atoi(optarg);
				break;
			case 'D':
				DIM = atoi(optarg);
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
	
	// Variables
	int i;
	long n_sites = intpow(params.system_size, DIM);
	
	double* phi; 														// Field on lattice NxNxN
	double* y_colloid; 													// Colloid position (DIM real numbers)
	long** neighbours; 													// i x j - table with j neighbours of site i
	
	ALLOC(phi, n_sites);
	CALLOC(y_colloid, DIM);
	ALLOC(neighbours, n_sites);

	//for(i = 0; i < n_sites; i++) ALLOC(neighbours[i], 2 * DIM); 		// At each index there are 2*D neighbours
	for(i = 0; i < n_sites; i++) {
		neighbours[i]=calloc( 2 * DIM , sizeof(long)); 
		if( neighbours[i] == NULL){printf("Allocation of A VERY POLITE ARRAY failed. Terminate. \n"); exit(2);}
	} 

	// Initialise fields
	initialise(&phi, &y_colloid, &neighbours, &params);

	print_params(&params);   // Print header with all parameters
	
	// Initialise observables
	observables obvs;	// Creates a pointer to an observables structure
	initialise_observable(&obvs, &params);

	// MC Iteration
	int mc_counter;
	for(mc_counter = 0; mc_counter < params.mc_runs; mc_counter++)
	{
		// Numerical integration of the dynamics
		//if(DEBUG){printf("#MCRUN %i\n", mc_counter);}
		if(MOD==0){															// Model A
			evolveA(&phi, &y_colloid, neighbours,&params, &obvs);
		}
		else{																// Model B
			evolveB(&phi, &y_colloid, neighbours, &params, &obvs);
		}

		wipe(&phi, &y_colloid, &params); // Reset to initial condition
	}

	print_observables(&obvs, &params);

	// Free memory and exit
	free(phi);
	free(y_colloid);
	for(i=0; i<n_sites; i++) free(neighbours[i]);
	return 0;
}


// DEFINITION OF FUNCTIONS

// Default parameters
void default_parameters(parameter* params){
	params->mass = 0.0;
	params->lambda = 0.1;
	params->quartic_u = 0.25;
	params->temperature = 1.0;
	params->relativeD = 1.0;
	params->trap_strength = 1.0;
	params->rng_seed = -1; 												// if seed is -1 (not given by user), it will be picked randomly in 'initialise'
	params->system_size = 32;
	params->delta_t = 0.001;
	params->n_timestep = 10;
	params->mc_runs = 10;
}

// All that needs to be done once
void initialise(double** phi, double** y_colloid, long*** neighbours, parameter* params){
	// i) GSL random number generator setup
	gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    if(params->rng_seed==-1) params->rng_seed = ((int) (((int) clock() ) % 100000)); // If no seed provided, draw a random one
    gsl_rng_set(r, params->rng_seed);
    printf("# RNG Seed %i\n",params->rng_seed);
	
	// ii) Initialise physical variables: phi, y
	long i, n_sites;																// We sample from an infinite temperature state. Maybe something else would be better.
	n_sites = intpow(params->system_size, DIM);
	for(i = 0; i < n_sites; i++){ 
		 (*phi)[i] = gsl_ran_gaussian_ziggurat(r,1.0); 
	}
	
	// y - colloid (initially in the middle of the lattice, where the harmonic well stands)
	long L = params->system_size;
	for(i=0; i<DIM; i++) (*y_colloid)[i] = ind2j(i, n_sites/2 ,L);
	
	// iii) What are each position's neighbours?
	neighborhood(neighbours,L);
}

// This initialises the observables structure later containing the measuremnrs
void initialise_observable(observables* obvs, parameter* params)
{
	obvs->write_time_delta = 0.2; // This is in physical time units, so writing occurs every (write_time_delta / n_timestep) integration step
	
	int writing_times = (int)(1 + (((params->n_timestep)*params->delta_t)/obvs->write_time_delta)); // This counts how many writing events will occur in time (including t=0, thus + 1)
	//if(DEBUG){printf("writing_times %i, n_timestep %lu, delta_t = %g, write_time_delta %g \n ", writing_times,params->n_timestep,  params->delta_t, obvs->write_time_delta);}
	CALLOC(obvs->colloid_msd, writing_times); // Save MSD vs time
	CALLOC(obvs->field_average, writing_times); // Prepare writing_times many arrays to store averages
	CALLOC(obvs->field_correlation, writing_times); // Prepare writing_times many arrays to store correlations
	int i;
	for(i = 0; i < writing_times; i++)
	{
		obvs->field_correlation[i] = calloc( params->system_size , sizeof(double) );
		if( (obvs->field_correlation[i]) == NULL){printf("Allocation of '(obvs->field_correlation[%i])'  failed. Terminate. \n", i); exit(2);} 
		obvs->field_average[i] = calloc( params->system_size , sizeof(double) );
		if( (obvs->field_average[i]) == NULL){printf("Allocation of '(obvs->field_average[%i])'  failed. Terminate. \n", i); exit(2);} 
	}
}

// Creates list of nearest neighbours in DIM dimensions
void neighborhood(long*** list, int L){
	int k,d;
	int vec[DIM], neigh[DIM];

	for(k=0; k<L*DIM; k++){												// Cycles over lattice sites
		for(d=0; d<DIM; d++){											// Finds the DIM coordinates of the current site
			vec[d] = ind2j(d,k,L);
			neigh[d] = vec[d];
		}
		for(d=0; d<DIM; d++){
			neigh[d] = (vec[d]+1)%L;									// Finds right neighbour
			(*list)[k][2*d] = vec2ind(neigh,L);							// Stores it in the neighbour list
			neigh[d] = (vec[d]+L-1)%L;									// Finds left neighbour
			(*list)[k][2*d+1] = vec2ind(neigh,L);						// Stores it in the neighbour list
			neigh[d] = vec[d];											// Restores local copy (prepares for next dimension)
		}
	}
}

// Reset field and colloid to init. cond.
void wipe(double** phi, double** y_colloid, parameter* params)
{
	int n_sites, i;
	n_sites = intpow(params->system_size, DIM);
	for(i = 0; i < n_sites; i++){ 
		 (*phi)[i] = gsl_ran_gaussian_ziggurat(r,1.0); 
	}
	
	// y - colloid (initially in the middle of the lattice, where the harmonic well stands)
	long L = params->system_size;
	for(i=0; i<DIM; i++) (*y_colloid)[i] = ind2j(i, n_sites/2 ,L);
}

// Model B evolution
void evolveB(double** phi, double** y_colloid, long** neighbours, parameter* params, observables* obvs){
	// Preparing variables for field evolution (Euler-Maruyama)
	long tstep, n_sites, n_timestep;
	n_sites = intpow(params->system_size, DIM);
	n_timestep = params->n_timestep;
	
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
	
	// Preparing variables for colloid evolution (Stochastic Runge-Kutta II)
	int i, y_site;
	double grad, noise, F2;															// There was also y_old but I think it's obsolete
	double w[DIM], F1[DIM], Y[3]={0};
	long L = params->system_size;
	double noise_intensity = sqrt(2.0 * params->temperature * params->relativeD * params->delta_t);
	for(i=0; i<DIM; i++) w[i] = ind2j(i,n_sites/2,L);								// Finds position of the harmonic well
	
	// Preparing measurement process
	int write_time_delta_step = (int) ((obvs->write_time_delta)/(params->delta_t));
	int next_writing_step = 0;
	obvs->write_count = 0;
	
	for(tstep = 0; tstep < n_timestep; tstep++){									// Time evolution
		
		// i) Create local copy of the colloid variable
		for(i=0; i<DIM; i++) Y[i] = (*y_colloid)[i];	
		
		// ii) Prediction step for the colloid
		y_site = closest_site(Y,L);													// Finds the lattice site closest to Y
		
		for(i=0; i<DIM; i++){ 														// Evolve each of the components separately
			// Compute gradient of the field under the colloid
			grad = 0.5 * ( (*phi)[y_site+intpow(L,i)] - (*phi)[y_site-intpow(L,i)] );
			
			// Compute temporary position
			noise = noise_intensity * gsl_ran_gaussian_ziggurat(r, 1.0);
			F1[i] =  params->relativeD * (params->lambda*grad -params->trap_strength*( Y[i]-w[i] ));
			(*y_colloid)[i] += params->delta_t * F1[i] + noise;
		}
		
		// iii) Evolve the field with the local copy of the colloid position	
		// I'm sure this can be optimised, but at least its explicit enough to see what's going on
		laplacian(&laplacian_phi, *phi, neighbours, n_sites);						// Write Laplacian of phi into laplacian_phi
		laplacian(&laplacian_square_phi, laplacian_phi, neighbours, n_sites);		// Write (D^2)^2 phi into laplacian_square_phi 
		laplacian_of_cube(&laplacian_phi_cubed, *phi, neighbours, n_sites);			// Write D2 [phi(x)^3] into laplacian_phi_cubed
		generate_noise_field(&noise_field, DIM*n_sites, params);					// Fill \vec{Lambda} with randomness
		gradient_field(&noise_gradient, noise_field, neighbours, n_sites);			// Compute gradient noise term
		

		// Add together to new step d/dt phi = -a * D2 phi - b D4 phi - u D2 (phi^3) + D * noise (D is nabla)
		printf("time: %g ", tstep*params->delta_t);
		phi_evolveB(phi, laplacian_phi, laplacian_square_phi, laplacian_phi_cubed, noise_gradient, n_sites, params, Y, neighbours);

		// iv) Correction step for the colloid (with the evolved field)
		y_site = closest_site(*y_colloid,L);										// Compute new site of the colloid
		
		for(i=0; i< DIM; i++){														// Evolve each of the components separately
			// Compute new gradient
			grad = 0.5 * ( (*phi)[y_site+intpow(L,i)] - (*phi)[y_site-intpow(L,i)] );
			
			// Compute corrected contribution
			F2 = params->relativeD * (params->lambda*grad -params->trap_strength*( (*y_colloid)[i]-w[i] ));
			
			// Sum the two contributions: y_n+1 = y_n + 1/2(F1+F2)*dt + noise
			noise = noise_intensity * gsl_ran_gaussian_ziggurat(r, 1.0);
			(*y_colloid)[i] = Y[i] + 0.5*(F1[i]+F2)*params->delta_t + noise;
		}

		// v) Measures
		if(tstep > next_writing_step)
		{
			//if(DEBUG){printf("Writing at t = %g, write count %i\n", tstep*params->delta_t, obvs->write_count);}
			measure(phi, y_colloid, tstep, params, obvs); 												// Evaluate all sorts of correlators etc. herea
			next_writing_step += write_time_delta_step;
			obvs->write_count++;
		}

	}
}

// Model A evolution
void evolveA(double** phi, double** y_colloid, long** neighbours, parameter* params, observables* obvs){
	// Preparing variables for field evolution (Euler-Maruyama)
	long tstep, n_sites, n_timestep;
	n_sites = intpow(params->system_size, DIM);
	n_timestep = params->n_timestep;
	
	double* laplacian_phi;
	ALLOC(laplacian_phi, n_sites);
	double* noise_field;
	ALLOC(noise_field, n_sites);

	// Preparing variables for colloid evolution (Stochastic Runge-Kutta II)
	int i, y_site;
	double grad, noise, F2;															// There was also y_old but I think it's obsolete
	double w[DIM], F1[DIM], Y[3]={0};
	long L = params->system_size;
	double noise_intensity = sqrt(2.0 * params->temperature * params->relativeD * params->delta_t);
	
	for(i=0; i<DIM; i++) w[i] = ind2j(i,n_sites/2,L);								// Finds position of the harmonic well
	
	for(tstep = 0; tstep < n_timestep; tstep++){									// Time evolution
									
		// i) Create local copy of the colloid variable
		for(i=0; i<DIM; i++) Y[i] = (*y_colloid)[i];	
		
		// ii) Prediction step for the colloid
		y_site = closest_site(Y,L);													// Finds the lattice site closest to Y
		
		for(i=0; i<DIM; i++){ 														// Evolve each of the components separately
			// Compute gradient of the field under the colloid
			grad = 0.5 * ( (*phi)[y_site+intpow(L,i)] - (*phi)[y_site-intpow(L,i)] );
			
			// Compute temporary position
			noise = noise_intensity * gsl_ran_gaussian_ziggurat(r, 1.0);
			F1[i] =  params->relativeD * (params->lambda*grad -params->trap_strength*( Y[i]-w[i] ));
			(*y_colloid)[i] += params->delta_t * F1[i] + noise;
		}
	
		// iii) Evolve the field with the local copy of the colloid position	
		laplacian(&laplacian_phi, *phi, neighbours, n_sites);						// Write Laplacian of phi into laplacian_phi
		generate_noise_field(&noise_field, n_sites, params);						// Fill Lambda with randomness
				
		// Add together to new step d/dt phi = -r * phi + \nabla^2 phi - u * (phi^3) - l * V(x-Y) + noise
		phi_evolveA(phi, laplacian_phi, noise_field, n_sites, params, Y);
		
		// iv) Correction step for the colloid (with the evolved field)
		y_site = closest_site(*y_colloid,L);										// Compute new site of the colloid
		 	
		for(i=0; i< DIM; i++){														// Evolve each of the components separately
			// Compute new gradient
			grad = 0.5 * ( (*phi)[y_site+intpow(L,i)] - (*phi)[y_site-intpow(L,i)] );
						
			// Compute corrected contribution
			F2 = params->relativeD * (params->lambda*grad -params->trap_strength*( (*y_colloid)[i]-w[i] ));
			
			// Sum the two contributions: y_n+1 = y_n + 1/2(F1+F2)*dt + noise
			noise = noise_intensity * gsl_ran_gaussian_ziggurat(r, 1.0);
			(*y_colloid)[i] = Y[i] + 0.5*(F1[i]+F2)*params->delta_t + noise;
		}

		// v) Measures
		measure(phi, y_colloid, tstep, params, obvs); 												// Evaluate all sorts of correlators etc. here
	}
}


// Field evolution - Model B
void phi_evolveB(double** phi, double* laplacian_phi, double* laplacian_square_phi, double* laplacian_phi_cubed, double* noise_gradient, long n_sites, parameter* params, double* Y, long ** neighbours){
	long i;
	long L = params->system_size;
	double delta_t = params->delta_t;
	printf("phi_old: %g, lapl_2_phi: %g, lapl_phi: %g, lapl_3_phi: %g, grad_noise: %g\n", (*phi)[0], laplacian_square_phi[0], laplacian_phi[0], laplacian_phi_cubed[0], noise_gradient[0]);
	for(i = 0; i < n_sites; i++){
		(*phi)[i] += (delta_t*( laplacian_square_phi[i] + params->mass*laplacian_phi[i] + params->quartic_u * laplacian_phi_cubed[i]) + noise_gradient[i]);
	}
	// Interaction with the colloid
	int y_site = closest_site(Y,L);													// Finds the lattice site closest to Y, the colloid
	(*phi)[y_site] += delta_t * params->lambda *2*DIM;
	for(i = 0; i < 2*DIM; i++){(*phi)[neighbours[y_site][i]] -= delta_t * params->lambda;}
}

// Field evolution - Model A
void phi_evolveA(double** phi, double* laplacian_phi, double* noise_field, long n_sites, parameter* params, double* Y){
	long i;
	double delta_t = params->delta_t;
	long L = params->system_size;
	for(i = 0; i < n_sites; i++)													// Notice noise_field already contains delta_t in its variance
	{
		(*phi)[i] += delta_t*(- params->mass*(*phi)[i] + laplacian_phi[i] + params->quartic_u * floatpow((*phi)[i],3)  ) + noise_field[i];
	}
	
	// Interaction with the colloid
	int y_site = closest_site(Y,L);													// Finds the lattice site closest to Y, the colloid
	(*phi)[y_site] += delta_t * params->lambda;										// Adds the interaction at that site
}

// Returns Laplacian as calculated from cubic neighbour cells in DIM dimensions
void laplacian(double** laplacian,  double* field, long** neighbours, long n_sites){
	double buffer;
	long pos, i;
	for(pos = 0; pos < n_sites; pos++){
		buffer = 0;
		for(i = 0; i < 2*DIM; i++) {buffer += field[neighbours[pos][i]]; }
		buffer -= (2*DIM*field[pos]);
		(*laplacian)[pos] = buffer;
	}
}

// Returns the Laplacian of field^3
void laplacian_of_cube(double** laplacian,  double* field, long** neighbours, long n_sites){
	double buffer;
	long pos, i;
	for(pos = 0; pos < n_sites; pos++){
		buffer = 0;
		for(i = 0; i < 2*DIM; i++){
			buffer += intpow( (field[neighbours[pos][i]]), 3); 
		}
		buffer -= (2*DIM*field[pos]*field[pos]*field[pos]);
		(*laplacian)[pos] = buffer;
	}
}

// This function generates a completely uncorrelated random field on a line
void generate_noise_field(double** noise_field, long length, parameter* params){
	long i;
	double noise_intensity = sqrt(2.0 * params->temperature * params->delta_t);
	for(i = 0; i < length; i++){
		(*noise_field)[i] = noise_intensity * gsl_ran_gaussian_ziggurat(r, 1.0);
	}
}

// Computes the gradient of a (noisy) field
void gradient_field(double** grad_noise, double* noise, long** neighbours, long n_sites){
	long i;
	int j; // neighbour of i
	double buffer;
	for(i = 0; i < n_sites; i++){
		buffer = 0;
		for(j = 0; j < DIM; j++){
			buffer += noise[neighbours[i][2*j]];
			buffer -= noise[neighbours[i][2*j+1]];
		}
		(*grad_noise)[i] = 0.5*buffer;
	}
}

// Perform measurements and print them out
void measure(double** phi, double** y_colloid, long tstep, parameter* params, observables* obvs){
	long i;
	//printf("%g\t", tstep*params->delta_t);
	for(i = 0; i < params->system_size; i++)
	{
		obvs->field_average[obvs->write_count][i] += (*phi)[i];
		obvs->field_correlation[obvs->write_count][i] += ((*phi)[0] * (*phi)[i]);
	}	
	
	long n_sites = intpow(params->system_size , DIM);
	for(i = 0; i < DIM; i++)
	{
		obvs->colloid_msd[obvs->write_count] += ((*y_colloid)[i] - ind2j(i, n_sites/2 , params->system_size))*((*y_colloid)[i] - ind2j(i, n_sites/2 , params->system_size));
	}
	//printf("\n");
}

// Translate from a list index (i) to lattice indicex j, and viceversa.
inline int ind2j(int j, int i, int L) {return ((int)(i/intpow(L,j)))%L ;}
int vec2ind(int *vec, int L){
	int i, res=0;
	for(i=0; i<DIM; i++) res += vec[i] * intpow(L,i);
	return res;
}

// Finds index of closest lattice site to the vector "vec" 
int closest_site(double *vec, int L){
	int i, site=0;
	for(i=0; i<DIM; i++)
	{	
		site += (int)(modulo( round(vec[i])  , L) * intpow(L,i)); // I'm using a special modulo function to avoid getting negative return values, eg -3 % 10 = -3, but modulo(-3,10) = 7.
	}
	return site;
}

// Factorial of a number
unsigned int factorial(unsigned int n){
    if (n == 0) return 1;
    return n * factorial(n - 1);
}

// The usual pow(a,b)=exp(log(a) * b) is slow AF
int intpow(int a, int b){
	int i, res=1;
	for(i=0; i<b; i++) res *= a;
	return res;
} 

// The usual pow(a,b)=exp(log(a) * b) is slow AF
double floatpow(double a, int b){
	int i;
	double res=1;
	for(i=0; i<b; i++) res *= a;
	return res;
} 

unsigned modulo( int value, unsigned m) {
// According to https://stackoverflow.com/questions/14997165/fastest-way-to-get-a-positive-modulo-in-c-c this is still fast
    int mod = value % (int)m;
    if (mod < 0) {
        mod += m;
    }
    return mod;
}

void print_observables(observables* obvs, parameter* params)
{
	int i,j;
	int system_size = params->system_size;
	int write_count = obvs->write_count;
	double weight = 1/((double) params->mc_runs);
	for(i = 0; i < write_count; i ++)
	{
		printf("#FIELDCORR %g", i * (obvs->write_time_delta));
		for(j = 0; j < system_size; j++)
		{
			printf("\t%g", weight*(obvs->field_correlation[i][j]- (obvs->field_average[i][0] * obvs->field_average[i][j]) ));
		}
		printf("\n");
	}

	// Output MSD of colloid
	int dim;
	double centre_square;	// This is y(0)^2, to be subtracted from the measured square displacement for centering.
	for(dim = 0; dim < DIM; dim++)
	{
		centre_square = intpow( ind2j(dim, ( intpow(system_size , DIM) )/2 , system_size) , 2);
	}

	for(i = 0; i < write_count; i++)
	{
		printf("# COLLOIDMSD %g", i * (obvs->write_time_delta));
		printf("\t%.3f",weight*(obvs->colloid_msd[i]));
		printf("\n");
	}
}

void print_params(parameter* params)
{
printf("# Parameters\n\
# MASS %g\n\
# LAMBDA (field-colloid-coupling) %g\n\
# U (quartic coupling) %g\n\
# TEMPERATURE %g\n\
# RELATIVE MOTILITY D %g\n\
# TRAP STRENGTH %g\n\
# RNG SEED %u\n\
# L (System Size) %i\n\
# DIM %i\n\
# DELTA T %g\n\
# TIMESTEPS %lu\n\
# MONTE CARLO RUNS %i\n",
		params->mass, params->lambda, params->quartic_u, params->temperature, params->relativeD, params->trap_strength, params->rng_seed, params->system_size, DIM, params->delta_t, params->n_timestep, params->mc_runs);
}


// Prints out the instructions
void printhelp(void){
	printf("# Colloid in Gaussian Field\n# '%s' built %s\n\
# Use with options flags\n\
# -L Length of three-dimensional lattice\n\
# -r Mass of Gaussian field\n\
# -l Voupling strength between colloid and field\n\
# -u Quartic coupling strength\n\
# -T Temperature of bath\n\
# -d Relative motility colloid/field\n\
# -k Strength of harmonic trap\n\
# -S Seed for RNG\n\
# -t Integration timestep\n\
# -D Dimension\n\
# -N Number of timesteps\n\
# -M Number of Monte Carlo samples\n\
# -h To see this helpscreen\n\
# -X Output source code\n",__FILE__,__DATE__);
}

// Prints out this whole code
void print_source(){
    printf("/* Source Code %s, created %s */\n",__FILE__,__DATE__);
    FILE *fp;
    int c;
   
    fp = fopen(__FILE__,"r");

    do {
     c = getc(fp); 
 	 putchar(c);
 	}
 	while(c != EOF);     
    fclose(fp);
}
