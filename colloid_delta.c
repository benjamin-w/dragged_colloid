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
#define DIM (3) 														// Dimension 
#define MOD 0															// 0 for Model A, 2 for Model B

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
} parameter;

// GSL RNG (it turns out that Ziggurat is faster than Box-Muller)
const gsl_rng_type *T;
gsl_rng *r;


// FUNCTION PROTOTYPES

void default_parameters(parameter*);
void initialise(double**, double**, long***, parameter*);
void evolveB(double**, double**, long**, parameter*);
void evolveA(double**, double**, long**, parameter*);
void laplacian(double**, double*, long**, long);
void laplacian_of_cube(double**, double*, long**, long);
void generate_noise_field(double**, long, parameter*);
void gradient_field(double**, double*, long**, long);
void phi_evolveB(double**, double*, double*, double*, double*, long, parameter*, double*, long **);
void phi_evolveA(double**, double*, double*, long, parameter*, double*);
void measure(double**, long, parameter*);
void printhelp(void);
void print_source(void);
void neighborhood(long***, int);
int ind2j(int, int, int);
int vec2ind(int*, int);
int closest_site(double *, int);
unsigned int factorial(unsigned int);
int intpow(int, int);
double floatpow(double, int);


// MAIN BODY

int main(int argc, char *argv[]){
	setlinebuf(stdout);
	printf("# RNG Seed %i\n",3%10);
	printf("# RNG Seed %i\n",(-3)%10);

	// Input
	parameter params;
	default_parameters(&params);

	opterr = 0;
	int c = 0;
	// ./colloid -r 2 -L 23 (order doesn't count)
    while( (c = getopt (argc, argv, "L:r:l:u:T:D:k:S:N:hX") ) != -1){
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
			case 'D':
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
		if( neighbours[i] == NULL){printf("Allocation of your mum failed. Terminate. \n"); exit(2);}
	} 

	// Numerical integration of the dynamics
	
	initialise(&phi, &y_colloid, &neighbours, &params);
	
	if(MOD==0){															// Model A
		evolveA(&phi, &y_colloid, neighbours,&params);
	}
	else{																// Model B
		evolveB(&phi, &y_colloid, neighbours, &params);
	}
		
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
	params->delta_t = 0.0001;
	params->n_timestep = 10;
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
	if(DEBUG)
	{
		for(k=0; k < 2*DIM; k++)
		{
			printf("# CHECKNEIGHBOUR \t%li\t%li\n", (*list)[0][k], (*list)[L*DIM-1][k]);
		}
	}
}

// Model B evolution
void evolveB(double** phi, double** y_colloid, long** neighbours, parameter* params){
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
		measure(phi, tstep, params); 												// Evaluate all sorts of correlators etc. here
	}
}

// Model A evolution
void evolveA(double** phi, double** y_colloid, long** neighbours, parameter* params){
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
		measure(phi, tstep, params); 												// Evaluate all sorts of correlators etc. here
	}
}

// Field evolution - Model B
void phi_evolveB(double** phi, double* laplacian_phi, double* laplacian_square_phi, double* laplacian_phi_cubed, double* noise_gradient, long n_sites, parameter* params, double* Y, long ** neighbours){
	long i;
	long L = params->system_size;
	double delta_t = params->delta_t;
	for(i = 0; i < n_sites; i++){
		(*phi)[i] += (delta_t*(-laplacian_square_phi[i] - params->mass*laplacian_phi[i] - params->quartic_u * laplacian_phi_cubed[i]) + noise_gradient[i]);
	}
	
	// Interaction with the colloid
	int y_site = closest_site(Y,L);													// Finds the lattice site closest to Y, the colloid
	(*phi)[y_site] += delta_t * params->lambda *2*DIM;
	for(i = 0; i < 2*DIM; i++) (*phi)[neighbours[y_site][i]] -= delta_t * params->lambda;
}

// Field evolution - Model A
void phi_evolveA(double** phi, double* laplacian_phi, double* noise_field, long n_sites, parameter* params, double* Y){
	long i;
	double delta_t = params->delta_t;
	long L = params->system_size;
	for(i = 0; i < n_sites; i++)													// Notice noise_field already contains delta_t in its variance
	{
		(*phi)[i] += delta_t*(- params->mass*(*phi)[i] + laplacian_phi[i] - params->quartic_u * floatpow((*phi)[i],3)  ) + noise_field[i];
	}
	
	// Interaction with the colloid
	int y_site = closest_site(Y,L);													// Finds the lattice site closest to Y, the colloid
	(*phi)[y_site] += delta_t * params->lambda;										// Adds the interaction at that site
}

// Returns Laplacian as calculated from cubic neighbour cells in DIM dimensions
void laplacian(double** laplacian,  double* field, long** neighbours, long n_sites){
	int buffer;
	long pos, i;
	for(pos = 0; pos < n_sites; pos++){
		buffer = 0;
		for(i = 0; i < 2*DIM; i++) buffer += field[neighbours[pos][i]]; 
		buffer -= (2*DIM*field[pos]);
		(*laplacian)[pos] = buffer;
	}
}

// Returns the Laplacian of field^3
void laplacian_of_cube(double** laplacian,  double* field, long** neighbours, long n_sites){
	int buffer;
	long pos, i;
	for(pos = 0; pos < n_sites; pos++){
		buffer = 0;
		for(i = 0; i < 2*DIM; i++){
			buffer += (field[neighbours[pos][i]])*(field[neighbours[pos][i]])*(field[neighbours[pos][i]]); 
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
void measure(double** phi, long tstep, parameter* params){
	long i;
	printf("%g\t", tstep*params->delta_t);
	for(i = 0; i < params->system_size; i++) printf("%g\t",((*phi)[0]) * (*phi)[i] );
	printf("\n");
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
	for(i=0; i<DIM; i++) site += (int)round(vec[i])%L * intpow(L,i);
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

// Prints out the instructions
void printhelp(void){
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
# -t Integration timestep\n\
# -N Number of timesteps\n\
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
