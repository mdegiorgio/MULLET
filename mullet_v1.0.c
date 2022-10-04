/**********************************
 ********** BALLETT v1.0 ***********
 **********************************

This software requires that you have GNU Scientific Library (GSL) installed
This software has only been tested on Linux

For speed of computation the program cycles over a predefined set of equilbrium values and recombination rates
The probality matrix for probabilities and substitutions is pre-computed prior to the start of any analysis
The probability matrix for frequency spectra is pre-computed using simulations

*/

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

double *rate_grid = NULL; //holds vector of grid of scaled recombination rates
double *equilibria_grid = NULL; //holds vector of grid of equilibria values
int num_sites = 0;
int num_equilibria = 0;
double e_min = 0.05; // Minimum equilibrium value
double e_max = 0.951; // Maximimum equilibrium value
double e_inc = 0.05; // Equilibrium increment
int num_rates = 0;

int num_species = 0; // Number of species K, as input by the user
int **topology = NULL;
double *branches = NULL;
int *n_min = NULL;
int *n_max = NULL;
double recratemin = 0.0;
double recratemax = 0.0;
//Variable stores probabilities of substitution and polymorphism over a grid of recombination values
//prob_poly_sub_trans[n1,...,nK][0][x][r] holds probability of polymorphism only in pop 1 in sample sizes n1, ..., nK that is recombination distances r at equilibrium frequency x
// ...
//prob_poly_sub_trans[n1,...,nK][K-1][x][r] holds probability of polymorphism only in pop K in sample sizes n1, ..., nK that is recombination distances r at equilibrium frequency x
//prob_poly_sub_trans[n1,...,nK][K][x][r] holds probability of subsitution with split 1 in sample sizes n1, ..., nK that is recombination distances r at equilibrium frequency x
// ...
//prob_poly_sub_trans[n1,...,nK][3K-4][x][r] holds probability of subsitution with split 2K-3 in sample sizes n1, ..., nK that is recombination distances r at equilibrium frequency x
//Variable is global as we will use this throughout whole program
//Variable is accessed using via prob_poly_sub_trans[N][c][x][r], where
//  N = unique integer from some function of N = (n[0]-n_min[0])*(n_max[0]-n_min[0]+1)^(K-1) + (n[1]-n_min[1])*(n_max[1]-n_min[1]+1)^(K-2) + ... + (n[K-1]-n_min[K-1])*(n_max[K-1]-n_min[K-1]+1)^(K-K)
//  c = sample configuration, taking on values 0 to 3K-4
//  x = equilibrium frequency
//  r = integer for recombination bin
double ****prob_poly_sub_trans = NULL;

double *log_fact = NULL;
double **log_binom = NULL;

double ***empiricalSpectra; // For T2 Transspecies and species 1, 2, ..., K
double *****simulated_spectrum = NULL; // For species 1, 2, ..., K

struct datatype {
	double loc;
    double rate;
	int *x;
	int *n;
    int config; // Index of configuration of data
    int N; // Integer indicated the sample index
} *data;

struct datatype *data=NULL;

void InitFact(void); 
void InitBinom(void);
//void InitFact_sample_size(int); 
//void InitBinom_sample_size(int);
double h(int, int, double, double, double);
double coal1(int, int, double, double, double);
double coal2(int, int, double, double, double);
double mig1(int, int, double, double, double);
double mig2(int, int, double, double, double);
gsl_vector* getDiagonal(int);
gsl_vector* getSubDiagonal(int, double, double, double);
gsl_vector* getSuperDiagonal(int, double, double, double);
double lnPower(double, double);
gsl_vector* buildConstHeightVector(int, double, double, double, gsl_vector*);
double expectedTreeHeightBalancing(int, double, double, double, double);
gsl_vector* buildConstLengthVector(int, double, double, double, gsl_vector*);
double expectedTreeLengthBalancing(int, double, double, double, double);
//void ProbPolyTrans(double*, double*, double, double, double, double, double, double, double);
void InitPolySubMatrixTrans(double, double);
double calc_ln_like_balsel_trans(int, int, int, double, double);
double calc_ln_like_background_trans(int, int, double**);
void scan_with_T1_trans(char*, int,double, double, char*);
double calc_ln_like_balsel_spectrum_trans(int, int, int, double, double);
double calc_ln_like_background_spectrum_trans(int, int, double**);
void scan_with_T2_trans(char*, int,double, double, char*);
void LoadMarginalSpectra(char*);
void get_min_max_sample_sizes(void);
void get_min_max_rho(void);
void read_tree_file(char*);
void assignConfig(void);
void readsnps(char*);
double** read_trans_configs(char*);
double **calcTransConfigFreq(void);
void printFractionTransConfig(char*, double**);
void freeTransConfig(double**);
//void ComputeSpectrum(char*, int, double, double, double, int);
void read_marginal_empirical_spectra(char**);
void obtain_marginal_empirical_spectra(void);
void write_empirical_spectra(char**, double***);
void free_empirical_spectra(double***);
void usage(void);

int main(int argc, char *argv[])
{
	srand(time(NULL));

	int gridsize;
	int sample_size;
	int num_reps;
	int totNumSites = 0;
	double equilibrium_frequency;
	double rho, theta1 = 0.05, theta2 = 0.05;
	char snpfn[1000], outfn[1000], recmapfn[1000], divfn[1000], polysubfn[1000], path[10000], transconfigfn[1000], treefn[1000];
    char **spectfn;
	double divergence;
	double seqDiff;
	double theta;
	double *polyFreq;
	double *subsFreq;
	double **transConfigFreqs;
    int i, j;

	if(argc < 3) {
		usage();
	}
		
    if(strcmp(argv[1], "-config") == 0) {
		if(argc!=6) {
			usage();
		}
		else {
			printf("You have chosen to get proportion of sites with a given transspecies configuration\n");
		}

        num_species = atoi(argv[2]);

        sprintf(treefn, "%s", argv[3]);
        read_tree_file(treefn);

		sprintf(snpfn, "%s", argv[4]);
        readsnps(snpfn);

        assignConfig();

		transConfigFreqs = calcTransConfigFreq();
		
        sprintf(outfn, "%s", argv[5]);
		printFractionTransConfig(outfn, transConfigFreqs);

		freeTransConfig(transConfigFreqs);
    }
    else if(strcmp(argv[1], "-spect") == 0) {
		if(argc<5) {
			usage();
		}
		
        num_species = atoi(argv[2]);

        if(argc!=(4+num_species)) {
			usage();
		}
        else {
			printf("You have chosen to get marginal empirical spectra\n");
		}

		sprintf(snpfn, "%s", argv[3]);
    
        spectfn = (char**)malloc(num_species * sizeof(char*));        
        for(i = 0; i < num_species; i++) {
            spectfn[i] = (char*)malloc(1000*sizeof(char));
            sprintf(spectfn[i], "%s", argv[4+i]);        
        }

		readsnps(snpfn);
    
        obtain_marginal_empirical_spectra();
        write_empirical_spectra(spectfn, empiricalSpectra);
        free_empirical_spectra(empiricalSpectra);
        
        for(i = 0; i < num_species; i++) {
            free(spectfn[i]);
        }
        free(spectfn);
	}
    else if(strcmp(argv[1], "-T1trans") == 0) {
		if(argc!=8) {
			usage();
		}
		else {
			printf("You have chosen to perform a scan using T1trans for transspecies polymorphisms\n");
		}

		gridsize = atoi(argv[2]);

		if(gridsize <=0) {
			printf("windowsize should be > 0!\n");
			usage();
		}

        num_species = atoi(argv[3]);

        sprintf(treefn, "%s", argv[4]);
        read_tree_file(treefn);

        sprintf(transconfigfn, "%s", argv[5]);

        sprintf(snpfn, "%s", argv[6]);
        readsnps(snpfn);

        assignConfig();        
		
		sprintf(outfn, "%s", argv[7]);
				
		InitFact();
		InitBinom();
		InitPolySubMatrixTrans(theta1, theta2);

		scan_with_T1_trans(outfn, gridsize, theta1, theta2, transconfigfn);
	}
    else if(strcmp(argv[1], "-T2trans") == 0) {
		if(argc<10) {
			usage();
		}
	
        num_species = atoi(argv[3]);

        if(argc!=(9+num_species)) {
			usage();
		}
        else {
			printf("You have chosen to perform a scan using T2trans for transspecies polymorphisms\n");
		}

		gridsize = atoi(argv[2]);

		if(gridsize <=0) {
			printf("windowsize should be > 0!\n");
			usage();
		}

        sprintf(treefn, "%s", argv[4]);
        read_tree_file(treefn);

        sprintf(transconfigfn, "%s", argv[5]);

        spectfn = (char**)malloc(num_species * sizeof(char*));        
        for(i = 0; i < num_species; i++) {
            spectfn[i] = (char*)malloc(1000*sizeof(char));
            sprintf(spectfn[i], "%s", argv[6+i]);        
        }

        sprintf(snpfn, "%s", argv[6+num_species]);
        readsnps(snpfn);

        assignConfig();        
		
        sprintf(path, "%s", argv[7+num_species]);
		sprintf(outfn, "%s", argv[8+num_species]);

        read_marginal_empirical_spectra(spectfn);
				
		InitFact();
		InitBinom();
		InitPolySubMatrixTrans(theta1, theta2);

        LoadMarginalSpectra(path);

		scan_with_T2_trans(outfn, gridsize, theta1, theta2, transconfigfn);
	}
	else if(strcmp(argv[1], "-SimSpect") == 0) {
		if(argc!=8) {
			usage();
		}
		else {
			printf("You have chosen to simulate a frequency spectrum\n");
		}

		sample_size = atoi(argv[2]);
		equilibrium_frequency = atof(argv[3]);
		theta1 = atof(argv[4]);
		theta2 = atof(argv[5]);
		num_reps = atoi(argv[6]);
		sprintf(outfn, "%s", argv[7]);
		
	//	InitFact_sample_size(sample_size);
	///	InitBinom_sample_size(sample_size);
	//	ComputeSpectrum(outfn, sample_size, equilibrium_frequency, theta1, theta2, num_reps);
	}
	else {
		usage();
	}

	return 0;
}

void InitFact(void) 
{	
	int n = 0;
    int n_max_all = 0;

    for(n = 0; n < num_species; n++) {
        if(n_max_all < n_max[n]) {
            n_max_all = n_max[n];
        }
    }
	
	log_fact = (double*)malloc((n_max_all+1)*sizeof(double));
	log_fact[0] = 0.0;
	for(n = 1; n <= n_max_all; n++) {
		log_fact[n] = log_fact[n-1] + log((double)n);
	}
}

void InitBinom(void)
{
	int n = 0;
	int k = 0;
    int n_min_all = n_min[0];
    int n_max_all = 0;

    for(n = 0; n < num_species; n++) {
        if(n_max_all < n_max[n]) {
            n_max_all = n_max[n];
        }

        if(n_min[n] < n_min_all) {
            n_min_all = n_min[n];
        }
    }
 
	printf("Initializing binomial coefficients\n");

	log_binom = (double**)malloc((n_max_all+1)*sizeof(double*));
	for(n = n_min_all; n <= n_max_all; n++) {
		log_binom[n] = (double*)malloc((n_max_all+1)*sizeof(double));
		
		for(k = 0; k <= n; k++) {
			log_binom[n][k] = log_fact[n] - log_fact[n-k] - log_fact[k];
		}
	}
}

double h(int i, int j, double x, double b1, double b2)
{
	if((i <= 1) && (j <= 1)) {
		return ((j*b1*x)/(1-x)) + ((i*b2*(1-x))/x);
	}
	else if(i <= 1) {
		return ((j*(j-1)/2)/(1-x)) + ((j*b1*x)/(1-x)) + ((i*b2*(1-x))/x);
	}
	else if(j <= 1) {
		return ((i*(i-1)/2)/x) + ((j*b1*x)/(1-x)) + ((i*b2*(1-x))/x);
	}
	else {
		return ((i*(i-1)/2)/x) + ((j*(j-1)/2)/(1-x)) + ((j*b1*x)/(1-x)) + ((i*b2*(1-x))/x);
	}
}

double coal1(int i, int j, double x, double b1, double b2)
{
	if(i <= 1) {
		return 0.0;
	}
	else {
		return (i*(i-1)/2)/(x*h(i,j,x,b1,b2));
	}
}

double coal2(int i, int j, double x, double b1, double b2)
{
	if(j <= 1) {
		return 0.0;
	}
	else {
		return (j*(j-1)/2)/((1-x)*h(i,j,x,b1,b2));
	}
}

double mig1(int i, int j, double x, double b1, double b2)
{
	return (i*b2*(1-x)) / (x * h(i, j, x, b1, b2));
}

double mig2(int i, int j, double x, double b1, double b2)
{
	return (j*b1*x) / ((1-x)*h(i, j, x, b1, b2));
}

gsl_vector* getDiagonal(int n)
{
	int i = 0;
	gsl_vector *diag = gsl_vector_alloc(n+1);

	for(i = 0; i < n+1; i++) {
		gsl_vector_set(diag, i, 1.0);
	}

	return diag;
}

gsl_vector* getSubDiagonal(int n, double x, double b1, double b2)
{
	int i = 0;
	gsl_vector *subDiag = gsl_vector_alloc(n);

	for(i = 0; i < n; i++) {
		gsl_vector_set(subDiag, i, -mig1(i+1, n - (i+1), x, b1, b2));
	}

	return subDiag;
}

gsl_vector* getSuperDiagonal(int n, double x, double b1, double b2)
{
	int i = 0;
	gsl_vector *superDiag = gsl_vector_alloc(n);

	for(i = 0; i < n; i++) {
		gsl_vector_set(superDiag, i, -mig2(i, n - i, x, b1, b2));
	}

	return superDiag;
}

double lnPower(double x, double n)
{
	if(x <= 0.0) {
		printf("ERROR: x <= 0.0 in lnPower()\n");//// MAKE THIS MESSAGE A LITTLE BETTER
		exit(-1);
	}

	return ((double)n)*log((double)x);
}

gsl_vector* buildConstHeightVector(int n, double x, double b1, double b2, gsl_vector* L)
{
	int k;
	gsl_vector *result = gsl_vector_alloc(n+1);

	// COMPUTES result[0] = (1 / h(0, n, x, b1, b2))  + coal2(0, n, x, b1, b2) * L[0]
	gsl_vector_set(result, 0, (1 / h(0, n, x, b1, b2))  + coal2(0, n, x, b1, b2) * gsl_vector_get(L,0));
	// COMPUTES result[n] = (1 / h(n, 0, x, b1, b2))  + coal1(n, 0, x, b1, b2) * L[n-1]
	gsl_vector_set(result, n, (1 / h(n, 0, x, b1, b2))  + coal1(n, 0, x, b1, b2) * gsl_vector_get(L,n-1));

	for(k = 1; k <= n -1; k++) {
		// COMPUTES result[k] = (1 / h(k, n - k, x, b1, b2)) + coal1(k, n - k, x, b1, b2) * L[k-1] + coal2(k, n - k, x, b1, b2) * L[k];
		gsl_vector_set(result, k, (1 / h(k, n - k, x, b1, b2)) + coal1(k, n - k, x, b1, b2) * gsl_vector_get(L,k-1) + coal2(k, n - k, x, b1, b2) * gsl_vector_get(L,k));
	}

	return result;
}

double expectedTreeHeightBalancing(int n, double x, double theta1, double theta2, double R)
{
	int k;
	int val = 0;
	double treeHeight = 0;
	double b1 = theta1 + R*(1 - x);
	double b2 = theta2 + R*x;
	gsl_vector *diag;
	gsl_vector *subDiag;
	gsl_vector *superDiag;
	gsl_vector *u;
	gsl_vector *v;

	v = gsl_vector_alloc(2);
	gsl_vector_set(v,0,0);
	gsl_vector_set(v,1,0);
    
	for(k = 2; k <= n; k++) {
		u = buildConstHeightVector(k, x, b1, b2, v);
		gsl_vector_free(v);
		v = gsl_vector_alloc(k+1);

		diag = getDiagonal(k);
		subDiag = getSubDiagonal(k, x, b1, b2);
		superDiag = getSuperDiagonal(k, x, b1, b2);

		val = gsl_linalg_solve_tridiag(diag, superDiag, subDiag, u, v);

		gsl_vector_free(u);
		gsl_vector_free(diag);
		gsl_vector_free(subDiag);
		gsl_vector_free(superDiag);
	}

	for(k = 0; k <= n; k++) {
		treeHeight = treeHeight + gsl_vector_get(v,k) * exp(log_binom[n][k] + lnPower(x, k) + lnPower(1.0-x, n - k));
	}

	gsl_vector_free(v);

	return treeHeight;
}

gsl_vector* buildConstLengthVector(int n, double x, double b1, double b2, gsl_vector* L)
{
	int k;
	gsl_vector *result = gsl_vector_alloc(n+1);

	// COMPUTES result[0] = (n / h(0, n, x, b1, b2))  + coal2(0, n, x, b1, b2) * L[0]
	gsl_vector_set(result, 0, (n / h(0, n, x, b1, b2))  + coal2(0, n, x, b1, b2) * gsl_vector_get(L,0));
	// COMPUTES result[n] = (n / h(n, 0, x, b1, b2))  + coal1(n, 0, x, b1, b2) * L[n-1]
	gsl_vector_set(result, n, (n / h(n, 0, x, b1, b2))  + coal1(n, 0, x, b1, b2) * gsl_vector_get(L,n-1));

	for(k = 1; k <= n -1; k++) {
		// COMPUTES result[k] = (n / h(k, n - k, x, b1, b2)) + coal1(k, n - k, x, b1, b2) * L[k-1] + coal2(k, n - k, x, b1, b2) * L[k];
		gsl_vector_set(result, k, (n / h(k, n - k, x, b1, b2)) + coal1(k, n - k, x, b1, b2) * gsl_vector_get(L,k-1) + coal2(k, n - k, x, b1, b2) * gsl_vector_get(L,k));
	}

	return result;
}

double expectedTreeLengthBalancing(int n, double x, double theta1, double theta2, double R)
{
	int k;
	int val = 0;
	double treeLen = 0;
	double b1 = theta1 + R*(1 - x);
	double b2 = theta2 + R*x;
	gsl_vector *diag;
	gsl_vector *subDiag;
	gsl_vector *superDiag;
	gsl_vector *u;
	gsl_vector *v;

	v = gsl_vector_alloc(2);
	gsl_vector_set(v,0,0);
	gsl_vector_set(v,1,0);

	for(k = 2; k <= n; k++) {
		u = buildConstLengthVector(k, x, b1, b2, v);
		gsl_vector_free(v);
		v = gsl_vector_alloc(k+1);

		diag = getDiagonal(k);
		subDiag = getSubDiagonal(k, x, b1, b2);
		superDiag = getSuperDiagonal(k, x, b1, b2);

		val = gsl_linalg_solve_tridiag(diag, superDiag, subDiag, u, v);

		gsl_vector_free(u);
		gsl_vector_free(diag);
		gsl_vector_free(subDiag);
		gsl_vector_free(superDiag);
	}

	for(k = 0; k <= n; k++) {
		treeLen = treeLen + gsl_vector_get(v,k) * exp(log_binom[n][k] + lnPower(x, k) + lnPower(1.0-x, n - k));
	}

	gsl_vector_free(v);

	return treeLen;
}

void InitPolySubMatrixTrans(double theta1, double theta2) 
{	
    int *nStack = NULL;
	int index, i, N, Nmax;
    int n = 0;
	int j = 0;
	int k = 0;
	int l = 0;
	int count = 0;
	double x = 0.0;	
	double r = 0.0;
	double q1 = 0.0;
	double q2 = 0.0;
	double ***treeHeight = NULL;
	double ***treeLen = NULL;
    int n_low = n_min[0];
    int n_high = n_max[0];
    double totTreeLen = 0.0;

    for(i = 0; i < num_species; i++) {
        if(n_min[i] < n_low) {
            n_low = n_min[i];
        }

        if(n_high < n_max[i]) {
            n_high = n_max[i];
        }
    }
	
	printf("Initializing matrix of transspecies probabilities\n");
	// Identify the number of equilibrium values
	for(x = e_min; x <= e_max; x = x + e_inc) {
		num_equilibria++;
	}

	// Initialize equilibria grid
	equilibria_grid = (double*)malloc(num_equilibria*sizeof(double));
	for(x = e_min; x <= e_max; x = x + e_inc) {
		equilibria_grid[count] = x;	
		count++;
	}

	// Identify the number of recombination values
	for(r = 0.0; r <= 1.0; r = r + 0.001) {
		num_rates++;
	}
	for(r = 1.1; r <= 10.0; r = r + 0.1) {
		num_rates++;
	}
	for(r = 20.0; r <= 100.0; r = r + 10.0) {
		num_rates++;
	}

	// Initialize recombination rate grid
	rate_grid = (double*)malloc(num_rates*sizeof(double));
	count = 0;	
	for(r = 0.0; r <= 1.0; r = r + 0.001) {
		rate_grid[count] = r;
		count++;
	}
	for(r = 1.1; r <= 10.0; r = r + 0.1) {
		rate_grid[count] = r;
		count++;
	}
	for(r = 20.0; r <= 100.0; r = r + 10.0) {
		rate_grid[count] = r;
		count++;
	}

	treeHeight = (double***)malloc((n_high+1)*sizeof(double**));
	treeLen = (double***)malloc((n_high+1)*sizeof(double**));

	for(n = n_low; n <= n_high; n++) {	
		treeHeight[n] = (double**)malloc(num_equilibria*sizeof(double*));
		treeLen[n] = (double**)malloc(num_equilibria*sizeof(double*));

		for(k = 0; k < num_equilibria; k++) {
			treeHeight[n][k] = (double*)malloc(num_rates*sizeof(double));
			treeLen[n][k] = (double*)malloc(num_rates*sizeof(double));

			for(l = 0; l < num_rates; l++) {
                treeHeight[n][k][l] = expectedTreeHeightBalancing(n, equilibria_grid[k], theta1, theta2, rate_grid[l]);
				treeLen[n][k][l] = expectedTreeLengthBalancing(n, equilibria_grid[k], theta1, theta2, rate_grid[l]);
			}
		}
	}

	// Cycle through to create all necessary values of sample sizes   
    nStack = (int*)malloc(num_species * sizeof(int));
   
    Nmax = 0;
    for(i = 0; i < num_species; i++) { // Get row in the transconfiguration array
        Nmax = Nmax +  ((int)pow(n_max[i] - n_min[i] + 1, num_species-(i+1))) * (n_max[i] - n_min[i]);
    }

    // Initialize probability matrix
	prob_poly_sub_trans = (double****)malloc((Nmax+1)*sizeof(double***));

    // Initialize a stack with minimum sample sizes
    for(index = 0; index < num_species; index++) {
        nStack[index] = n_min[index];
    }

    index = num_species-1;
    while(1) {
        N = 0;
        for(i = 0; i < num_species; i++) { // Get row in the transconfiguration array
            N = N + ((int)pow(n_max[i] - n_min[i] + 1, num_species-(i+1))) * (nStack[i] - n_min[i]);
        }   
      
        prob_poly_sub_trans[N] = (double***)malloc((3*num_species - 1)*sizeof(double**));

        for(j = 0; j < 3*num_species - 3; j++) {
            prob_poly_sub_trans[N][j] = (double**)malloc(num_equilibria*sizeof(double*));

            for(k = 0; k < num_equilibria; k++) {
                prob_poly_sub_trans[N][j][k] = (double*)malloc(num_rates*sizeof(double));

                for(l = 0; l < num_rates; l++) {
                    totTreeLen = 0.0;
                    for(i = 0; i < 2*num_species - 1; i++) {
                        totTreeLen = totTreeLen + branches[i];
                    }
                    for(i = 0; i < num_species; i++) {
                        totTreeLen = totTreeLen + treeLen[nStack[i]][k][l] - treeHeight[nStack[i]][k][l];
                    }
                    
                    // COMPUT THE PROBABILITIES
                    if(j < num_species) {
                        prob_poly_sub_trans[N][j][k][l] = treeLen[nStack[j]][k][l] / totTreeLen; // Polymorphism
                    }
                    else if(num_species <= j && j < 2*num_species) {
                        prob_poly_sub_trans[N][j][k][l] = (branches[j-num_species] - treeHeight[nStack[j-num_species]][k][l]) / totTreeLen; //  Substitution on external branch
                    }
                    else {
                        prob_poly_sub_trans[N][j][k][l] = branches[j-num_species] / totTreeLen; //  Substitution on external branch
                    }
                
        //            printf("N=%d", N);
        //            for(i = 0; i < num_species; i++) {
        //                printf(" %d", nStack[i]);
        //            }
        //            printf(" %lf %lf  CONFIG(%d)  %lf\n", equilibria_grid[k], rate_grid[l], j, prob_poly_sub_trans[N][j][k][l]);

                }
            }
        }
        
        while(index >= 0 && nStack[index] == n_max[index]) {
            nStack[index] = n_min[index];
            index--;
        }
        
        if(index >= 0) {
            nStack[index]++;
            index++;
        }
        
        while(index < num_species-1) {
            nStack[index] = n_min[index];
            index++;
        }
        
        if(N == Nmax) {
            break;
        }
    }
    
	for(n = n_low; n <= n_high; n++) {	
		for(k = 0; k < num_equilibria; k++) {
            free(treeHeight[n][k]);
            free(treeLen[n][k]);
		}
    
        free(treeHeight[n]);
        free(treeLen[n]);
	}
    free(treeHeight);
    free(treeLen);
}

void LoadMarginalSpectra(char *path)
{
	FILE *infile;
	int n, m, i, j, k, l;	
	double total;
	char filename[1000];
	
    simulated_spectrum = (double*****)malloc(num_species * sizeof(double****));

    for(i = 0; i < num_species; i++) {
        simulated_spectrum[i] = (double****)malloc(n_max[i] * sizeof(double***));    

        for(n = n_min[i]; n <= n_max[i]; n++) {
            simulated_spectrum[i][n] = (double***)malloc(num_equilibria * sizeof(double**));
            
            for(k = 0; k < num_equilibria; k++) {
                simulated_spectrum[i][n][k] = (double**)malloc(num_rates * sizeof(double*));
            
                sprintf(filename, "%s/n%d/spectrum_n%d_x%d.txt", path, n, n, 5*k + 5);
            //	printf("Opening %s\n", filename);
				
	            infile = fopen(filename, "r");

                for(l = 0; l < num_rates; l++) {
    				simulated_spectrum[i][n][k][l] = (double*)malloc(n * sizeof(double));     

                    for(m = 0; m < n; m++) {
				        fscanf(infile, "%lf", &simulated_spectrum[i][n][k][l][m]); // Read a line from the file
			        }

			        total = 0.0;
			        for(m = 1; m < n; m++) {
				        total = total + simulated_spectrum[i][n][k][l][m];	
			        }
			
			        for(m = 1; m < n; m++) {
				        simulated_spectrum[i][n][k][l][m] = simulated_spectrum[i][n][k][l][m] / total;	
			        }   
                }

                fclose(infile);
            }        
        }
    }
}

double calc_ln_like_balsel_trans(int marker, int windowRadius, int x, double theta1, double theta2)
{
	double ln_like = 0.0;
	double recombinationRate = 0.0;
	//int N = 0;
	int i = 0, j = 0;
	int r = 0; // VARIABLE FOR RECOMBINATION RATE INDEX
	double interval_factor = 0.0;
    double check;
    int flag_rec = 1;
    
	for(i = marker - windowRadius; i <= marker + windowRadius; i++) {
        flag_rec = 1;
		recombinationRate = fabs(data[i].rate - data[marker].rate);
				
		// If the recombination rate is higher than the number of rates in the file, just use the spectrum for the largest rate.
		// This should be fine because the largest rate should be large enough such that you are close to neutrality. 
		if(rate_grid[num_rates-1] < recombinationRate) { 
			r = num_rates - 2;
            flag_rec = 0;
		}
				
        if(flag_rec) {
		    // Find recombination rate index that fits recombination rate
		    for(r = 0; r < num_rates-1; r++) {
			    // If found index
			    if(rate_grid[r] <= recombinationRate && recombinationRate <= rate_grid[r + 1]) {
				    break;				
			    }
		    }
        }
       
		interval_factor = (recombinationRate - rate_grid[r]) / (rate_grid[r + 1] - rate_grid[r]);

//        N = 0;
//        for(j = 0; j < num_species; j++) { // Get row in the transconfiguration array
//            N = N + ((int)pow(n_max[j] - n_min[j] + 1, num_species-(j+1))) * (data[i].n[j] - n_min[j]);
//        }

        // Get the likelihood of the configuration under balancing selection
        ln_like = ln_like + log(prob_poly_sub_trans[data[i].N][data[i].config][x][r] + interval_factor * (prob_poly_sub_trans[data[i].N][data[i].config][x][r + 1] - prob_poly_sub_trans[data[i].N][data[i].config][x][r]));
	}

	return ln_like;
}

double calc_ln_like_background_trans(int marker, int windowRadius, double **transConfigFreqs)
{
	double ln_like = 0.0;
	int i = 0, j = 0;
//    int N = 0;

	for(i = marker - windowRadius; i <= marker + windowRadius; i++) {
//        N = 0;
//        for(j = 0; j < num_species; j++) { // Get row in the transconfiguration array
//            N = N + ((int)pow(n_max[j] - n_min[j] + 1, num_species-(j+1))) * (data[i].n[j] - n_min[j]);
//        }   

        ln_like = ln_like + log(transConfigFreqs[data[i].N][data[i].config]); // Get the likelihood of the configuration from whole genome diversity
	}

	return ln_like;
}

void scan_with_T1_trans(char *outfn, int windowRadius, double theta1, double theta2, char *transconfigfn)
{
	int startPos = windowRadius;
	int endPos = num_sites - windowRadius - 1;
	int i = 0;
	FILE* outfile;
	double ln_like_balsel, ln_like_background, ln_lr;
	double ln_like_temp;	
	int x_max, x;	
	double **transConfigFreqs;
	
	transConfigFreqs = read_trans_configs(transconfigfn);
	
	outfile=fopen(outfn, "w");

	printf("Performing scan using T1 for transspecific polymorphisms and writing results to %s...\n", outfn);

	for(i = startPos; i <= endPos; i++) {
		x_max = equilibria_grid[0];

		ln_like_balsel = calc_ln_like_balsel_trans(i, windowRadius, x_max, theta1, theta2);
		for(x = 1; x < num_equilibria; x++) {
			ln_like_temp = calc_ln_like_balsel_trans(i, windowRadius, x, theta1, theta2);		
			
			if(ln_like_balsel < ln_like_temp) {
				ln_like_balsel = ln_like_temp;
				x_max = x;
			}		
		}

		ln_like_background = calc_ln_like_background_trans(i, windowRadius, transConfigFreqs);

		ln_lr = 2*(ln_like_balsel - ln_like_background);

		fprintf(outfile, "%lf\t%lf\n", data[i].loc, ln_lr);
	}

	fclose(outfile);

	freeTransConfig(transConfigFreqs);
}





double calc_ln_like_balsel_spectrum_trans(int marker, int windowRadius, int x, double theta1, double theta2)
{
	double ln_like = 0.0;
	double recombinationRate = 0.0;
	//int N = 0;
	int i = 0, j = 0;
	int r = 0; // VARIABLE FOR RECOMBINATION RATE INDEX
	double interval_factor = 0.0;
    double probability1 = 0.0;
	double probability2 = 0.0;
    double check;
    int flag_rec = 1;
    
	for(i = marker - windowRadius; i <= marker + windowRadius; i++) {
        flag_rec = 1;
		recombinationRate = fabs(data[i].rate - data[marker].rate);
				
		// If the recombination rate is higher than the number of rates in the file, just use the spectrum for the largest rate.
		// This should be fine because the largest rate should be large enough such that you are close to neutrality. 
		if(rate_grid[num_rates-1] < recombinationRate) { 
			r = num_rates - 2;
            flag_rec = 0;
		}
				
        if(flag_rec) {
		    // Find recombination rate index that fits recombination rate
		    for(r = 0; r < num_rates-1; r++) {
			    // If found index
			    if(rate_grid[r] <= recombinationRate && recombinationRate <= rate_grid[r + 1]) {
				    break;				
			    }
		    }
        }
       
		interval_factor = (recombinationRate - rate_grid[r]) / (rate_grid[r + 1] - rate_grid[r]);

 //       N = 0;
 //       for(j = 0; j < num_species; j++) { // Get row in the transconfiguration array
 //           N = N + ((int)pow(n_max[j] - n_min[j] + 1, num_species-(j+1))) * (data[i].n[j] - n_min[j]);
 //       }

        // Get the likelihood of the configuration under balancing selection
        if(data[i].config < num_species) { // Get the likelihood of the configuration under balancing selection (polymorphism)
            probability1 = prob_poly_sub_trans[data[i].N][data[i].config][x][r] * simulated_spectrum[data[i].config][data[i].n[data[i].config]][x][r][data[i].n[data[i].config]-data[i].x[data[i].config]];           
            probability2 = prob_poly_sub_trans[data[i].N][data[i].config][x][r + 1] * simulated_spectrum[data[i].config][data[i].n[data[i].config]][x][r + 1][data[i].n[data[i].config]-data[i].x[data[i].config]];    
            ln_like = ln_like + log( probability1 + interval_factor * (probability2 - probability1) );       
        }
        else { // Get the likelihood of the configuration under balancing selection (substitution)
            probability1 = prob_poly_sub_trans[data[i].N][data[i].config][x][r];           
            probability2 = prob_poly_sub_trans[data[i].N][data[i].config][x][r + 1];    
            ln_like = ln_like + log( probability1 + interval_factor * (probability2 - probability1) );
        }
	}

	return ln_like;
}

double calc_ln_like_background_spectrum_trans(int marker, int windowRadius, double **transConfigFreqs)
{
	double ln_like = 0.0;
	int i = 0, j = 0;
//    int N = 0;

	for(i = marker - windowRadius; i <= marker + windowRadius; i++) {
//        N = 0;
//        for(j = 0; j < num_species; j++) { // Get row in the transconfiguration array
//            N = N + ((int)pow(n_max[j] - n_min[j] + 1, num_species-(j+1))) * (data[i].n[j] - n_min[j]);
//        }   

        if(data[i].config < num_species) { // Get the likelihood of the configuration from whole genome diversity (polymorphism)
            ln_like = ln_like + log(transConfigFreqs[data[i].N][data[i].config] * empiricalSpectra[data[i].config][data[i].n[data[i].config]][data[i].x[data[i].config]]); 
        }
        else { // Get the likelihood of the configuration from whole genome diversity (substitution)
            ln_like = ln_like + log(transConfigFreqs[data[i].N][data[i].config]);
        }   
	}

	return ln_like;
}

void scan_with_T2_trans(char *outfn, int windowRadius, double theta1, double theta2, char *transconfigfn)
{
	int startPos = windowRadius;
	int endPos = num_sites - windowRadius - 1;
	int i = 0;
	FILE* outfile;
	double ln_like_balsel, ln_like_background, ln_lr;
	double ln_like_temp;	
	int x_max, x;	
	double **transConfigFreqs;
	
	transConfigFreqs = read_trans_configs(transconfigfn);
	
	outfile=fopen(outfn, "w");

	printf("Performing scan using T2 for transspecific polymorphisms and writing results to %s...\n", outfn);

	for(i = startPos; i <= endPos; i++) {
		x_max = equilibria_grid[0];

		ln_like_balsel = calc_ln_like_balsel_spectrum_trans(i, windowRadius, x_max, theta1, theta2);
		for(x = 1; x < num_equilibria; x++) {
			ln_like_temp = calc_ln_like_balsel_spectrum_trans(i, windowRadius, x, theta1, theta2);		
			
			if(ln_like_balsel < ln_like_temp) {
				ln_like_balsel = ln_like_temp;
				x_max = x;
			}		
		}

		ln_like_background = calc_ln_like_background_spectrum_trans(i, windowRadius, transConfigFreqs);

		ln_lr = 2*(ln_like_balsel - ln_like_background);

		fprintf(outfile, "%lf\t%lf\n", data[i].loc, ln_lr);
	}

	fclose(outfile);

	freeTransConfig(transConfigFreqs);
}




















void get_min_max_sample_sizes(void)
{
	int i, j;

    n_min = (int*)malloc(num_species * sizeof(int));
    n_max = (int*)malloc(num_species * sizeof(int));

	if(num_sites==0) {
        for(i = 0; i < num_species; i++) {
		    n_min[i] = 0;
		    n_max[i] = 0;
        }
		return;
	}

    for(i = 0; i < num_species; i++) {
	    n_min[i] = data[0].n[i];
	    n_max[i] = data[0].n[i];
    }

	for(i = 0; i < num_sites; i++) {
        for(j = 0; j < num_species; j++) {
		    if(data[i].n[j] > n_max[j]) {
			    n_max[j] = data[i].n[j];
		    }

		    if(data[i].n[j] < n_min[j]) {
			    n_min[j] = data[i].n[j];
		    }
        }
	}
}

void get_min_max_rho(void)
{
	int i, j;
    
	if(num_sites==0) {
        recratemin = 0.0;
        recratemax = 0.0;
		return;
	}

    recratemin = data[0].rate;
    recratemax = data[0].rate;
	
	for(i = 0; i < num_sites; i++) {
        if(data[i].rate > recratemax) {
		    recratemax = data[i].rate;
	    }

	    if(data[i].rate < recratemin) {
		    recratemin = data[i].rate;
	    }
	}
}

void read_tree_file(char *treefn)
{
	FILE *infile;
    int i = 0, j = 0;
    	   
    branches = (double*)malloc((2*num_species-3) * sizeof(double));
	topology = (int**)malloc((2*num_species - 3) * sizeof(int*));

    for(i = 0; i < 2*num_species - 3; i++) {
        topology[i] = (int*)malloc(num_species * sizeof(int));

        for(j = 0; j < 2*num_species - 3; j++) {
            topology[i][j] = 0;
        }    
        branches[i] = 0.0;
    }
 
	infile=fopen(treefn, "r");

    for(i = 0; i < 2*num_species - 3; i++) { 
        // Read in species presence
        for(j = 0; j < num_species; j++) {
			fscanf(infile, "%d", &topology[i][j]);
        }
    
        // Read in branch length
        fscanf(infile, "%lf", &branches[i]);
	}

	fclose(infile);
}

void assignConfig(void)
{
	int i = 0, j = 0, k = 0;
    int N = 0;
    int flag_poly = 1;
    int flag_sub = 1;
    int x_sub; // If first species is in clade, then x_sub = 0 if x=0 and x_sub=1 if x=n

   	for(i = 0; i < num_sites; i++) {
        N = 0;
        for(j = 0; j < num_species; j++) { // Get row in the transconfiguration array
            N = N + ((int)pow(n_max[j] - n_min[j] +1, num_species-(j+1))) * (data[i].n[j] - n_min[j]);
        }       
        
        data[i].N = N; // Find the sample size index

        flag_poly = 1;
        for(j = 0; j < num_species; j++) { // Check to see if polymorphism and, if so, then what type
            if(0 < data[i].x[j] && data[i].x[j] < data[i].n[j]) { // Polymorphism in species j+1
                data[i].config = j; // Find the sample configuration (here being one of K polymorphic configurations)
                flag_poly = 0;
                break;
            } 
        }

        // Check to see type of substitution
        if(flag_poly) { // If did not find a polymorphic site            
            if(data[i].x[0] == 0) {
                x_sub = 0;
            }
            else if(data[i].x[0] == data[i].n[0]) {
                x_sub = 1;
            }
            else {
                printf("ERROR: Expecting substitution. Strange sample configuration (x[0],n[0]) = (%d,%d)\n", data[i].x[0], data[i].n[0]);
                exit(-1);
            }
    
            for(j = 0; j < 2*num_species - 3; j++) {
                flag_sub = 0;
                for(k = 0; k < num_species; k++) {
                    if(topology[j][k] == topology[j][0]) {
                        if( (x_sub == 0 && data[i].x[k] > 0) || (x_sub == 1 && data[i].x[k] < data[i].n[k]) ) {  // Does not match split
                            flag_sub = 1;                          
                            break;                        
                        }
                    }
                    else if(topology[j][k] != topology[j][0]) {
                        if( (x_sub == 0 && data[i].x[k] < data[i].n[k]) || (x_sub == 1 && data[i].x[k] > 0) ) {  // Does not match split
                            flag_sub = 1;                        
                            break;                        
                        }
                    }
                    else {
                        printf("ERROR: There should not be any other cases\n");
                        exit(-1);
                    }
                }
            
                if(flag_sub == 0) { // Matches split j
                    data[i].config = num_species + j; // Find the sample configuration (here being one of 2K-3 substitution configurations)
                    break;
                }
            }
        }

        if(flag_poly == 1 && flag_sub == 1) {
            printf("ERROR: Data configuration does not match tree\n");
            printf("x_sub = %d\n", x_sub);
            printf("%lf", data[i].loc);            
            for(j = 0; j < num_species; j++) {
                printf("\t%d\t%d", data[i].x[j], data[i].n[j]);
            }
            printf("\n");
            exit(-1);
        }
	}
}

void readsnps(char *infn)
{
	FILE *infile=fopen(infn, "r");
	int pos=0, col=0, i, j;
	char c, str[1000];

	printf("Reading joint genetic variation file %s\n", infn);

	if(num_sites > 0) {
		free(data);
	}

	num_sites=0;
	c=fgetc(infile);
    
	while(EOF!=(c=fgetc(infile))) {
		if(c=='\n') {
			num_sites++;
		}
	}

	fclose(infile);
	
	data = malloc(num_sites*sizeof(struct datatype));
    for(i=0; i<num_sites; i++) {
        data[i].x = malloc(num_species*sizeof(int));
        data[i].n = malloc(num_species*sizeof(int));    
    }

    infile=fopen(infn, "r");

	for(i=0; i<num_sites; i++) {
        fscanf(infile, "%lf", &data[i].loc);
        fscanf(infile, "%lf", &data[i].rate);
	      	
        for(j=0; j< num_species; j++) {
			fscanf(infile, "%i", &data[i].x[j]);
            fscanf(infile, "%i", &data[i].n[j]);
		}
	}

	fclose(infile);

	get_min_max_sample_sizes();
    get_min_max_rho();

	printf("\tDone reading SNPs (num_sites = %d)\n", num_sites);
    for(i = 0; i < num_species; i++) {
        printf("\tDone reading SNPs (min_sample_size[%d], max_sample_size[%d]) = (%d, %d)\n", i+1, i+1, n_min[i], n_max[i]);
    }
    printf("\tDone reading SNPs (recratemin, recratemax) = (%lf, %lf)\n", recratemin, recratemax);
}

double** read_trans_configs(char *transconfigfn)
{
	FILE *infile;
    int i = 0, j = 0;
    double **freq = NULL;
	int sample_size;
    int dimFreq = 1;	

    dimFreq = 1;
    for(i = 0; i < num_species; i++) {
        dimFreq = dimFreq * (n_max[i] - n_min[i] + 1);
    }   

	freq = (double**)malloc(dimFreq * sizeof(double*)); // One row per set of sample sizes

    for(i = 0; i < dimFreq; i++) {
        freq[i] = (double*)malloc((3*num_species-3) * sizeof(double)); // K polymorphic and 2K-3 substitution
        for(j = 0; j < 3*num_species - 3; j++) {
            freq[i][j] = 0.0;
        }    
    }
 
	infile=fopen(transconfigfn, "r");

    for(i = 0; i < dimFreq; i++) { 
        // Read in sample sizes
        for(j = 0; j < num_species; j++) {
			fscanf(infile, "%d", &sample_size);
        }
    
        // Read in 3K-3 configuration values
        for(j = 0; j < 3*num_species - 3; j++) {
            fscanf(infile, "%lf", &freq[i][j]);	
		}
	}

	fclose(infile);

	return freq;
}

double **calcTransConfigFreq(void)
{
	int i = 0, j = 0, k = 0;
	double **freq = NULL;
    int dimFreq = 1;
	int **num = NULL;
	int *denom = NULL;
    int N = 0;
    int flag_poly = 1;
    int flag_sub = 1;
    int x_sub; // If first species is in clade, then x_sub = 0 if x=0 and x_sub=1 if x=n

    dimFreq = 1;
    for(i = 0; i < num_species; i++) {
        dimFreq = dimFreq * (n_max[i] - n_min[i] + 1);
    }   

    printf("dimFreq = %d\n", dimFreq);
    
	freq = (double**)malloc(dimFreq * sizeof(double*)); // One row per set of sample sizes
    num = (int**)malloc(dimFreq * sizeof(int*)); // One row per set of sample sizes
    denom = (int*)malloc(dimFreq * sizeof(int)); // One element per set of sample sizes

    for(i = 0; i < dimFreq; i++) {
        freq[i] = (double*)malloc((3*num_species-3) * sizeof(double)); // K polymorphic and 2K-3 substitution
        num[i] = (int*)malloc((3*num_species-3) * sizeof(int)); // K polymorphic and 2K-3 substitution
                
        for(j = 0; j < 3*num_species - 3; j++) {
            freq[i][j] = 0.0;
            num[i][j] = 0;
        }    
        denom[i] = 0;
    }

	for(i = 0; i < num_sites; i++) {
        num[data[i].N][data[i].config]++;
        denom[data[i].N]++;
	}

	for(i = 0; i < dimFreq; i++) {
		if(denom[i] > 0) {
			for(j = 0; j < 3*num_species - 3; j++) {
				freq[i][j] = ((double) num[i][j]) / ((double) denom[i]);
			}
		}
		else {
			for(j = 0; j < 3*num_species - 3; j++) {
				freq[i][j] = 0.0;
			}
		}
	}

    for(i = 0; i < dimFreq; i++) {
        free(num[i]);
    }
    free(num);
    free(denom);

	return freq;
}

/*double **calcTransConfigFreq(void)
{
	int i = 0, j = 0, k = 0;
	double **freq = NULL;
    int dimFreq = 1;
	int **num = NULL;
	int *denom = NULL;
    int N = 0;
    int flag_poly = 1;
    int flag_sub = 1;
    int x_sub; // If first species is in clade, then x_sub = 0 if x=0 and x_sub=1 if x=n

    dimFreq = 1;
    for(i = 0; i < num_species; i++) {
        dimFreq = dimFreq * (n_max[i] - n_min[i] + 1);
    }   

    printf("dimFreq = %d\n", dimFreq);
    
	freq = (double**)malloc(dimFreq * sizeof(double*)); // One row per set of sample sizes
    num = (int**)malloc(dimFreq * sizeof(int*)); // One row per set of sample sizes
    denom = (int*)malloc(dimFreq * sizeof(int)); // One element per set of sample sizes

    for(i = 0; i < dimFreq; i++) {
        freq[i] = (double*)malloc((3*num_species-3) * sizeof(double)); // K polymorphic and 2K-3 substitution
        num[i] = (int*)malloc((3*num_species-3) * sizeof(int)); // K polymorphic and 2K-3 substitution
                
        for(j = 0; j < 3*num_species - 3; j++) {
            freq[i][j] = 0.0;
            num[i][j] = 0;
        }    
        denom[i] = 0;
    }
 
	for(i = 0; i < num_sites; i++) {
        N = 0;
        for(j = 0; j < num_species; j++) { // Get row in the transconfiguration array
            N = N + ((int)pow(n_max[j] - n_min[j] +1, num_species-(j+1))) * (data[i].n[j] - n_min[j]);
        }       
        
        data[i].N = N; // Find the sample size index

        flag_poly = 1;
        for(j = 0; j < num_species; j++) { // Check to see if polymorphism and, if so, then what type
            if(0 < data[i].x[j] && data[i].x[j] < data[i].n[j]) { // Polymorphism in species j+1
                num[N][j]++;
                data[i].config = j; // Find the sample configuration (here being one of K polymorphic configurations)
                flag_poly = 0;
                break;
            } 
        }

        // Check to see type of substitution
        if(flag_poly) { // If did not find a polymorphic site            
            if(data[i].x[0] == 0) {
                x_sub = 0;
            }
            else if(data[i].x[0] == data[i].n[0]) {
                x_sub = 1;
            }
            else {
                printf("ERROR: Expecting substitution. Strange sample configuration (x[0],n[0]) = (%d,%d)\n", data[i].x[0], data[i].n[0]);
                exit(-1);
            }
    
            for(j = 0; j < 2*num_species - 3; j++) {
                flag_sub = 0;
                for(k = 0; k < num_species; k++) {
                    if(topology[j][k] == topology[j][0]) {
                        if( (x_sub == 0 && data[i].x[k] > 0) || (x_sub == 1 && data[i].x[k] < data[i].n[k]) ) {  // Does not match split
                            flag_sub = 1;                          
                            break;                        
                        }
                    }
                    else if(topology[j][k] != topology[j][0]) {
                        if( (x_sub == 0 && data[i].x[k] < data[i].n[k]) || (x_sub == 1 && data[i].x[k] > 0) ) {  // Does not match split
                            flag_sub = 1;                        
                            break;                        
                        }
                    }
                    else {
                        printf("ERROR: There should not be any other cases\n");
                        exit(-1);
                    }
                }
            
                if(flag_sub == 0) { // Matches split j
                    num[N][num_species + j]++;
                    data[i].config = num_species + j; // Find the sample configuration (here being one of 2K-3 substitution configurations)
                    break;
                }
            }
        }

        if(flag_poly == 1 && flag_sub == 1) {
            printf("ERROR: Data configuration does not match tree\n");
            printf("x_sub = %d\n", x_sub);
            printf("%lf", data[i].loc);            
            for(j = 0; j < num_species; j++) {
                printf("\t%d\t%d", data[i].x[j], data[i].n[j]);
            }
            printf("\n");
            exit(-1);
        }
    
        denom[N]++;
	}

	for(i = 0; i < dimFreq; i++) {
		if(denom[i] > 0) {
			for(j = 0; j < 3*num_species - 3; j++) {
				freq[i][j] = ((double) num[i][j]) / ((double) denom[i]);
			}
		}
		else {
			for(j = 0; j < 3*num_species - 3; j++) {
				freq[i][j] = 0.0;
			}
		}
	}

    for(i = 0; i < dimFreq; i++) {
        free(num[i]);
    }
    free(num);
    free(denom);

	return freq;
}*/

void printFractionTransConfig(char *outfn, double **transConfigFreq)
{
	int *nStack = NULL;  
    int index, i, N, Nmax;
	FILE *outfile = fopen(outfn, "w");
    
    nStack = (int*)malloc(num_species * sizeof(int));
   
    Nmax = 0;
    for(i = 0; i < num_species; i++) { // Get row in the transconfiguration array
        Nmax = Nmax +  ((int)pow(n_max[i] - n_min[i] + 1, num_species-(i+1))) * (n_max[i] - n_min[i]);
    }

    // Initialize a stack with minimum sample sizes
    for(index = 0; index < num_species; index++) {
        nStack[index] = n_min[index];
    }

    index = num_species-1;
    while(1) {
        N = 0;
        for(i = 0; i < num_species; i++) { // Get row in the transconfiguration array
            N = N + ((int)pow(n_max[i] - n_min[i] + 1, num_species-(i+1))) * (nStack[i] - n_min[i]);
        }   
      
        // Print out sample sizes
        fprintf(outfile, "%d", nStack[0]);
        for(i = 1; i < num_species; i++) {
            fprintf(outfile, "\t%d", nStack[i]);
        }

        // Print out transconfig frequencies
        for(i = 0; i < 3*num_species - 3; i++) {
            fprintf(outfile, "\t%lf", transConfigFreq[N][i]);
        }
        fprintf(outfile, "\n");

        while(index >= 0 && nStack[index] == n_max[index]) {
            nStack[index] = n_min[index];
            index--;
        }
        
        if(index >= 0) {
            nStack[index]++;
            index++;
        }

        while(index < num_species-1) {
            nStack[index] = n_min[index];
            index++;
        }

        if(N == Nmax) {
            break;
        }
    }

	fclose(outfile);
}

void freeTransConfig(double **transConfigFreq)
{
	int i, Nmax;
	
    Nmax = 0;
    for(i = 0; i < num_species; i++) { // Get row in the transconfiguration array
        Nmax = Nmax +  ((int)pow(n_max[i] - n_min[i] + 1, num_species-(i+1))) * (n_max[i] - n_min[i]);
    }

    for(i = 0; i <= Nmax; i++) {
        free(transConfigFreq[i]);
    }
    free(transConfigFreq);
}

void read_marginal_empirical_spectra(char **spectfn)
{
	FILE* infile;
	int i, j, k, minSample, maxSample;
	double total = 0.0;
    

    empiricalSpectra = (double ***)malloc(num_species * sizeof(double **));
    for(i = 0; i < num_species; i++) {
        infile = fopen(spectfn[i], "r");
        
        fscanf(infile, "%d", &minSample);
	    fscanf(infile, "%d", &maxSample);
        empiricalSpectra[i] = (double **)malloc((maxSample+1) * sizeof(double *));
     
        for(j = minSample; j <= maxSample; j++) {
            empiricalSpectra[i][j] = (double *)malloc((j+1) * sizeof(double));
            
            total = 0.0;
            for(k = 1; k <= j - 1; k++) {
                fscanf(infile, "%lf", &empiricalSpectra[i][j][j - k]);
			    total = total + empiricalSpectra[i][j][j - k];	
            }

            for(k = 1; k <= j - 1; k++) {
                empiricalSpectra[i][j][k] = empiricalSpectra[i][j][k] / total;
            }
        }
        fclose(infile);     
    }
}


void obtain_marginal_empirical_spectra(void)
{
	int i, j, k;
	int ***counts;
	int sum = 0;

	empiricalSpectra = (double ***)malloc(num_species * sizeof(double **));
    counts = (int ***)malloc(num_species * sizeof(int **));
    for(i = 0; i < num_species; i++) {
        empiricalSpectra[i] = (double **)malloc((n_max[i]+1) * sizeof(double *));
        counts[i] = (int **)malloc((n_max[i]+1) * sizeof(int *)); 

        for(j = n_min[i]; j <= n_max[i]; j++) {
            empiricalSpectra[i][j] = (double *)malloc((j+1) * sizeof(double));
            counts[i][j] = (int *)malloc((j+1) * sizeof(int));            
        
            for(k = 1; k <= j - 1; k++) {
                empiricalSpectra[i][j][k] = 0.0;
                counts[i][j][k] = 0;
            }
        }     
    }
   
    for(i = 0; i < num_sites; i++) {
        for(j = 0; j < num_species; j++) {
		    if((data[i].x[j] > 0) && (data[i].x[j] < data[i].n[j])) {
			    counts[j][data[i].n[j]][data[i].n[j] - data[i].x[j]]++;	// CONVERT TO DERIVED FREQUENCY SPECTRUM
		    }
        }
	}
	
    for(i = 0; i < num_species; i++) {
	    for(j = n_min[i]; j <= n_max[i]; j++) {
		    sum = 0;

		    for(k = 1; k <= j - 1; k++) {
			    sum = sum + counts[i][j][k];
		    }

		    for(k = 1; k <= j - 1; k++) {
			    if(sum > 0) {
				    empiricalSpectra[i][j][k] = counts[i][j][k] / ((double) sum);
			    }
		    }
	    }
    }

    for(i = 0; i < num_species; i++) {	
        for(j = n_min[i]; j <= n_max[i]; j++) {
	        free(counts[i][j]);
	    }

        free(counts[i]);
    }

	free(counts);
}

void write_empirical_spectra(char **outfn, double ***spectra)
{
	FILE* outfile;
	int i, j, k;

    for(i = 0; i < num_species; i++) {
        outfile=fopen(outfn[i], "w");

    	fprintf(outfile, "%d %d\n", n_min[i], n_max[i]);
	    for(j = n_min[i]; j <= n_max[i]; j++) {
		    for(k = 1; k <= j - 1; k++) {
			    fprintf(outfile, "%lf ", spectra[i][j][k]);
		    }
		    fprintf(outfile, "\n");
	    }
        
        fclose(outfile);        
    }	
}

void free_empirical_spectra(double ***spectra)
{
	int i, j;
    
    for(i = 0; i < num_species; i++) {
       	for(j = n_min[i]; j <= n_max[i]; j++) {
		    free(spectra[i][j]);
	    }
    
        free(spectra[i]);
    }

	free(spectra);
}

void usage() 
{
	printf("\n\n*****************************************\n");
	printf("***************** USAGE *****************\n");
	printf("*****************************************\n\n");
	
	printf("********** Prepare scan **********\n");
    printf("Get proportion of transspecies configuration sites: ./MULLET -config NumSpecies TreeFile CombinedSNPFile ConfigFile\n");	
    printf("Get marginal empirical spectra: ./MULLET -spect NumSpecies CombinedSNPFile SpectFile1 SpectFile2 ... SpectFileK\n");

    printf("\n********** Shared balancing selection **********\n");
	//printf("Generate spectrum: ./MULLET -SimSpect n x theta1 theta2 num_replicates OutFile\n");
	printf("Scan with T1trans: ./MULLET -T1trans WINDOWSIZE NumSpecies TreeFile ConfigFile SNPFile OutFile\n");
    printf("Scan with T2trans: ./MULLET -T2trans WINDOWSIZE NumSpecies TreeFile ConfigFile SpectFile1 ... SpectFileK SNPFile PATH OutFile\n");	
   
	exit(-1);
}
