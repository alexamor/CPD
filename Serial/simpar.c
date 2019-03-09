//SIMPAR.C

// Libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


// Constant variables and functions
#define RND0_1 ((double) random() / ((long long)1<<31))
#define G 6.67408e-11
#define EPSLON 0.01
#define NR_OF_SIDES 9

// struct that stores the information of each particle
typedef struct particle_t{
	double x, y; // Position
	double vx, vy; // Velocity
	double m; // Mass
	double fx, fy; // force
} particle_t;


// struct that stores the information of each mass center in a cell
typedef struct cell{
	double x, y; // Position
	double m; // Mass
} cell;


int main{
	long seed, ncside;
	long long n_part, n_tstep;
	particle_t *par = NULL;
	cell **cell_mat = NULL;
	long long i, j, k, t;
	long long F, d2;

	int column, row, adj_column, adj_row, adj_x, adj_y;

	// check the correct number of arguments
	if( argc != 5){
		printf("Number of arguments invalid. Please check the way you're running the program, as follows: ./simpar seed ncside n_part n_tstep \n" );
		exit(0);
	}

	// saving the input parameters
	seed = atol(argv[1]);
	ncside = atol(argv[2]);
	if (ncside < 3){
		printf("Number of sides must be equal or larger than 3.\n");
		exit(0);
	}
	n_part = atoll(argv[3]);
	n_tstep = atoll(argv[4]);

	// Get random numbers from seed
	srandom(seed);

	// allocation of memory for the particle_t array
	par = (particle_t *) malloc( sizeof(particle_t) * n_part);
	if (par == NULL){
		printf("No memory available for the number of particles required.\n");
		exit(0);
	}

	// creation of particles - initialization of their position, velocity and mass
	for(i = 0; i < n_part; i++)
    {
    	// position
        par[i].x = RND0_1;
        par[i].y = RND0_1;

        // velocity
        par[i].vx = RND0_1 / ncside / 10.0;
        par[i].vy = RND0_1 / ncside / 10.0;

        // mass
        par[i].m = RND0_1 * ncside / (G * 1e6 * n_part);

        // zero the force
        par[i].fx = 0;
        par[i].fy = 0;
    }

    // memory allocation for each cell
    cell_mat = (cell **) malloc( sizeof(cell*) * ncside);
    if (cell_mat == NULL){
    	printf("No memory available for the number of cells required.\n");
    	exit(0);
    }
    for(i = 0; i < ncside; i++){

    	// initialization to zero 
    	cell_mat[i] = (cell *) calloc( ncside, sizeof(cell));

    	if (cell_mat[i] == NULL){
    		printf("No memory available for the number of cells required.\n");
    		exit(0);
    	}
    }

    
    // cycle of time step iterations
    for(t = 0; t < n_tstep; t++){
    	
    	// calculation of the mass center
    	for(i = 0; i < n_part; i++){

    		// get location in grid from the position - truncating the float value 
    		column = (int) par[i].x * ncside;
    		row = (int) par[i].y * ncside;

    		// average calculated progressively without needing to store every x and y value
    		cell_mat[column][row].x = (cell_mat[column][row].x*cell_mat[column][row].m + par[i].m * par[i].x) / (cell_mat[column][row].m + par[i].m);
    		cell_mat[column][row].y = (cell_mat[column][row].y*cell_mat[column][row].m + par[i].m * par[i].y) / (cell_mat[column][row].m + par[i].m);

    		// total mass
    		cell_mat[column][row].m += par[i].m; 
    	}

    	// calculation of the gravitational force
    	for(i = 0; i < n_part; i++){
    		
    		// get location in grid from the position - truncating the float value 
    		column = (int) par[i].x * ncside;
    		row = (int) par[i].y * ncside;

			for(j = -1; j < 2; j++){
				
				adj_column = column + j;

				// TODO ver se conseguimos transformar numa operação matematica atraves do modulo
				// check if adjacent cells are out of the border and rectify
				if(adj_column >= ncside){
					adj_column = 0;
					adj_x = cell_mat[adj_column]
				}
				else if(adj_column < 0){
					adj_column = ncside - 1;
				}

				for(k = -1; k < 2; j++){

					adj_row = row + k;

					// check if adjacent cells are out of the border and rectify
					if(adj_row >= ncside){
						adj_row = 0;
					}
					else if(adj_row < 0 ){
						adj_row = ncside - 1;
					}

					// calculate distance and gravitational force
					d2 = pow( cell_mat[adj_column][adj_row].x - par[i].x, 2) + pow( cell_mat[adj_column][adj_row].y - par[i].y, 2) ;

					// check to see if distance is minimal
					if (d2 < EPSLON)
						F = 0;
					else
						F = G *par[i].m*cell_mat[adj_column][adj_row].m / d2;

					// get cartesian components of the gravitational force
					if ( j == -1 && k == -1 ){
						par[i].fx = 0;
						par[i].fy = 0;
					}

					par[i].fx += F*cos()

				}



			}    		

    	}


    	// zeroing the mass centers 
    	for(i = 0; i < ncside; i++)
    		for(j = 0; j < ncside; j++){
    			cell_mat[i][j].x = 0;
    			cell_mat[i][j].y = 0;
    			cell_mat[i][j].m = 0;
    		}

    
    }



    // freeing the allocated memory
    free(par);
    for (i = 0; i < ncside; i++){
    	free(cell_mat[i]);
    }
    free(cell_mat);

	return 0;
}