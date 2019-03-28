#ifdef _POMP
#  undef _POMP
#endif
#define _POMP 200110

#include "simpar.c.opari.inc"
#line 1 "simpar.c"
//SIMPAR.C

// Libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>



// Constant variables and functions
#define RND0_1 ((double) random() / ((long long)1<<31))
#define G 6.67408e-11
#define EPSLON 0.0005


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


int main(int argc, char *argv[]){
	long seed, ncside;
	long long n_part, n_tstep;
	particle_t *par = NULL;
	cell **cell_mat = NULL;
	int i, j, k, t;
	double F, d2, dx, dy;
	double ax, ay;
	int column, row, adj_column, adj_row, adj_x, adj_y;
	double m = 0, mx = 0, my = 0;

	//matrix of locks
	omp_lock_t **my_lock_matrix = NULL;



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


	// NON PARALLELIZABLE
	// creation of particles - initialization of their position, velocity and mass
	for(i = 0; i < n_part; i++)
    {
    	// position
        par[i].x = RND0_1;
        par[i].y = RND0_1;

        //printf("par x %lf par y %lf \n", par[i].x, par[i].y);

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



POMP_Parallel_fork(&omp_rd_1);
#line 110 "simpar.c"
#pragma omp parallel    
{ POMP_Parallel_begin(&omp_rd_1);
POMP_For_enter(&omp_rd_1);
#line 110 "simpar.c"
#pragma omp          for nowait
    for(i = 0; i < ncside; i++){

    	// initialization to zero
    	cell_mat[i] = (cell *) calloc( ncside, sizeof(cell));

    	if (cell_mat[i] == NULL){
    		printf("No memory available for the number of cells required.\n");
    		exit(0);
    	}
    }
POMP_Barrier_enter(&omp_rd_1);
#pragma omp barrier
POMP_Barrier_exit(&omp_rd_1);
POMP_For_exit(&omp_rd_1);
POMP_Parallel_end(&omp_rd_1); }
POMP_Parallel_join(&omp_rd_1);
#line 121 "simpar.c"



    // memory allocation of the locks matrix
    my_lock_matrix = (omp_lock_t **) malloc( sizeof(omp_lock_t *)* ncside);
    if (my_lock_matrix == NULL){
    	printf("No memory available for the number of locks required.\n");
    	exit(0);
    }
POMP_Parallel_fork(&omp_rd_2);
#line 130 "simpar.c"
#pragma omp parallel    
{ POMP_Parallel_begin(&omp_rd_2);
POMP_For_enter(&omp_rd_2);
#line 130 "simpar.c"
#pragma omp          for nowait
    for(i = 0; i < ncside; i++){

    	my_lock_matrix[i] = (omp_lock_t *) malloc( sizeof(omp_lock_t) * ncside);

	}
POMP_Barrier_enter(&omp_rd_2);
#pragma omp barrier
POMP_Barrier_exit(&omp_rd_2);
POMP_For_exit(&omp_rd_2);
POMP_Parallel_end(&omp_rd_2); }
POMP_Parallel_join(&omp_rd_2);
#line 136 "simpar.c"

POMP_Parallel_fork(&omp_rd_3);
#line 137 "simpar.c"
#pragma omp parallel     private(i,j)
{ POMP_Parallel_begin(&omp_rd_3);
POMP_For_enter(&omp_rd_3);
#line 137 "simpar.c"
#pragma omp          for              nowait
	for (i = 0; i < ncside; i++)
		for (j = 0; j < ncside; j++)
		{
			POMP_Init_lock(&my_lock_matrix[i][j]);
		}
POMP_Barrier_enter(&omp_rd_3);
#pragma omp barrier
POMP_Barrier_exit(&omp_rd_3);
POMP_For_exit(&omp_rd_3);
POMP_Parallel_end(&omp_rd_3); }
POMP_Parallel_join(&omp_rd_3);
#line 143 "simpar.c"

    // cycle of time step iterations
    for(t = 0; t < n_tstep; t++){

    	// zeroing the mass centers and initating locks
POMP_Parallel_fork(&omp_rd_4);
#line 148 "simpar.c"
    #pragma omp parallel     private(i,j)
{ POMP_Parallel_begin(&omp_rd_4);
POMP_For_enter(&omp_rd_4);
#line 148 "simpar.c"
    #pragma omp          for              nowait
    	for(i = 0; i < ncside; i++)
    		for(j = 0; j < ncside; j++){
    			cell_mat[i][j].x = 0;
    			cell_mat[i][j].y = 0;
    			cell_mat[i][j].m = 0;
    		}
POMP_Barrier_enter(&omp_rd_4);
#pragma omp barrier
POMP_Barrier_exit(&omp_rd_4);
POMP_For_exit(&omp_rd_4);
POMP_Parallel_end(&omp_rd_4); }
POMP_Parallel_join(&omp_rd_4);
#line 155 "simpar.c"

    	// calculation of the mass center
POMP_Parallel_fork(&omp_rd_5);
#line 157 "simpar.c"
    #pragma omp parallel     private(i,column,row)
{ POMP_Parallel_begin(&omp_rd_5);
POMP_For_enter(&omp_rd_5);
#line 157 "simpar.c"
    #pragma omp          for                       nowait
    	for(i = 0; i < n_part; i++){

    		// get location in grid from the position - truncating the float value
    		column = (int) (par[i].x * ncside);
    		row = (int) (par[i].y * ncside);



    		// average calculated progressively without needing to store every x and y value
POMP_Atomic_enter(&omp_rd_6);
#line 167 "simpar.c"
				#pragma omp atomic
    		cell_mat[column][row].y += par[i].m * par[i].y;
POMP_Atomic_exit(&omp_rd_6);
#line 169 "simpar.c"
POMP_Atomic_enter(&omp_rd_7);
#line 169 "simpar.c"
				#pragma omp atomic
				cell_mat[column][row].x += par[i].m * par[i].x;
POMP_Atomic_exit(&omp_rd_7);
#line 171 "simpar.c"

    		// total mass
POMP_Atomic_enter(&omp_rd_8);
#line 173 "simpar.c"
				#pragma omp atomic
    		cell_mat[column][row].m += par[i].m;
POMP_Atomic_exit(&omp_rd_8);
#line 175 "simpar.c"

    	}
POMP_Barrier_enter(&omp_rd_5);
#pragma omp barrier
POMP_Barrier_exit(&omp_rd_5);
POMP_For_exit(&omp_rd_5);
POMP_Parallel_end(&omp_rd_5); }
POMP_Parallel_join(&omp_rd_5);
#line 177 "simpar.c"

POMP_Parallel_fork(&omp_rd_9);
#line 178 "simpar.c"
			#pragma omp parallel     private(i, j)
{ POMP_Parallel_begin(&omp_rd_9);
POMP_For_enter(&omp_rd_9);
#line 178 "simpar.c"
   #pragma omp          for               nowait
			for(i = 0; i < ncside; i++)
				for(j = 0; j < ncside; j++){
					cell_mat[i][j].x /= cell_mat[i][j].m;
					cell_mat[i][j].y /= cell_mat[i][j].m;
				}
POMP_Barrier_enter(&omp_rd_9);
#pragma omp barrier
POMP_Barrier_exit(&omp_rd_9);
POMP_For_exit(&omp_rd_9);
POMP_Parallel_end(&omp_rd_9); }
POMP_Parallel_join(&omp_rd_9);
#line 184 "simpar.c"

    	// calculation of the gravitational force
POMP_Parallel_fork(&omp_rd_10);
#line 186 "simpar.c"
#pragma omp parallel     private(i, j, k, column, row, adj_column, adj_row, dx, dy, d2, F)
{ POMP_Parallel_begin(&omp_rd_10);
POMP_For_enter(&omp_rd_10);
#line 186 "simpar.c"
#pragma omp          for                                                                   nowait
    	for(i = 0; i < n_part; i++){

    		// get location in grid from the position - truncating the float value
    		column = (int) (par[i].x * ncside);
    		row = (int) (par[i].y * ncside);

    		//printf("column %d   row %d\n", column, row);
    		//printf("x %.2lf y %.2lf\n", par[i].x, par[i].y );

			for(j = -1; j < 2; j++){

				//printf("j = %d\n", j);

				adj_column = column + j;

				// TODO ver se conseguimos transformar numa operação matematica atraves do modulo
				// check if adjacent cells are out of the border and rectify
				if(adj_column >= ncside){
					adj_column = 0;
				}
				else if(adj_column < 0){
					adj_column = ncside - 1;
				}

				for(k = -1; k < 2; k++){

					//printf("k = %d\n", k);

					adj_row = row + k;

					// check if adjacent cells are out of the border and rectify
					if(adj_row >= ncside){
						adj_row = 0;
					}
					else if(adj_row < 0 ){
						adj_row = ncside - 1;
					}

					// calculate usual distances
					dx = cell_mat[adj_column][adj_row].x - par[i].x;
					dy = cell_mat[adj_column][adj_row].y - par[i].y;

					//printf("mat.y %lf  par.y %lf\n", cell_mat[adj_column][adj_row].y, par[i].y);

					// calculate distances when the cells are out of borders
					if(j == -1 && adj_column == ncside - 1)
						dx =  cell_mat[adj_column][adj_row].x - par[i].x - 1;

					if(j == 1 && adj_column == 0)
						dx = 1 + (cell_mat[adj_column][adj_row].x - par[i].x);

					if(k == 1 && adj_row == 0)
						dy = cell_mat[adj_column][adj_row].y - par[i].y + 1;

					if(k == -1 && adj_row == ncside - 1)
						dy = 1 + (cell_mat[adj_column][adj_row].y - par[i].y);


					//printf("dx %lf     dy %lf\n", dx, dy);


					// calculate distance and gravitational force
					d2 = dx*dx + dy*dy;

					// check to see if distance is minimal
					if (d2 < EPSLON)
						F = 0;
					else
						F = G *par[i].m*cell_mat[adj_column][adj_row].m / d2;

					//Divide by 0 error check
					if(F == 0)
						dx = 1;

					//printf("%lf\n", d2 );

					// get cartesian components of the gravitational force
					if ( j == -1 && k == -1 ){
						par[i].fx = 0;
						par[i].fy = 0;
					}

					//printf("atan  %lf\n", atan(dy/dx));

					// calculate force
					par[i].fx += F*cos(atan(dy/dx));
					par[i].fy += F*sin(atan(dy/dx));

					//printf("fx %lf   fy %lf\n", par[i].fx, par[i].fy);
				}

			}


    	}
POMP_Barrier_enter(&omp_rd_10);
#pragma omp barrier
POMP_Barrier_exit(&omp_rd_10);
POMP_For_exit(&omp_rd_10);
POMP_Parallel_end(&omp_rd_10); }
POMP_Parallel_join(&omp_rd_10);
#line 282 "simpar.c"

    	// calculation of new velocity and position
POMP_Parallel_fork(&omp_rd_11);
#line 284 "simpar.c"
#pragma omp parallel     private(i, ax, ay)
{ POMP_Parallel_begin(&omp_rd_11);
POMP_For_enter(&omp_rd_11);
#line 284 "simpar.c"
#pragma omp          for                    nowait
    	for(i = 0; i < n_part; i++){

    		// get acceleration
    		ax = par[i].fx / par[i].m;
    		ay = par[i].fy / par[i].m;


    		//printf("fx %lf   fy %lf\n", par[i].fx, par[i].fy);
    		//printf("ax %lf ay %lf \n", ax, ay);

    		// position
    		par[i].x += par[i].vx + ax;
    		par[i].y += par[i].vy + ay;

    		if(par[i].x < 0){
    			par[i].x += 1;
    		}
    		if(par[i].x > 1){
    			par[i].x -= 1;
    		}

    		if(par[i].y < 0){
    			par[i].y += 1;
    		}
    		if(par[i].y > 1){
    			par[i].y -= 1;
    		}

    		// velocity
    		par[i].vx += ax;
    		par[i].vy += ay;


    	}
POMP_Barrier_enter(&omp_rd_11);
#pragma omp barrier
POMP_Barrier_exit(&omp_rd_11);
POMP_For_exit(&omp_rd_11);
POMP_Parallel_end(&omp_rd_11); }
POMP_Parallel_join(&omp_rd_11);
#line 319 "simpar.c"

    	printf("Finished iteration: %d\n", t);

    }


    // calculate final mass center
POMP_Parallel_fork(&omp_rd_12);
#line 326 "simpar.c"
#pragma omp parallel     reduction(+:mx,my,m)
{ POMP_Parallel_begin(&omp_rd_12);
POMP_For_enter(&omp_rd_12);
#line 326 "simpar.c"
#pragma omp          for                      nowait
    for(i = 0; i < n_part; i++){

    	    // average calculated progressively without needing to store every x and y value
    		mx += par[i].m * par[i].x;
    		my += par[i].m * par[i].y;

    		//printf("mx %lf my %lf\n", par[i].x, par[i].y );

    		// total mass
    		m += par[i].m;

    }
POMP_Barrier_enter(&omp_rd_12);
#pragma omp barrier
POMP_Barrier_exit(&omp_rd_12);
POMP_For_exit(&omp_rd_12);
POMP_Parallel_end(&omp_rd_12); }
POMP_Parallel_join(&omp_rd_12);
#line 339 "simpar.c"

	//Note: Altrough is faster to average the cells the result ends up poorly rounded, failing the correct answer for +-0.01
    /*
    for(i = 0; i < ncside; i++)
    {
    	for(j = 0; j < ncside; j++)
    	{
    		mx += cell_mat[i][j].x * cell_mat[i][j].m;
    		my += cell_mat[i][j].y * cell_mat[i][j].m;

    		m += cell_mat[i][j].m;
    	}
    }
    */

    mx = mx/m;
    my = my/m;

    // output
    printf("%.2lf %.2lf\n", par[0].x, par[0].y );
    printf("%.2lf %.2lf\n", mx, my );

POMP_Parallel_fork(&omp_rd_13);
#line 361 "simpar.c"
#pragma omp parallel     private(i,j)
{ POMP_Parallel_begin(&omp_rd_13);
POMP_For_enter(&omp_rd_13);
#line 361 "simpar.c"
#pragma omp          for              nowait
	for (i = 0; i < ncside; i++)
		for (j = 0; j < ncside; j++)
		{
			POMP_Destroy_lock(&my_lock_matrix[i][j]);
		}
POMP_Barrier_enter(&omp_rd_13);
#pragma omp barrier
POMP_Barrier_exit(&omp_rd_13);
POMP_For_exit(&omp_rd_13);
POMP_Parallel_end(&omp_rd_13); }
POMP_Parallel_join(&omp_rd_13);
#line 367 "simpar.c"

    // freeing the allocated memory
    free(par);
POMP_Parallel_fork(&omp_rd_14);
#line 370 "simpar.c"
#pragma omp parallel    
{ POMP_Parallel_begin(&omp_rd_14);
POMP_For_enter(&omp_rd_14);
#line 370 "simpar.c"
#pragma omp          for nowait
    for (i = 0; i < ncside; i++){
    	free(cell_mat[i]);
    	free(my_lock_matrix[i]);
    }
POMP_Barrier_enter(&omp_rd_14);
#pragma omp barrier
POMP_Barrier_exit(&omp_rd_14);
POMP_For_exit(&omp_rd_14);
POMP_Parallel_end(&omp_rd_14); }
POMP_Parallel_join(&omp_rd_14);
#line 375 "simpar.c"
    free(cell_mat);
    free(my_lock_matrix);

	return 0;
}
