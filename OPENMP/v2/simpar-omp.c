//SIMPAR.C

// Libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>



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
    double aux_atan = 0;



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



#pragma omp parallel for
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

    	// zeroing the mass centers and initating locks
    #pragma omp parallel for private(i,j)
    	for(i = 0; i < ncside; i++)
    		for(j = 0; j < ncside; j++){
    			cell_mat[i][j].x = 0;
    			cell_mat[i][j].y = 0;
    			cell_mat[i][j].m = 0;
    		}

    	// calculation of the mass center
    #pragma omp parallel for private(i,column,row)
    	for(i = 0; i < n_part; i++){

    		// get location in grid from the position - truncating the float value
    		column = (int) (par[i].x * ncside);
    		row = (int) (par[i].y * ncside);



    		// average calculated progressively without needing to store every x and y value
		#pragma omp atomic
    		cell_mat[column][row].y += par[i].m * par[i].y;
		#pragma omp atomic
				cell_mat[column][row].x += par[i].m * par[i].x;

    		// total mass
		#pragma omp atomic
    		cell_mat[column][row].m += par[i].m;

    	}

		#pragma omp parallel for private(i, j)
			for(i = 0; i < ncside; i++)
				for(j = 0; j < ncside; j++){
					cell_mat[i][j].x /= cell_mat[i][j].m;
					cell_mat[i][j].y /= cell_mat[i][j].m;
				}

    	// calculation of the gravitational force
	#pragma omp parallel for private(i, j, k, column, row, adj_column, adj_row, dx, dy, d2, F, aux_atan)
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

					//printf("d2 %lf\n", d2 );
                    //printf("m %lf\n", cell_mat[adj_column][adj_row].m);

					// get cartesian components of the gravitational force
					if ( j == -1 && k == -1 ){
						par[i].fx = 0;
						par[i].fy = 0;
					}

					//printf("atan  %lf, dx %lf, dy %lf\n", atan(dy/dx), dx, dy);
                    aux_atan = atan(dy/dx);


                    if (!isnan(aux_atan)){

                        // calculate force
                        par[i].fx += F*cos(aux_atan);
                        par[i].fy += F*sin(aux_atan);

                    }



					//printf("fx %lf   fy %lf\n", par[i].fx, par[i].fy);
				}

			}


    	}

    	// calculation of new velocity and position
#pragma omp parallel for private(i, ax, ay)
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

    	//printf("Finished iteration: %d\n", t);

    }


    // calculate final mass center
#pragma omp parallel for reduction(+:mx,my,m)
    for(i = 0; i < n_part; i++){

    	    // average calculated progressively without needing to store every x and y value
    		mx += par[i].m * par[i].x;
    		my += par[i].m * par[i].y;

    		//printf("mx %lf my %lf\n", par[i].x, par[i].y );

    		// total mass
    		m += par[i].m;

    }

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

    // freeing the allocated memory
    free(par);
#pragma omp parallel for
    for (i = 0; i < ncside; i++){
    	free(cell_mat[i]);
    }
    free(cell_mat);

	return 0;
}
