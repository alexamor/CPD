//SIMPAR.C

// Libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>


// Constant variables and functions
#define RND0_1 ((double) random() / ((long long)1<<31))
#define G 6.67408e-11
#define EPSLON 0.0005

#define BLOCK_LOW(id,p,n) ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id) + 1, p, n) - 1)
#define BLOCK_SIZE(id, p, n) (BLOCK_HIGH(id, p, n) - BLOCK_LOW(id, p, n) + 1)
#define BLOCK_OWNER(index, p, n) (((p)*((index) + 1) - 1) / (n))

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


//void scatterPart(particle_t * par, int * counts, int * disp, particle_t * sub_par, int size);


int main(int argc, char *argv[]){
	long seed, ncside;
	long long n_part, n_tstep;
	particle_t *par = NULL, *sub_par = NULL;
	int *counts, *disp;
	int size;
	cell *cell_mat = NULL, *sub_cell_mat;
	int i, j, k, t;
	double F, d2, dx, dy;
	double ax, ay;
	int column, row, adj_column, adj_row;
	double m = 0, sub_m = 0, mx = 0, sub_mx = 0, my = 0, sub_my = 0;
	int id, nprocs;
	double *y_cell, *y_cell_sub, *x_cell, *x_cell_sub, *m_cell, *m_cell_sub;

	// check the correct number of arguments
	if( argc != 5){
		printf("Number of arguments invalid. Please check the way you're running the program, as follows: ./simpar seed ncside n_part n_tstep \n" );
		MPI_Finalize();
		exit(0);
	}

	// initialize MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

	printf("--Inits\n");
	fflush(stdout);


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
		MPI_Finalize();
		exit(0);
	}

	// allocation of size of sub array of particles
	counts = (int *) malloc( sizeof(int) * nprocs);
	if (counts == NULL){
		printf("No memory available for the number of particles required.\n");
		MPI_Finalize();
		exit(0);
	}

	// allocation of offsets of sizes
	disp = (int *) malloc( sizeof(int) * nprocs);
	if (disp == NULL){
		printf("No memory available for the number of particles required.\n");
		MPI_Finalize();
		exit(0);
	}

	if(!id){
		for(i = 0; i < nprocs; i++){
			counts[i] = BLOCK_SIZE(i, nprocs, n_part);
		}

		disp[0] = 0;
		for(i = 1; i < nprocs; i++){
			disp[i] = disp[i-1] + counts[i-1];
		}
	}

	size = BLOCK_SIZE(id, nprocs, n_part);


	// allocation of memory for the particle_t sub array
	sub_par = (particle_t *) malloc( sizeof(particle_t) * size);
	if (par == NULL){
		printf("No memory available for the number of particles required.\n");
		MPI_Finalize();
		exit(0);
	}

	printf("--Memories\n");
	fflush(stdout);



	// creation of particles - initialization of their position, velocity and mass
	if(!id)
	{
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

	}

	////Create MPI type of particle
    const int par_nitems = 7;
    int par_blocklengths[7] = {1,1,1,1,1,1,1,};
    MPI_Datatype par_types[7] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Datatype mpi_par_type;
    MPI_Aint par_offsets[7];

    par_offsets[0] = offsetof(particle_t, x);
    par_offsets[0] = offsetof(particle_t, y);
    par_offsets[0] = offsetof(particle_t, vx);
    par_offsets[0] = offsetof(particle_t, vy);
    par_offsets[0] = offsetof(particle_t, m);
    par_offsets[0] = offsetof(particle_t, fx);
    par_offsets[0] = offsetof(particle_t, fy);

    MPI_Type_create_struct(par_nitems, par_blocklengths, par_offsets, par_types, &mpi_par_type);
    MPI_Type_commit(&mpi_par_type);
    
    //scatterPart(par, counts, disp, sub_par, size);	

    ////Create MPI type of cell
    const int cell_nitems = 3;
    int cell_blocklengths[3] = {1,1,1};
    MPI_Datatype cell_types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Datatype mpi_cell_type;
    MPI_Aint cell_offsets[3];

    cell_offsets[0] = offsetof(cell, x);
	cell_offsets[1] = offsetof(cell, y);
    cell_offsets[2] = offsetof(cell, m);

    MPI_Type_create_struct(cell_nitems, cell_blocklengths, cell_offsets, cell_types, &mpi_cell_type);
    MPI_Type_commit(&mpi_cell_type);

    printf("--Types\n");
	fflush(stdout);

	printf("--size -n%d - %d\n", id, size);
	fflush(stdout);

	MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

    int error_code = MPI_Scatterv(par, counts, disp, mpi_par_type, sub_par, size, mpi_par_type, 0, MPI_COMM_WORLD);
  
	char error_string[BUFSIZ];    
	int length_of_error_string;    
	MPI_Error_string(error_code, error_string, &length_of_error_string);    
	//fprintf(stdout, "%3d: %s\n", id, error_string);

	printf("++ Scatter  %s\n", error_string);
	fflush(stdout);  

    printf("--Scatter -n%d\n", id);
	fflush(stdout);


    // memory allocation for each cell
    cell_mat = (cell *) calloc( ncside * ncside, sizeof(cell));
    if (cell_mat == NULL){
    	printf("No memory available for the number of cells required.\n");
    	MPI_Finalize();
    	exit(0);
    }
    /*
    for(i = 0; i < ncside; i++){

    	// initialization to zero
    	cell_mat[i] = (cell *) calloc( ncside, sizeof(cell));

    	if (cell_mat[i] == NULL){
    		printf("No memory available for the number of cells required.\n");
    		MPI_Finalize();
    		exit(0);
    	}
    }
    */

    sub_cell_mat = (cell *) calloc( ncside * ncside, sizeof(cell));
    if (sub_cell_mat == NULL){
    	printf("No memory available for the number of cells required.\n");
    	MPI_Finalize();
    	exit(0);
    }

    y_cell = (double *) calloc( ncside * ncside, sizeof(double));
    if (y_cell == NULL){
    	printf("No memory available for the number of cells required.\n");
    	MPI_Finalize();
    	exit(0);
    }
    y_cell_sub = (double *) calloc( ncside * ncside, sizeof(double));
    if (y_cell_sub == NULL){
    	printf("No memory available for the number of cells required.\n");
    	MPI_Finalize();
    	exit(0);
    }

    x_cell = (double *) calloc( ncside * ncside, sizeof(double));
    if (x_cell == NULL){
    	printf("No memory available for the number of cells required.\n");
    	MPI_Finalize();
    	exit(0);
    }
    x_cell_sub = (double *) calloc( ncside * ncside, sizeof(double));
    if (x_cell_sub == NULL){
    	printf("No memory available for the number of cells required.\n");
    	MPI_Finalize();
    	exit(0);
    }
    m_cell = (double *) calloc( ncside * ncside, sizeof(double));
    if (m_cell == NULL){
    	printf("No memory available for the number of cells required.\n");
    	MPI_Finalize();
    	exit(0);
    }
    m_cell_sub = (double *) calloc( ncside * ncside, sizeof(double));
    if (m_cell_sub == NULL){
    	printf("No memory available for the number of cells required.\n");
    	MPI_Finalize();
    	exit(0);
    }


    // cycle of time step iterations
    for(t = 0; t < n_tstep; t++){

    	// zeroing the mass centers
    	for(i = 0; i < ncside; i++)
    		for(j = 0; j < ncside; j++){
    			sub_cell_mat[i * ncside + j].x = 0;
    			sub_cell_mat[i * ncside + j].y = 0;
    			sub_cell_mat[i * ncside + j].m = 0;
    			y_cell[i * ncside + j] = 0;
    			x_cell[i * ncside + j] = 0;
    			m_cell[i * ncside + j] = 0;
    		}

		printf("--Zeroes -n%d -i%d\n", id, t);
		fflush(stdout);


    	// calculation of the mass center
		for(i = 0; i < size; i++){

			// get location in grid from the position - truncating the float value
			column = (int) (sub_par[i].x * ncside);
			row = (int) (sub_par[i].y * ncside);

			printf("--row columns -n%d -i%d\n", id, i);
			fflush(stdout);

			y_cell[column * ncside + row] += sub_par[i].m * sub_par[i].y;
			x_cell[column * ncside + row] += sub_par[i].m * sub_par[i].x;
			m_cell[column * ncside + row] += sub_par[i].m;
			/*// average calculated progressively without needing to store every x and y value
			sub_cell_mat[column * ncside + row].y += sub_par[i].m * sub_par[i].y;
			sub_cell_mat[column * ncside + row].x += sub_par[i].m * sub_par[i].x;

			// total mass
			sub_cell_mat[column * ncside + row].m += sub_par[i].m;*/

		}

		printf("--Mass center-n%d -i%d\n", id, t);
		fflush(stdout);

		//Perform redution on matrix of cells calculations
		MPI_Reduce(x_cell, x_cell_sub, ncside * ncside, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(y_cell, y_cell_sub, ncside * ncside, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		error_code = MPI_Reduce(m_cell, m_cell_sub, ncside * ncside, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

		error_string[BUFSIZ];    
		length_of_error_string;    
		MPI_Error_string(error_code, error_string, &length_of_error_string);    
		//fprintf(stdout, "%3d: %s\n", id, error_string);

		printf("++ Reduce  %s\n", error_string);
		fflush(stdout);

		printf("--Reduce-n%d\n", id);
		fflush(stdout);

		if(!id){
			for(i = 0; i < ncside; i++)
				for(j = 0; j < ncside; j++){
					cell_mat[i * ncside + j].x = x_cell_sub[i* ncside + j];
					cell_mat[i * ncside + j].y = y_cell_sub[i* ncside + j];
					cell_mat[i * ncside + j].m = m_cell_sub[i* ncside + j];

					cell_mat[i * ncside + j].x /= cell_mat[i * ncside + j].m;
					cell_mat[i * ncside + j].y /= cell_mat[i * ncside + j].m;
				}			
		}

		printf("--Division-n%d\n", id);
		fflush(stdout);

		//send cell data to everyone
		error_code = MPI_Bcast(cell_mat, ncside* ncside, mpi_cell_type, 0, MPI_COMM_WORLD);

		error_string[BUFSIZ];    
		length_of_error_string;    
		MPI_Error_string(error_code, error_string, &length_of_error_string);    
		//fprintf(stdout, "%3d: %s\n", id, error_string);

		printf("++ Broad %s\n", error_string);
		fflush(stdout);

		printf("--Broadcast-n%d\n", id);
		fflush(stdout);

    	// calculation of the gravitational force
    	for(i = 0; i < size; i++){

    		// get location in grid from the position - truncating the float value
    		column = (int) (sub_par[i].x * ncside);
    		row = (int) (sub_par[i].y * ncside);

    		//printf("column %d   row %d\n", column, row);
    		//printf("x %.2lf y %.2lf\n", par[i].x, par[i].y );

			for(j = -1; j < 2; j++){

				//printf("j = %d\n", j);

				adj_column = column + j;

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
					dx = cell_mat[adj_column * ncside + adj_row].x - sub_par[i].x;
					dy = cell_mat[adj_column * ncside + adj_row].y - sub_par[i].y;

					//printf("mat.y %lf  par.y %lf\n", cell_mat[adj_column * ncside + adj_row].y, par[i].y);

					// calculate distances when the cells are out of borders
					if(j == -1 && adj_column == ncside - 1)
						dx =  cell_mat[adj_column * ncside + adj_row].x - sub_par[i].x - 1;
					else if(j == 1 && adj_column == 0)
						dx = 1 + (cell_mat[adj_column * ncside + adj_row].x - sub_par[i].x);

					if(k == 1 && adj_row == 0)
						dy = cell_mat[adj_column * ncside + adj_row].y - sub_par[i].y + 1;
					else if(k == -1 && adj_row == ncside - 1)
						dy = 1 + (cell_mat[adj_column * ncside + adj_row].y - sub_par[i].y);


					//printf("dx %lf     dy %lf\n", dx, dy);


					// calculate distance and gravitational force
					d2 = dx*dx + dy*dy ;

					// check to see if distance is minimal
					if (d2 < EPSLON)
						F = 0;
					else
						F = G * sub_par[i].m * cell_mat[adj_column * ncside + adj_row].m / d2;


					//printf("%lf\n", d2 );

					// get cartesian components of the gravitational force
					if ( j == -1 && k == -1 ){
						sub_par[i].fx = 0;
						sub_par[i].fy = 0;
					}

					// calculate force
					if(d2 != 0 && !isnan(d2)){

						sub_par[i].fx += (F*dx)/sqrt(d2);
						sub_par[i].fy += (F*dy)/sqrt(d2);


					}

				}

			}


   		}

   		// get new positions
    	for(i = 0; i < size; i++){

			// get acceleration
			ax = sub_par[i].fx / sub_par[i].m;
			ay = sub_par[i].fy / sub_par[i].m;

			// position
			sub_par[i].x += sub_par[i].vx + ax;
			sub_par[i].y += sub_par[i].vy + ay;

			if(sub_par[i].x < 0){
				sub_par[i].x += 1;
			}
			else if(sub_par[i].x > 1){
				sub_par[i].x -= 1;
			}

			if(sub_par[i].y < 0){
				sub_par[i].y += 1;
			}
			else if(sub_par[i].y > 1){
				sub_par[i].y -= 1;
			}


			// velocity
			sub_par[i].vx += ax;
			sub_par[i].vy += ay;



			//printf("fx %lf   fy %lf\n", par[i].fx, par[i].fy);
			//printf("ax %lf ay %lf \n", ax, ay);
			fflush(stdout);


    	}

    	//Wait for everyone before starting new iteration
    	error_code = MPI_Barrier(MPI_COMM_WORLD);

    	error_string[BUFSIZ];    
		length_of_error_string;    
		MPI_Error_string(error_code, error_string, &length_of_error_string);    
		//fprintf(stdout, "%3d: %s\n", id, error_string);

		printf("++ Barrier -n%d %s\n", id, error_string);
		fflush(stdout);	

    }





    // calculate final mass center
    for(i = 0; i < size; i++){

    	    // average calculated progressively without needing to store every x and y value
    		sub_mx += sub_par[i].m * sub_par[i].x;
    		sub_my += sub_par[i].m * sub_par[i].y;

    		//printf("mx %lf my %lf\n", par[i].x, par[i].y );

    		// total mass
    		sub_m += sub_par[i].m;

    }

    MPI_Reduce(&sub_mx, &mx, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&sub_my, &my, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&sub_m, &m, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	//Note: Altrough is faster to average the cells the result ends up poorly rounded, failing the correct answer for +-0.01
    /*
    for(i = 0; i < ncside; i++)
    {
    	for(j = 0; j < ncside; j++)
    	{
    		mx += cell_mat[i * ncside + j].x * cell_mat[i * ncside + j].m;
    		my += cell_mat[i * ncside + j].y * cell_mat[i * ncside + j].m;

    		m += cell_mat[i * ncside + j].m;
    	}
    }
    */
	if(!id){
	    mx = mx/m;
	    my = my/m;

	    // output
	    printf("%.2lf %.2lf\n", par[0].x, par[0].y );
	    printf("%.2lf %.2lf\n", mx, my );	
	}


    // freeing the allocated memory
    free(par);
    free(sub_par);
    free(counts);
    free(disp);
    free(y_cell);
    free(y_cell_sub);
    free(x_cell);
    free(x_cell_sub);
    free(m_cell);
    free(m_cell_sub);
    /*
    for (i = 0; i < ncside; i++){
    	free(cell_mat[i]);
    }
    */
    free(cell_mat);
    free(sub_cell_mat);

    MPI_Finalize();

	return 0;
}
/*
void scatterPart(particle_t * par, int * counts, int * disp, particle_t * sub_par, int size){
	MPI_Scatterv(par.x, counts, disp, MPI_DOUBLE, sub_par.x, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatterv(par.y, counts, disp, MPI_DOUBLE, sub_par.y, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatterv(par.vx, counts, disp, MPI_DOUBLE, sub_par.vx, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatterv(par.vy, counts, disp, MPI_DOUBLE, sub_par.vy, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatterv(par.m, counts, disp, MPI_DOUBLE, sub_par.m, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}
*/
