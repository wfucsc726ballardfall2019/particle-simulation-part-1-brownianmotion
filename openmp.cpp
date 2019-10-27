#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include "common.h"
#include "omp.h"

//
//  benchmarking program
//
int main( int argc, char **argv )
{   
    int navg,nabsavg=0,numthreads; 
    double dmin,absmin=1.0,davg,absavg=0.0;
	
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" ); 
        printf( "-no turns off all correctness checks and particle output\n");   
        return 0;
    }

    int n = read_int( argc, argv, "-n", 1000 );
    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );

    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;      

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );


    //once again, we need our buckets

    double size = sqrt( 0.0005 * n );
    int buckets_dim = ceil(size/(0.01));
    std::vector<int> bucket_to_particle[buckets_dim][buckets_dim];
    //and our inverse operator

    int particle_to_bucket[n];


    //recruit of team of (somehow slow!?!!?) threads
    #pragma omp parallel  private(dmin)
    {
    numthreads = omp_get_num_threads();
    for( int step = 0; step < NSTEPS; step++ )
    {
        navg = 0;
        davg = 0.0;
	    dmin = 1.0;
        int xbucket;
        int ybucket;
        int xmin,xmax,ymin,ymax;
        

        //we can completely parallelize the binning process...
        #pragma omp for
        for ( int i = 0; i < n; i++)
        {   
            particles[i].ax = particles[i].ay = 0;
            xbucket = floor(particles[i].x/(0.01));
            ybucket = floor(particles[i].y/(0.01));
            particle_to_bucket[i] = xbucket*buckets_dim + ybucket;
            //except when we write to the linked lists; need to be careful not to double write.
            #pragma omp critical
            bucket_to_particle[xbucket][ybucket].push_back(i);
        }   

        //we can also completely parallelize the update loop; they are all independent
        #pragma omp for reduction (+:navg) reduction(+:davg)
        for( int i=0; i< n; i++){
                xbucket = particle_to_bucket[i] / buckets_dim;
                ybucket = particle_to_bucket[i] % buckets_dim;     
                    
                    if(xbucket > 0){xmin = -1;}
                    else{ xmin = 0;}
                    if(xbucket < buckets_dim-1){ xmax = 1;}
                    else{ xmax = 0;}

                    if(ybucket > 0){ ymin = -1;}
                    else{ ymin = 0;}
                    if(ybucket < buckets_dim-1){ ymax = 1;}
                    else{ ymax = 0;}

                    for( int nbucketx = xmin; nbucketx <= xmax; nbucketx++){
                        for( int nbuckety = ymin; nbuckety <= ymax; nbuckety++){

                            for(int l = 0; l < bucket_to_particle[xbucket+nbucketx][ybucket+nbuckety].size(); l++){
                                    
                                    apply_force(particles[i],particles[bucket_to_particle[xbucket+nbucketx][ybucket+nbuckety].at(l)],&dmin,&davg,&navg);                            
                            }

                        }
                    }


            }

        //and the bucket-emptying is parallelizable as well.
        #pragma omp for
        for( int bucket_id=0; bucket_id < buckets_dim * buckets_dim; bucket_id ++){
                xbucket = bucket_id / buckets_dim;
                ybucket = bucket_id % buckets_dim;
                bucket_to_particle[xbucket][ybucket].clear();
            }



        //from here down, the code is unchanged.
        //where is the slowdown???? It still runs, just not fast...


		
        //
        //  move particles
        //
        #pragma omp for
        for( int i = 0; i < n; i++ ) 
            move( particles[i] );
  
        if( find_option( argc, argv, "-no" ) == -1 ) 
        {
          //
          //  compute statistical data
          //
          #pragma omp master
          if (navg) { 
            absavg += davg/navg;
            nabsavg++;
          }

          #pragma omp critical
	  if (dmin < absmin) absmin = dmin; 
		
          //
          //  save if necessary
          //
          #pragma omp master
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
    }
}



    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d,threads = %d, simulation time = %g seconds", n,numthreads, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
      if (nabsavg) absavg /= nabsavg;
    // 
    //  -The minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation where particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //
    printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
    if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
    if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");
    
    //
    // Printing summary data
    //
    if( fsum)
        fprintf(fsum,"%d %d %g\n",n,numthreads,simulation_time);

    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );

    free( particles );

    if( fsave )
        fclose( fsave );
    
    return 0;
}

