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

    double size = sqrt( 0.0005 * n );
    int buckets_dim = ceil(size/(2*0.01));
    std::vector<int> buckets[buckets_dim][buckets_dim];


    #pragma omp parallel  private(dmin)
    {
    numthreads = omp_get_num_threads();
    for( int step = 0; step < NSTEPS; step++ )
    {
        navg = 0;
        davg = 0.0;
	    dmin = 1.0;
        double temp = read_timer();
        //
        //  compute all forces
        //
        // #pragma omp for reduction (+:navg) reduction(+:davg)
        // for( int i = 0; i < n; i++ )
        // {
        //     particles[i].ax = particles[i].ay = 0;
        //     for (int j = 0; j < n; j++ )
        //         apply_force( particles[i], particles[j],&dmin,&davg,&navg);
        // }
        


        #pragma omp for
        for ( int i = 0; i < n; i++)
        {   
            
            particles[i].ax = particles[i].ay = 0;
            int xbucket = floor(particles[i].x/(2*0.01));
            int ybucket = floor(particles[i].y/(2*0.01));
            #pragma omp critical
            buckets[xbucket][ybucket].push_back(i);
           
        }

 
        #pragma omp for collapse(2) reduction (+:navg) reduction(+:davg)
        for( int xbucket = 0; xbucket < buckets_dim; xbucket ++){
            for( int ybucket = 0; ybucket < buckets_dim; ybucket ++){
                for(int k = 0; k < buckets[xbucket][ybucket].size(); k++){
               
                    int i = buckets[xbucket][ybucket].at(k);
                    int xmin,xmax,ymin,ymax;
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

                            for(int l = 0; l < buckets[xbucket+nbucketx][ybucket+nbuckety].size(); l++){
                                if(l != k || nbucketx!= 0 || nbuckety != 0){
                                    
                                    
                                    
                                    apply_force(particles[i],particles[buckets[xbucket+nbucketx][ybucket+nbuckety].at(l)],&dmin,&davg,&navg);
                                }
                                
                            }

                        }
                    }

                }
            }
            
        }

        #pragma omp for collapse(2)
        for( int xbucket = 0; xbucket < buckets_dim; xbucket ++){
                for( int ybucket = 0; ybucket < buckets_dim; ybucket ++){  
                    buckets[xbucket][ybucket].clear();
                }
            }



		
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

