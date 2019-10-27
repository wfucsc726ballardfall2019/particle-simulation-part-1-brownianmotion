#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include "common.h"

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    int navg,nabsavg=0;
    double davg,dmin, absmin=1.0, absavg=0.0;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
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






    ///we need to make some buckets to hold our particles. Use a 2D array of vectors of appropriate size.

    double size = sqrt( 0.0005 * n );
    int buckets_dim = ceil(size/(0.01));
    std::vector<int> bucket_to_particle[buckets_dim][buckets_dim];

    //we are also going to need the inverse; takes partice number to bucket it belongs to. (we do natural linear indexing on the 2D bucket array)
    int particle_to_bucket[n];


    for( int step = 0; step < NSTEPS; step++ )
    {
	 navg = 0;
     davg = 0.0;
	 dmin = 1.0;

        //
        //  compute forces
        // Here we use our buckets, to update the particles
        

        //first, bin the particles
        for ( int i = 0; i < n; i++)
        {   
            particles[i].ax = particles[i].ay = 0;
            int xbucket = floor(particles[i].x/(0.01));
            int ybucket = floor(particles[i].y/(0.01));
            //add to bucket
            bucket_to_particle[xbucket][ybucket].push_back(i);
            // and record bucket number
            particle_to_bucket[i] = xbucket*buckets_dim + ybucket;
        }

        //Now, go through each particle and apply force from particles in neighboring buckets
            for( int i=0; i< n; i++){
                int xbucket = particle_to_bucket[i] / buckets_dim;
                int ybucket = particle_to_bucket[i] % buckets_dim;     
                // really annoying edge-case stuff. We can only access neighbor buckets if we are in the interior
                    int xmin,xmax,ymin,ymax;
                    if(xbucket > 0){xmin = -1;}
                    else{ xmin = 0;}
                    if(xbucket < buckets_dim-1){ xmax = 1;}
                    else{ xmax = 0;}

                    if(ybucket > 0){ ymin = -1;}
                    else{ ymin = 0;}
                    if(ybucket < buckets_dim-1){ ymax = 1;}
                    else{ ymax = 0;}
                    //for all the neighboring buckets
                    for( int nbucketx = xmin; nbucketx <= xmax; nbucketx++){
                        for( int nbuckety = ymin; nbuckety <= ymax; nbuckety++){
                            //for all the particles in those buckets
                            for(int l = 0; l < bucket_to_particle[xbucket+nbucketx][ybucket+nbuckety].size(); l++){
                                    //aply force
                                    apply_force(particles[i],particles[bucket_to_particle[xbucket+nbucketx][ybucket+nbuckety].at(l)],&dmin,&davg,&navg);                            
                            }

                        }
                    }

                
            }
            
        //then go back and empty those buckets
        for( int xbucket = 0; xbucket < buckets_dim; xbucket ++){
                for( int ybucket = 0; ybucket < buckets_dim; ybucket ++){  
                    bucket_to_particle[xbucket][ybucket].clear();
                }
            }


        //from here down, the code is unchanged

        for( int i = 0; i < n; i++ ) 
            move( particles[i] );		

        if( find_option( argc, argv, "-no" ) == -1 )
        {
          //
          // Computing statistical data
          //
          if (navg) {
            absavg +=  davg/navg;
            nabsavg++;
          }
          if (dmin < absmin) absmin = dmin;
		
          //
          //  save if necessary
          //
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
    }
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d, simulation time = %g seconds", n, simulation_time);

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
        fprintf(fsum,"%d %g\n",n,simulation_time);
 
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
