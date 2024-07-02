// MIT License
// 
// Copyright (c) 2023 Trevor Bakker 
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <getopt.h>
#include <stdint.h>
#include "utility.h"
#include "star.h"
#include "float.h"

#include <pthread.h>

#define NUM_STARS 30000 
#define MAX_LINE 1024
#define DELIMITER " \t\n"

struct Star star_array[ NUM_STARS ];
uint8_t   (*distance_calculated)[NUM_STARS];

double  min  = FLT_MAX;
double  max  = FLT_MIN;
//mutex for locking critical regions
pthread_mutex_t mutex;
pthread_mutex_t mutex2;
//index to find the size each thread will get for the star_array
int thread_index=0;
//assume 1 thread is -t # is not used
int num_threads=1;

void showHelp()
{
  printf("Use: findAngular [options]\n");
  printf("Where options are:\n");
  printf("-t          Number of threads to use\n");
  printf("-h          Show this help\n");
}

// 
// Embarassingly inefficient, intentionally bad method
// to calculate all entries one another to determine the
// average angular separation between any two stars 
float determineAverageAngularDistance( struct Star arr[], int start, int end)
{
  double mean = 0;

  uint32_t i, j;
  uint64_t count = 0;
  
  for (i = start; i < end; i++)
  {
    for (j = 0; j < NUM_STARS; j++)
    {
      if( i!=j && distance_calculated[i][j] == 0 )
      {
        double distance = calculateAngularDistance( arr[i].RightAscension, arr[j].Declination, arr[j].RightAscension, arr[j].Declination ) ;
        distance_calculated[i][j] = 1;
        distance_calculated[j][i] = 1;
pthread_mutex_lock(&mutex);
        count++;

        if( min > distance )
        {
          min = distance;
        }

        if( max < distance )
        {
          max = distance;
        }
        mean = mean + (distance-mean)/count;
pthread_mutex_unlock(&mutex);
      }
    }
  }
  return mean;
}

void *thread_main(void* arg)
{
  int start=0;
  int end=0;
  struct Star *starArr=(struct Star *)arg;
  //will hold the value of distance but wont be accurate since it will keep the most recent value of last thread
  double* distance = malloc(sizeof(double));
  //holds the size that the star array will be divided depemding on the number of threads
  int arr_size = NUM_STARS / num_threads;
  pthread_mutex_lock(&mutex2);
  {
    start=thread_index * arr_size;
    end=(thread_index+1) * arr_size;
    thread_index++;
  }
  pthread_mutex_unlock(&mutex2);
  *distance=determineAverageAngularDistance( starArr , start, end);
  pthread_exit((void*)distance);
}

int main( int argc, char * argv[] )
{
  FILE *fp;
  uint32_t star_count = 0;

  uint32_t n;

  distance_calculated = malloc(sizeof(uint8_t[NUM_STARS][NUM_STARS]));

  if( distance_calculated == NULL )
  {
    uint64_t num_stars = NUM_STARS;
    uint64_t size = num_stars * num_stars * sizeof(uint8_t);
    printf("Could not allocate %ld bytes\n", size);
    exit( EXIT_FAILURE );
  }

  int i, j;
  // default every thing to 0 so we calculated the distance.
  // This is really inefficient and should be replace by a memset
  for (i = 0; i < NUM_STARS; i++)
  {
    for (j = 0; j < NUM_STARS; j++)
    {
      distance_calculated[i][j] = 0;
    }
  }

  for( n = 1; n < argc; n++ )          
  {
    if( strcmp(argv[n], "-help" ) == 0 )
    {
      showHelp();
      exit(0);
    }
    else if( strcmp(argv[n], "-t" ) == 0 && argv[n+1] != NULL )
    {
      //printf("\nusing %d threads\n",atoi(argv[n+1]));
      num_threads=atoi(argv[n+1]);
    }
  }

  fp = fopen( "data/tycho-trimmed.csv", "r" );

  if( fp == NULL )
  {
    printf("ERROR: Unable to open the file data/tycho-trimmed.csv\n");
    exit(1);
  }

  char line[MAX_LINE];
  while (fgets(line, 1024, fp))
  {
    uint32_t column = 0;

    char* tok;
    for (tok = strtok(line, " ");
            tok && *tok;
            tok = strtok(NULL, " "))
    {
       switch( column )
       {
          case 0:
              star_array[star_count].ID = atoi(tok);
              break;
       
          case 1:
              star_array[star_count].RightAscension = atof(tok);
              break;
       
          case 2:
              star_array[star_count].Declination = atof(tok);
              break;

          default: 
             printf("ERROR: line %d had more than 3 columns\n", star_count );
             exit(1);
             break;
       }
       column++;
    }
    star_count++;
  }
  fclose(fp);
  
  printf("%d records read\n", star_count );

  //start of multithreading
  //vars for the method of timing the program
  long start_t,end_t;
  double total_t;
  //var to hold the mean of the distances (this is not how it would be actually calculated)
  double* distance;
  //array of threads of which will use the thread main function in order to find they part of the stars arr and then call the deterMine.. function
  pthread_t tid[num_threads];

  int index;
  start_t=clock();
  for(index=0;index<num_threads;index++)
  {
    pthread_create(&tid[index],NULL, thread_main,(void*)star_array);
  }
  for(index=0;index<num_threads;index++)
  {
    pthread_join(tid[index],(void**)&distance);
  }
  end_t=clock();

  // start_t=clock();
  // // Find the average angular distance in the most inefficient way possible
  // double distance =  determineAverageAngularDistance( star_array );
  // end_t=clock();

  printf("Average distance found is %lf\n", *distance );
  printf("Minimum distance found is %lf\n", min );
  printf("Maximum distance found is %lf\n", max );
  // //printing the time taken to calc distance
  total_t=(double)(end_t-start_t)/CLOCKS_PER_SEC;
  printf("Time taken: %lf seconds\n",total_t);

  return 0;
}

