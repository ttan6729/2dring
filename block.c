#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/types.h>

//each block can be considered information of a set of vectors
typedef struct Blocks
{
	int vec_num;
	int vec_len;
	double *data;
	double *mv;
	double *mean;
	double *var;
} Block;

void delete_block(Block *block)
{
	free(block->data);
	free(block->mv);
	free(block);
	return;
}

Block *build_block(int num,int len)
{
	Block *block = (Block *)malloc( sizeof(Block) );
	block-> vec_num = num;
	block-> vec_len = len;
	block-> data = malloc( sizeof(double) * num * len);
	block-> mv = malloc( 2 * sizeof(double) * num );
	block-> mean = block->mv+0;
	block-> var = block->mv+num;
	return block;
}

void read_rows(char *src, Block *block, int start_id)
{
	int len = block->vec_len;
	int num = block->vec_num;
	int i,j;
	double buff;
	FILE *file;
	file = fopen(src,"r");
	for( i = 0; i < start_id * len; i++)
	{
		if (!fscanf(file, "%lf", &buff))
			break;
	}

	for ( i = 0; i < len * num; i++)
	{
		if(!fscanf(file,"%lf", (block->data)+i))
			break;
	}
/**
    printf("mean of block %d from process %d:",py,myid);
    print_vector(block1_mean,block_size);
    printf("var of block %d:",py);
    print_vector(block1_variance,block_size);
**/	
	fclose(file);
	return;
}

/**calculate mean and variance of given vector with given length**/
void vector_mean_variance(double *vector,int len, double *i_mean,double *variance)
{
    int i;
    double var = 0,mean = 0;
    for(i=0; i < len; i++)
        mean += vector[i];

    mean = mean/len;

    for(i = 0; i < len; i ++)
        var += (vector[i]-mean)*(vector[i]-mean);

    var = var/(len-1);

    (*i_mean) = mean;
    *variance = var;
    return;
}

void calculate_mean_var(Block *block)
{
	int len = block->vec_len;
	int i;
	for (i = 0; i < block->vec_num; i++)
		vector_mean_variance( block->data + i*len,len,block->mean+i, block->var+i);

	return;
}



