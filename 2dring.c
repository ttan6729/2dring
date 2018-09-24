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
#include "block.c"

#define _BAD_SOURCES
#define max(a,b) ( (a>=b) ? a:b)
#define min(a,b) ( (a<b) ? a:b)

//tags used for MPI communcations
#define  X_result_Tag 0
#define  Y_result_Tag 1
#define  data_tag 2

char *src = NULL;
char *output_fp = NULL;
/** length of each vector**/
/** length of each interval in pairwise algorithm**/
int interval = 10000;
int thread_number;
int steps;
int block_size;
int	M,N;  // size of each sequence and number of sequences
int x_d=1, y_d;  // dimension of 2d ring
int x_comm, y_comm;
//x and y id of processor in designed 2d processor matrix
int px,py;
int chunk_size = 1;
int myid,numprocs;
struct timeval t0, t1, dt,  t_start, t_end, t_pass;
double computation_time = 0;

MPI_Status stat1,stat2;
MPI_Request req1,req2;

void read_arg(int argc, char **argv)
{
    if(myid == 0 && argc < 3)
    {
        printf("2d Ring all-piras calculation v. 1.0"
               "Usage: mpirun -np <numprocs> 2dring <N> [number of vectors] <M> [length of each vector] <output> [path of output] <src> [path of src] <x_d> [x dimension of 2dring]\n"
               "example: mpirun -np 32 2dring -N 1024 -M 1000 -output correaltion_matrix.txt -src sample_data.txt -x 2\n"
               "or srun -n 16 -d 24 2dring -N 1024 -M 1000 -output correaltion_matrix.txt -src sample_data.txt -x 2\n"
               "notice, currently numprocs should be divisible by xd, and numprocs/xd should be even\n"
               );
    }


	int count = 1;
	char *buffer;
	while( count < argc )
	{
		buffer = argv[count];
		if ( strcmp(buffer,"-N") == 0 )
			N = atoi(argv[++count]);
		else if ( strcmp(buffer,"-M") == 0 )
			M = atoi(argv[++count]);
		else if ( strcmp(buffer,"-t") == 0 )
			thread_number = atoi(argv[++count]);
		else if ( strcmp(buffer,"-output") == 0 )
			output_fp = argv[++count] ;
		else if ( strcmp(buffer,"-src") == 0 )
	   		src = argv[++count] ;
		else if ( strcmp(buffer,"-x") == 0 ) 
			x_d = atoi( argv[++count] );

		count++;
	}

	if ( M == -1 || N == -1 || src == NULL)
	{
		printf("warnning, input paramter is incmoplete\n");
		exit(0);
	}


	y_d = numprocs/x_d;
	if ( (numprocs % x_d) != 0)
	{
		printf("node number can't be divided by x dimeension\n");
		exit(0);
	}
	if(myid == 0)
	{
    	printf("x_d: %d, y_d: %d, numprocs:%d \n",x_d,y_d,numprocs);
    	printf("src: %s, N: %d, M: %d \n",src,N,M);
	}
    
    if(x_d*y_d != numprocs)
    {
        printf("product of matrix dimesnion should equal number of processors\n");
        exit(0);
    }
	return;
}


/** covariance between two vectors**/
double update_covariance(double *x, double *y,int len,double x_mean,double y_mean)
{
    double cov = 0;
    int i;
    for( i = 0; i < len ; i++)
        cov += (x[i]-x_mean)*(y[i]-y_mean);

    return cov;
}
//send data in ring order between processors
void left_ring_send(int substeps,MPI_Comm comm,int rank, int p_id ,double *send,double *recv, int len)
{
    int recevier_id = (p_id + rank -1 - substeps)%rank;
    MPI_Sendrecv(send, len, MPI_DOUBLE, recevier_id, data_tag, recv,len,MPI_DOUBLE,
              (p_id+1+substeps)%rank,data_tag,comm,&stat1);
    return;
}


//exchange data in ring order between processors
void left_ring_exchange(int sweep_range,MPI_Comm comm,int rank, int p_id ,double *exchange_data,double *buffer,int len)
{
    int recevier_id = (p_id + rank -1 - sweep_range)%rank;
    MPI_Sendrecv(exchange_data, len, MPI_DOUBLE, recevier_id,data_tag,
    buffer, len,MPI_DOUBLE, (p_id+1+sweep_range)%rank,data_tag,comm,&stat1);
    printf("exhange, id:%d, recevier_id:%d\n",p_id,recevier_id);

    double *temp = buffer;
    buffer = exchange_data;
    exchange_data = temp;
    return;
}


double correlation(double *x, double *y, int len, double x_mean,double x_var, double y_mean, double y_var)
{
    int i;
    double cov=0,corr =0;
    for(i=0; i < len; i++)
    {
        cov += (x[i]-x_mean)*(y[i]-y_mean);
    }

    cov  = cov/(len-1);
    corr = cov/(sqrt(x_var)*sqrt(y_var));

    if(myid == 5)
    {
        printf("");
    } 
        //printf("checked correlation result: %lf, %lf, var: %lf, %lf , corr is %lf\n",x_mean,y_mean,x_var,y_var,corr);
    
    return corr;
}

void single_block_calculation(Block *b,double *result,int result_count)
{
    int i,j,k,k_m,index = 0;
    #pragma omp parallel for schedule(dynamic,chunk_size) private(i,j,index)
    for(i =0; i < b->vec_num -1; i++)
    {
        for(j = i+1; j < b->vec_num; j++)
        {
            result[result_count+index] = correlation(b->data+i*M,b->data+j*M,M,b->mean[i],b->var[i],b->mean[j],b->var[j]);        
            if(myid == 12)
                printf("checked single result:%lf, index %d\n",result[result_count+index],result_count+index);  
            index++;
        }
    }
}


void pair_block_calculation(Block *b1,Block *b2,double *result, int result_count)
{
    gettimeofday(&t_start, NULL);
    int i,j,k,k_m,index=0;
    int interval_number = M/interval;

    #pragma omp parallel for schedule(dynamic,chunk_size) private(i,j,index)
    for(i = 0; i < b1->vec_num; i++)
    {
        for(j = 0; j < b2->vec_num; j++)
        {
            result[result_count+index] = correlation(b1->data+i*M,b2->data+j*M,M,b1->mean[i],b1->var[i],b2->mean[j],b2->var[j]);
            index++;
        }
    }

    gettimeofday(&t_end, NULL);
    timersub(&t_end, &t_start, &t_pass);
    computation_time += (int)t_pass.tv_sec + ((int)t_pass.tv_usec)/1000000.0;
}



void distributed_work(int px,int py,int substeps,Block *b1,Block *b2, double *intermidate_result)
{
    printf("begin distribtued work (%d,%d)\n",px,py);
    int i,j,k,iter,tid,index;
    double *b1_mean = b1 -> mean;
    double *b1_var = b1 -> var;
    double *b2_mean = b2 -> mean;
    double *b2_var = b2 -> var;

    int result_count = 0;
    int b1_id = py;
    int b2_id = b1_id + px*substeps;

    gettimeofday(&t_start,NULL);
    for(iter =0; iter < substeps-1; iter++)    
    {
        b2_id = (b2_id+1)%y_d;
        pair_block_calculation(b1,b2,intermidate_result, result_count);
    
        if(myid == 0)
            printf("iter %d\n",iter);

        result_count += b1->vec_num * b2->vec_num;
        int receiver_yid = (b1_id + y_d-px*substeps-iter-1-1)%y_d;
        int receiver_id = px*y_d + receiver_yid;
        int sender_id = px*y_d+ (b2_id+1)%y_d;
        MPI_Sendrecv(b1->data,b1->vec_num*M,MPI_DOUBLE,receiver_id,data_tag,
                     b2->data,b2->vec_num*M,MPI_DOUBLE,sender_id,data_tag,MPI_COMM_WORLD,&stat1);
        MPI_Sendrecv(b1->mv,b1->vec_num*2,MPI_DOUBLE,receiver_id,data_tag,
                     b2->mv,b2->vec_num*2,MPI_DOUBLE,sender_id,data_tag,MPI_COMM_WORLD,&stat1);                  
        b2->mean = b2->mv+0;
        b2->var = b2->mv+b2->vec_num;
    }

    if(y_d%2 == 0 && px == x_d-1 && py >= y_d/2)
    {
        single_block_calculation(b1,intermidate_result,result_count);
        result_count += block_size*(block_size-1)/2;
        single_block_calculation(b2,intermidate_result,result_count);

    }
    else
    {
        pair_block_calculation(b1,b2,intermidate_result, result_count);
    }
    return;
}

int main(int argc, char *argv[])
{
    int x,y,i,j,k,l,count=0;
    MPI_Init(NULL,NULL);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);

    read_arg(argc, argv);
    char p_name[1024];

    if(myid == 0)
        puts("warning: this is a version used to test performance of 2d jacobi ring for large number of seqeuence(1*10e6), the resultant correlation matrix is not kept due to its size\n");

    thread_number =  omp_get_num_threads();

    px = myid/y_d;
    py = myid%y_d;
    // contruct two new groups of comunnicators
    MPI_Comm_split(MPI_COMM_WORLD,px,py,&x_comm);
    MPI_Comm_split(MPI_COMM_WORLD,py,px,&y_comm);

    if(numprocs%2 == 1)
        steps = (y_d-1)/2;
    else
        steps = y_d/2;

    if(steps%x_d !=0)
    {
        puts("number of step should be divisible by x dimension of processor group");
        exit(0);
    }

    printf("construct blocks, %d\n",myid);
    int substeps = steps/x_d;
    block_size = N/y_d;

    Block *block1 = build_block(block_size,M);
    Block *block2 = build_block(block_size,M);
    int partial_bs = block_size/x_d;
    Block *partial_block1 = build_block(partial_bs,M);
    
    printf("read rows from src %s, processor %d\n",src,myid);    
    read_rows(src,partial_block1,py*block_size+px*partial_bs);

    //generate_synthetic_data(partial_block1, py*block_size+px*partial_vb+1,partial_vb,M);
    calculate_mean_var(partial_block1);
    MPI_Allgather(partial_block1->data,partial_bs*M,MPI_DOUBLE,block1->data,partial_bs*M,MPI_DOUBLE,y_comm);
    MPI_Allgather(partial_block1->mean,partial_bs,MPI_DOUBLE,block1->mean,partial_bs,MPI_DOUBLE,y_comm);
    MPI_Allgather(partial_block1->var,partial_bs,MPI_DOUBLE,block1->var,partial_bs,MPI_DOUBLE,y_comm);


    //int receiver_id =  py==0 ? (px+1)*y_d-1 : myid-1;
    //int sender_id = px*y_d +  (py+1)%y_d ;
    int receiver_id = (py+ y_d - px*substeps -1)%y_d + px * y_d;
    int b2_id = (py + px*substeps +1)%y_d;    
    int sender_id = b2_id + px*y_d; 

    MPI_Sendrecv(block1->data,block_size*M,MPI_DOUBLE,receiver_id,data_tag,
                 block2->data,block_size*M,MPI_DOUBLE,sender_id,data_tag,MPI_COMM_WORLD,&stat1);
    MPI_Sendrecv(block1->mv,block_size*2,MPI_DOUBLE,receiver_id,data_tag,
                 block2->mv,block_size*2,MPI_DOUBLE,sender_id,data_tag,MPI_COMM_WORLD,&stat1);
    block2->mean = block2->mv+0;
    block2->var = block2->mv+block_size;    

    int result_size = block_size * block_size * substeps;
    double *intermidate_result = malloc(sizeof(double) * result_size);
    if(y_d%2 == 0)
        distributed_work(px,py,substeps,block1,block2,intermidate_result);
    MPI_Barrier(MPI_COMM_WORLD);

    if(numprocs%2 ==1)
    {
        printf("currently only even version is available");
        exit(0);
    }

    if(myid != 0)
    {
        MPI_Send(intermidate_result,result_size,MPI_DOUBLE,0,X_result_Tag,MPI_COMM_WORLD);
    }

    if(myid == 0)
    {
        MPI_Request *reqs = malloc(sizeof(MPI_Request)*numprocs);
        MPI_Status *stats = malloc(sizeof(MPI_Status)*numprocs);

        double **results = malloc(sizeof(double *) * numprocs);
        results[0] = intermidate_result;
        for(i = 1; i < numprocs; i++ )
        {
          results[i] = malloc(sizeof(double) * result_size);
          MPI_Irecv(results[i],result_size,MPI_DOUBLE,i,X_result_Tag,MPI_COMM_WORLD,reqs+i-1);
          printf("receive intermediate result from processor %d\n",i);
        }


        MPI_Waitall(numprocs-1,reqs, stats);
        printf("finish receiving resultat data\n");

        double *correlations = malloc(sizeof(double) * N*(N-1)/2);
        double **correlation_matrix = malloc(sizeof(double *) * (N-1));
        count = 0;
        for(i = 0; i < N-1; i++)
        {
            correlation_matrix[i] = correlations + count;
            count += (N-i-1);
        }
        int block2_id;
        count = 0;
        
        if ( output_fp == NULL )
            output_fp = "correlation_matrix.txt";        

        for(i=0; i < numprocs; i++)
        {
            printf("sort resutlt %d\n",i);
            int i_px = i/y_d;
            int i_py = i%y_d;
            int block1_id = i_py;

            count = 0;
            int checked_substeps = substeps;
            if(y_d%2 == 0 && i_px == x_d-1 && i_py >= y_d/2)
                checked_substeps -= 1;

            for(j=0; j <checked_substeps; j++)
            {
                block2_id = (block1_id + i_px*substeps+j+1)%y_d;
                for(k = 0; k < block_size; k++)
                {
                    for(l = 0; l < block_size; l++)
                    {
                        x = min(block1_id*block_size+k,block2_id*block_size+l);
                        y = max(block1_id*block_size+k,block2_id*block_size+l);
                        y -= (x+1);
                        correlation_matrix[x][y] = results[i][count++];
                    }
                }
            }

            if(y_d%2 == 0 && i_px == x_d-1 && i_py >= y_d/2)
            {
                block2_id = (block1_id + (steps-1)+1)%y_d;              
                for(j=0; j < block_size; j++)
                {
                    for(k=0; k < block_size-1-j; k++)
                    {
                        correlation_matrix[block1_id*block_size+j][k]= results[i][count++];
                    }
                }
                for(j=0; j < block_size; j++)
                {
                    for(k=0; k < block_size-1-j; k++)
                    {                    
                        correlation_matrix[block2_id*block_size+j][k]= results[i][count++];
                    }
                }

            }

        }
        
        //write result to result file
        FILE *fp = fopen(output_fp,"w");
        for( i = 0; i < N-1; i++)
        {
            for( j = 0; j < N-1-i; j++)
                fprintf(fp,"%.3f ",correlation_matrix[i][j]);
            fprintf(fp,"\n");
        }

        /**
        for(i=1; i < numprocs; i++)
            free(results[i]);

        free(correlation_matrix);
        free(correlations);
        free(results);
        free(reqs);
        free(stats); **/
    }    

    printf("finish, %d\n",myid);
    MPI_Barrier(MPI_COMM_WORLD);    
    MPI_Finalize();
//    delete_block(partial_block1);
//    delete_block(block1);
//    delete_block(block2);    
    return 0;

}

