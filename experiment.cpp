#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
/*data for processing
steps as genes
| 1 | 2 | 3 | 4 |
each gene has its own mutation points available for chosing
muta[Ngene][Nmuta]
*/
typedef struct {
	int N_gene;//number of experimental steps N_gene=4
	int N_method;//number of method in each step N_method=16
	int N_try;//number of trying in one batch N_try=10 ?
	int N_iteration;//number of batches for history check
	int N_iteration_max;//number of batches for history check
	int ***gene;//content of all steps. size=N_iter*N_steps * N_try
	double** fitness;//fitness for each try: size=N_iter*N_try
	int **muta;//avilable operation for each step max size= N_method*N_gene
	int *N_muta;//effective number for muta
	int mutation_size_max;//changes for each optimization
	int **gene_label;//label for gene samples which is new [0] or old [iter_num]
	double survivor_factor;
	int N_survive;//survived sample in each round
}data;
//initialize data structure 

void init_data(data &d)
{
	FILE* fp_input = fopen("parameter.dat", "r");
	//initialize constant
	d.N_gene = 4;//total steps in experiment[0~3]
	d.N_method = 8;//maximum available methods/states for each gene site
	d.N_try = 10;//parallel try for each iteration (1~10)
	d.N_iteration = 1;//current iteration number starting from 1
	d.N_iteration_max = 10;//max iteration number for trying
	d.mutation_size_max=2;//step size for mutation
	d.survivor_factor = 0.2;//
	char buff[1024];
	fgets(buff, 1024, fp_input);//title
	fgets(buff, 1024, fp_input);//data
	sscanf(buff, "%d%d%d%d%lf", &d.N_gene, &d.N_try, &d.N_iteration_max, &d.mutation_size_max, &d.survivor_factor);
	fclose(fp_input);
	d.N_survive=int(d.N_try * d.survivor_factor);
	//initialize array
	//d.muta = new int[d.N_method];
	d.gene = new int** [d.N_try+1];
	for (int i = 0; i < d.N_iteration_max+1; i++)
	{
		d.gene[i] = new int* [d.N_try+1];
		for (int j = 0; j < d.N_try+1; j++)
		{
			d.gene[i][j] = new int[d.N_gene];
		}
	}

	d.fitness = new double* [d.N_iteration_max+1];
	for (int i = 0; i < d.N_iteration_max+1; i++)
	{
		d.fitness[i] = new double[d.N_try+1];
	}

	d.gene_label = new int* [d.N_iteration_max+1];
	for (int i = 0; i < d.N_iteration_max+1; i++)
	{
		d.gene_label[i] = new int[d.N_try+1];
	}

	d.muta = new int* [d.N_gene];
	for (int i = 0; i < d.N_gene; i++)
	{
		d.muta[i] = new int[d.N_method];
	}
	d.N_muta = new int[d.N_gene];

	//initialize value for some arrays
	for (int i = 0; i < d.N_iteration_max+1; i++)
	{
		for (int j = 0; j < d.N_try; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				d.gene[i][j][k] = 0;
			}
			d.fitness[i][j] = 0;
			d.gene_label[i][j] = 0;
		}
	}
	//initialize muta list
	d.N_muta[0] = 5; d.N_muta[1] = 8;
	d.N_muta[2] = 8; d.N_muta[3] = 5;//?bug
	for (int i = 0; i < d.N_muta[0]; i++) { d.muta[0][i] = 3 + i; }
	for (int i = 0; i < d.N_muta[1]; i++) { d.muta[1][i] = i; }
	for (int i = 0; i < d.N_muta[2]; i++) { d.muta[2][i] = i; }
	for (int i = 0; i < d.N_muta[3]; i++) { d.muta[3][i] = i; }
}
int validate(data& d, int* gene)
{
	//basic requirement : N_method methods form a loop for chosing from
	//if no restriction condition:
	//more restrictions could be added here:
	//r(1) first step only allow for method 3 4 5 6 7
	//r(2) last step only allow for method 1 2 3 4
	if (gene[0] < 0 || gene[0] > 7) return -1;
	if (gene[1] < 0 || gene[1] > 7) return -1;
	if (gene[2] < 0 || gene[2] > 7) return -1;
	if (gene[3] < 0 || gene[3] > 7) return -1;

	if (gene[0] < 3) return -1;
	if (gene[3] < 1 || gene[3] > 4 ) return -1;
	if (gene[1] == 0 && gene[2] != 0) return -1;
	

	//requirement fulfilled, now check history record.
	int found_in = 0;//if found_in == -1 then invalid direct return the value
	for (int i = 1; i <= d.N_iteration-1; i++)
	{
		for (int j = 0; j < d.N_try; j++)
		{
			if (d.gene[i][j][0] == gene[0] &&
				d.gene[i][j][1] == gene[1] &&
				d.gene[i][j][2] == gene[2] &&
				d.gene[i][j][3] == gene[3])
			{
				found_in = i;
				return found_in;
			}
		}
	}
	return found_in;
}
void gen_random(data& d)
{
	//generate sample based on random euclidean distance
	int gene[4];
	int nn = 1;
	while (nn <= d.N_try)
	{
		for (int i = 0; i < d.N_gene; i++)
		{
			int r = rand() % d.N_gene;//0 1 2 3 : gene site
			int m = rand() % d.N_method;
			gene[i] = m;
		}
		int label = validate(d, gene);
		if (label !=-1)//validated
		{
			for (int i = 0; i < 4; i++)
			{
				d.gene[1][nn][i] = gene[i];
			}
			d.gene_label[1][nn] = label;
			nn++;
		}
	}
}

//output next experimental gene selection
void output_new_gene(data& d)
{
	char nm[128];
	sprintf(nm, "gene_%d.dat", d.N_iteration);
	FILE* fp = fopen(nm, "w");
	for (int i = 1; i <= d.N_try; i++)
	{
		fprintf(fp, "%d %d %d %d %d\n", d.gene[d.N_iteration][i][0], d.gene[d.N_iteration][i][1],
			d.gene[d.N_iteration][i][2], d.gene[d.N_iteration][i][3], d.gene_label[d.N_iteration][i]);
	}
	fclose(fp);

}
void output_fitness(data& d)
{
	char nm[128];
	sprintf(nm, "fitness_%d.dat", d.N_iteration);
	FILE* fp = fopen(nm, "w");
	for (int i = 1; i <= d.N_try; i++)
	{
		fprintf(fp, "%lf\n", d.fitness[d.N_iteration][i]);
	}
	fclose(fp);
}
//update fitness from input
void update_fitness(data &d)
{
	FILE* fp_fitness = fopen("fitness.dat", "r");
	char buffer[128];
	for (int i = 1; i <= d.N_try; i++)
	{
		fgets(buffer, 128, fp_fitness);
		sscanf(buffer, "%lf", &d.fitness[d.N_iteration][i]);
	}
	output_fitness(d);//record to history
	//sort genes via fitness array:
	for (int i = 1; i <= d.N_try-1; i++)
	{
		for (int j = i+1; j <= d.N_try; j++)
		{
			if (d.fitness[d.N_iteration][i] < d.fitness[d.N_iteration][j])
			{
				//exchange d.fitness[i] and d.fitness[j]
				double tempf = d.fitness[d.N_iteration][i];
				d.fitness[d.N_iteration][i] = d.fitness[d.N_iteration][j];
				d.fitness[d.N_iteration][j] = tempf;


				//exchange d.gene[i] and d.gene[j]
				for (int k = 0; k < 4; k++)
				{
					int tempi = d.gene[d.N_iteration][i][k];
					d.gene[d.N_iteration][i][k] = d.gene[d.N_iteration][j][k];
					d.gene[d.N_iteration][j][k] = tempi;
				}
			}
		}
	}
}

//based on fitness results generate new genes
//also check the mutation history where the fitness is already avaiable
//mutation must fulfill requirement.

void gen_mutation(data &d)
{
	//generate mutations based on the survived genes (50% of all samples in previous iteration) 
	int gene1[4];//mutated gene
	int nn = 1;
	//based on survived samples, gen next generation via crossgene and mutation algorithm:
	while (nn<=d.N_try)
	{
		if (nn < d.N_survive)
		{
			for (int m = 0; m < d.N_gene; m++)
			{
				gene1[m] = d.gene[d.N_iteration - 1][nn][m];
			}
			int label = validate(d, gene1);
			if (label!=-1)
			{
				for (int i = 0; i < d.N_gene; i++)
				{
					d.gene[d.N_iteration][nn][i] = gene1[i];
				}
				d.gene_label[d.N_iteration][nn] = label;
				nn++;
			}
		}
		else
		{
			//50% for crossgene algorithm, 50% for single site mutation:
			int r_operation = rand() % 2;
			if (r_operation==0)//do crossgene algorithm: 
			{
				int r_gene1 = (rand() % d.N_survive) + 1;//id starts from 1
				int r_gene2 = (rand() % d.N_survive) + 1;//id starts from 1
				int cnt = 0;
				while (r_gene1 == r_gene2)//make sure select two different genes
				{
					r_gene1 = (rand() % d.N_survive) + 1;//id starts from 1
					r_gene2 = (rand() % d.N_survive) + 1;//id starts from 1
					cnt++;
					if (cnt > 10000) {
						printf("it seems like all genes are pretty much identical already, you could stop the calculation now");
						exit(0);
					}
					
				}
				//cross gene between two genes
				int r_dir = rand() % 2;//determine direction
				int r_pos = rand() % (d.N_gene-1);//deteemine gene cutting position
				if (r_dir==0)
				{
					for (int m = 0; m < d.N_gene; m++)
					{
						if (m<r_pos)
						{
							gene1[m] = d.gene[d.N_iteration - 1][r_gene1][m];
						}
						else 
						{
							gene1[m] = d.gene[d.N_iteration - 1][r_gene2][m];
						}
					}
				}
				else 
				{
					for (int m = 0; m < d.N_gene; m++)
					{
						if (m < r_pos)
						{
							gene1[m] = d.gene[d.N_iteration - 1][r_gene2][m];
						}
						else
						{
							gene1[m] = d.gene[d.N_iteration - 1][r_gene1][m];
						}
					}
				}
			}
			else {//do single site mutation algorithm:
				int r_gene =(rand() % d.N_survive)+1;//id starts from 1
				int r_pos = rand() % d.N_gene;//0 1 2 3 : gene site
				int r_muta = rand() % (d.mutation_size_max * 2 + 1) - d.mutation_size_max;//-m ~ +m offset
				gene1[0] = d.gene[d.N_iteration - 1][r_gene][0];
				gene1[1] = d.gene[d.N_iteration - 1][r_gene][1];
				gene1[2] = d.gene[d.N_iteration - 1][r_gene][2];
				gene1[3] = d.gene[d.N_iteration - 1][r_gene][3];
				gene1[r_pos] = gene1[r_pos] + r_muta;
			}
			int label = validate(d, gene1);
			if (label!=-1)
			{
				for (int i = 0; i < d.N_gene; i++)
				{
					d.gene[d.N_iteration][nn][i] = gene1[i];
				}
				d.gene_label[d.N_iteration][nn] = label;
				nn++;
			}
		}
	}
}

void restart(data& d, int num)
{
	char nm[128];
	char buff[1024];
	FILE* fp;
	//history trajectory start from 0 to num-1
	for (int i = 1; i <= num; i++)
	{
		sprintf(nm, "gene_%d.dat", i);
		fp = fopen(nm, "r");
		for (int j = 1; j <= d.N_try; j++)
		{
			fgets(buff, 1024, fp);
			sscanf(buff, "%d%d%d%d%d", &d.gene[i][j][0], &d.gene[i][j][1], 
				&d.gene[i][j][2], &d.gene[i][j][3],&d.gene_label[i][j]);
		}
		fclose(fp);
	}

	for (int i = 1; i <= num; i++)
	{
		sprintf(nm, "fitness_%d.dat", i);
		fp = fopen(nm, "r");
		for (int j = 1; j <= d.N_try; j++)
		{
			fgets(buff, 1024, fp);
			sscanf(buff, "%lf", &d.fitness[i][j]);
		}
		fclose(fp);
	}
	d.N_iteration = num;
}
void precalculation(data& d)
{
	int is_restart = false;
	printf("Is this run start from previous calculation? [1]yes[0]no");
	scanf("%d", &is_restart);
	init_data(d);
	if (is_restart == 0)
	{
		gen_random(d);
		output_new_gene(d);
		printf("Initial genes are prepared please check [gene_1.dat]\n");
	}
	else if (is_restart == 1)
	{
		int tmp;
		printf("Please input iteration number you want to start");
		scanf("%d", &tmp);
		restart(d, tmp);
	}
}
int main(int argc, char const *argv[])
{
	data data1;
	srand(time(NULL));
	precalculation(data1);
	bool flag=true;
	
	while(flag)
	{
		printf("\n=================\nItereation [%d]\n=================\n", data1.N_iteration);
		printf("Now please prepare your fitness file with name [fitness.dat]\n");

		system("PAUSE");
		update_fitness(data1);
		
		printf("\n-----------------\nPreparing the gene%d\n-----------\n",data1.N_iteration);
		data1.N_iteration++;
		bool again = true;
		while (again)
		{
			gen_mutation(data1);
			output_new_gene(data1);
			printf("\n?? What do you think about the setup from gene%d.dat ??\n[GOOD: 1 ] [BAD: 0]", data1.N_iteration);
			int good=0;
			scanf("%d", &good);
			while (good != 0 && good != 1)
			{
				printf("\nSorry, I don't understand your input, Could you try again:[GOOD: 1 ] [BAD: 0]\n");
				scanf("%d", &good);
			}
			if (good == 1) {
				again = false;
				printf("\nGOOD, let's move onto the next one\n");
			}
			else 
			{
				printf("\nOps!, it's bad, Let me try again\n");
			}

		}
		
		
		
		if (data1.N_iteration == data1.N_iteration_max)
		{
			flag = false;
		}
	}
	return 0;
}
