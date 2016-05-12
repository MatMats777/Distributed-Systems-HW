#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>

// Definir o tamanho da matriz aqui
#define N 8

int calculaFator(int n){
	int max_fact = 	1;

	while (n/max_fact > max_fact) { 
		max_fact++; 
		while(n%max_fact!=0) max_fact++;
	}

	return max_fact;
}

int main(int argc, char *argv[]) {
// Variáveis do MPI
	int procs;
	int myId;

// Outras variáveis
	int Nx_sub, Ny_sub;
	int X_i, X_f;
	int Y_i, Y_f;

	int M[N][N];
	int **sub_M;
	int diag[N];
	int *sub_diag;
	int soma_p[N], soma_recv[N], soma;

	int i, j;

	MPI_Status status;

// Parâmtros para a topologia cartesiana
	MPI_Comm CART_COMM;
	int dim[2], period[2], reorder;
	int coords[2];

// Inicia o MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myId);

// Definir a matriz aqui:
	if (myId == 0){
// . x ---->
/* y */ M[0][0] = 00; M[0][1] = 01; M[0][2] = 02; M[0][3] = 03; M[0][4] = 04; M[0][5] = 05; M[0][6] = 06; M[0][7] = 07;
/* | */ M[1][0] = 10; M[1][1] = 11; M[1][2] = 12; M[1][3] = 13; M[1][4] = 14; M[1][5] = 15; M[1][6] = 16; M[1][7] = 17;
/* | */ M[2][0] = 20; M[2][1] = 21; M[2][2] = 22; M[2][3] = 23; M[2][4] = 24; M[2][5] = 25; M[2][6] = 26; M[2][7] = 27;
/* V */ M[3][0] = 30; M[3][1] = 31; M[3][2] = 32; M[3][3] = 33; M[3][4] = 34; M[3][5] = 35; M[3][6] = 36; M[3][7] = 37;
		M[4][0] = 40; M[4][1] = 41; M[4][2] = 42; M[4][3] = 43; M[4][4] = 44; M[4][5] = 45; M[4][6] = 46; M[4][7] = 47;
		M[5][0] = 50; M[5][1] = 51; M[5][2] = 52; M[5][3] = 53; M[5][4] = 54; M[5][5] = 55; M[5][6] = 56; M[5][7] = 57;
		M[6][0] = 60; M[6][1] = 61; M[6][2] = 62; M[6][3] = 63; M[6][4] = 64; M[6][5] = 65; M[6][6] = 66; M[6][7] = 67;
		M[7][0] = 70; M[7][1] = 71; M[7][2] = 72; M[7][3] = 73; M[7][4] = 74; M[7][5] = 75; M[7][6] = 76; M[7][7] = 77;

		for (i = 0; i < N; i++) diag[i] = M[i][i];
	} else {
		memset(M, 0, sizeof(int)*N*N);
	}

// Cria o comunicador cartesiano
	reorder = 0;
	dim[1] = calculaFator(procs);
	dim[0] = procs/dim[1];
	
	period[0] = period[1] = 0;

	MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &CART_COMM);

// Captura as coordenadas do processo
	MPI_Cart_coords(CART_COMM, myId, 2, coords);

// Define o tamanho das sub-matrizes
	Nx_sub = N/dim[0];   // X
	Ny_sub = N/dim[1];   // Y

// Define as coordenadas das bordas do bloco (sub-matriz)
	X_i = coords[0] * Nx_sub;
	X_f = X_i + Nx_sub - 1;

	Y_i = coords[1] * Ny_sub;
	Y_f = Y_i + Ny_sub - 1;

// Distribui os blocos da matriz
	int *p;
	p = (int *)malloc(Nx_sub*Ny_sub*sizeof(int)); 

	sub_M = (int **)malloc(Nx_sub * sizeof(int*)); 
	for ( i = 0; i < Ny_sub; ++i) {
		sub_M[i] = &(p[i*Nx_sub]);
	}

// Definição de novos tipos para o envio da matriz
	int sizes[2]    = {N, N};         
	int subsizes[2] = {Nx_sub, Ny_sub};     
	int starts[2]   = {X_i,Y_i};                        
	MPI_Datatype type, subarrtype;
	MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_INT, &type);
	MPI_Type_create_resized(type, 0, Ny_sub*sizeof(int), &subarrtype);
	MPI_Type_commit(&subarrtype);

	int *globalptr=NULL;
	if (myId == 0)
		globalptr = &(M[0][0]);

	int sendcounts[procs];
	int displs[procs];

	if (myId == 0) {
		for (i=0; i<procs; i++)
			sendcounts[i] = 1;

		int disp = 0;
		
		for (i=0; i<dim[1]; i++) {
			for (j=0; j<dim[0]; j++) {
				displs[i*dim[0]+j] = disp;
				disp += 1;
			}
			disp += ((Ny_sub)-1)*dim[0];
		}
	}

// Envio da matriz
	MPI_Scatterv(globalptr, sendcounts, displs, subarrtype, &(sub_M[0][0]),
		Nx_sub*Ny_sub, MPI_INT,
		0, MPI_COMM_WORLD);

// Imprime as submatrizes e matriz principal    
	for (i = 0; i < procs; ++i) {
		if (myId == i){
			printf("Rank = %d\n", myId);
			if (myId == 0) {
				printf("dim[0] = %d\ndim[1] = %d\n", dim[0], dim[1]);
				printf("Global matrix: \n");
				for (int ii=0; ii<N; ii++) {
					for (int jj=0; jj<N; jj++) {
						printf("%3d ", M[ii][jj]);
					}
					printf("\n");
				}
			}
			printf("Local Matrix:\n");
			for (int ii=0; ii<Ny_sub; ii++) {
				for (int jj=0; jj<Nx_sub; jj++) {
					printf("%3d ", sub_M[ii][jj]);
				}
				printf("\n");
			}
			printf("\n");

		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

// Envia a diagonal
	MPI_Bcast(diag, N, MPI_INT, 0, MPI_COMM_WORLD);

// Multiplica cada elemento de cada linha pelo elemento da diagonal da linha correspondente
	for (i = 0; i < Nx_sub; ++i) {
		for (j = 0; j < Ny_sub; ++j) {
			sub_M[i][j] *= diag[i+X_i];
		}
	}

// Imprime as submatrizes e matriz principal novamente
	for (i = 0; i < procs; ++i) {
		if (myId == i){
			printf("Rank = %d\n", myId);
			if (myId == 0) {
				printf("Global matrix: \n");
				for (int ii=0; ii<N; ii++) {
					for (int jj=0; jj<N; jj++) {
						printf("%3d ", M[ii][jj]);
					}
					printf("\n");
				}
			}
			printf("Local Matrix:\n");
			for (int ii=0; ii<Ny_sub; ii++) {
				for (int jj=0; jj<Nx_sub; jj++) {
					printf("%3d ", sub_M[ii][jj]);
				}
				printf("\n");
			}
			printf("\n");

		}
		MPI_Barrier(MPI_COMM_WORLD);
	}


// Soma das colunas em cada processo
	memset(soma_p, 0, sizeof(int)*N);
	for (i = 0; i < Nx_sub; ++i) {
		for (j = 0; j < Ny_sub; ++j) {
			soma_p[i + X_i] += sub_M[i][j];
		}
	}
	
// Junta as somas parciais
	soma = 0;
	MPI_Reduce(soma_p, soma_recv, N, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	if (myId == 0) {
		for (i = 0; i < N; ++i) {
			soma += soma_recv[i];
		}

		printf("Resultado: %d\n", soma);
	}


	MPI_Finalize();
	return 0;
	}
