#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h> 

// Definir o tamanho da matriz aqui
#define N 8

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

 	int i, j, aux;

	MPI_Status status;

// Parâmtros para a topologia cartesiana
	MPI_Comm CART_COMM;
	int dim[2], period[2], reorder;
	int tamDims
	int coods[2]

// Inicia o MPI
	MPI_Init(&argc, &argv);
  	MPI_Comm_size(MPI_COMM_WORLD, &procs);
  	MPI_Comm_rank(MPI_COMM_WORLD, &myId);

// Definir a matriz aqui:
  	if (myId == 0) 
  	{
  		M = {
 // . x ---->
 /* y */ { 00 , 01 , 02 , 03 , 04 , 05 , 06 , 07 },
 /* | */ { 10 , 11 , 12 , 13 , 14 , 15 , 16 , 17 },
 /* | */ { 20 , 21 , 22 , 23 , 24 , 25 , 26 , 27 },
 /* V */ { 30 , 31 , 32 , 33 , 34 , 35 , 36 , 37 },
 		 { 40 , 41 , 42 , 43 , 44 , 45 , 46 , 47 },
 		 { 50 , 51 , 52 , 53 , 54 , 55 , 56 , 57 },
 		 { 60 , 61 , 62 , 63 , 64 , 65 , 66 , 67 },
 		 { 70 , 71 , 72 , 73 , 74 , 75 , 76 , 77 }};
 		
 		for (i = 0; i < N; i++) diag[i] = M[i][i];
  	} else {
  		memset(M, 0, sizeof(int)*N*N);
  	}

// Cria o comunicador cartesiano
  	reorder = 1;
  	tamDims = sqrt(procs);
  	for (i = 0; i < 2; ++i) {
  		period[i] = 0;
  		dim[i] = tamDims;
  	}

  	MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &CART_COMM);

// Captura as coordenadas do processo
  	MPI_Cart_coords(CART_COMM, myId, 2, coords);

// Define o tamanho das sub-matrizes
  	Nx_sub = (coords[0] < N%tamDims) ? N/tamDims + 1 : N/tamDims;   // X
  	Ny_sub = (coords[1] < N%tamDims) ? N/tamDims + 1 : N/tamDims;   // Y

// Define as coordenadas das bordas do bloco (sub-matriz)
  	X_i = (coords[0] < N%tamDims) ? coords[0] * Nx_sub : N%tamDims * (Nx_sub + 1) + (coords[0] - N%tamDims) * Nx_sub;
	X_f = X_i + Nx_sub - 1;

	Y_i = (coords[1] < N%tamDims) ? coords[1] * Ny_sub : N%tamDims * (Ny_sub + 1) + (coords[1] - N%tamDims) * Ny_sub;
  	Y_f = Y_i + Ny_sub - 1;	

// Distribui os blocos da matriz
    if ((sub_M = malloc(Nx_sub * Ny_sub * sizeof(int))) == NULL) { exit(1); }

    MPI_Datatype blocktype;
    MPI_Datatype blocktype2;

    MPI_Type_vector(Ny_sub, Nx_sub, N, MPI_CHAR, &blocktype2);
    MPI_Type_create_resized( blocktype2, 0, sizeof(char), &blocktype);
    MPI_Type_commit(&blocktype);

    int disps[dim[0]*dim[1]];
    int counts[dim[0]*dim[1]];
    for (int ii=0; ii<dim[0]; ii++) {
        for (int jj=0; jj<dim[1]; jj++) {
            disps[ii*dim[1]+jj] = ii*N*Ny_sub+jj*Nx_sub;
            counts [ii*dim[1]+jj] = 1;
        }
    }
 MPI_Scatterv(M, counts, disps, blocktype, sub_M, Ny_sub*Nx_sub, MPI_CHAR, 0, MPI_COMM_WORLD);

// Envia a parcela da diagonal
    if ((sub_diag = malloc(Ny_sub * sizeof(int))) == NULL) { exit(1); }

    for (int proc=0; proc<procs; proc++) {
        if (proc == myId) {
            printf("Rank = %d\n", myId);
            if (rank == 0) {
                printf("Global matrix: \n");
                for (int ii=0; ii<N; ii++) {
                    for (int jj=0; jj<N; jj++) {
                        printf("%3d ",(int)M[ii*N+jj]);
                    }
                    printf("\n");
                }
            }
            printf("Local Matrix:\n");
            for (int ii=0; ii<Ny_sub; ii++) {
                for (int jj=0; jj<Nx_sub; jj++) {
                    printf("%3d ",(int)b[ii*Nx_sub+jj]);
                }
                printf("\n");
            }
            printf("\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }


// Multiplica cada elemento de cada linha pelo elemento da diagonal da linha correspondente
/*  	for (j = 0; j < Ny_sub; ++j) {
  		sub_diag[j];
  		for (i = 0; i < Nx_sub; ++i) {
  			sub_M[i][j] *= aux;
  		}
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
  	}*/


	MPI_Finalize();
 	return 0;
}