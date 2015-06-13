#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include "mpi.h"
#include <sys/time.h>


int main(int argc, char *argv[])
{
    //variables pour le programme
    int problemeExecuter, valeurInitiale, nombreIteration;

    int idProc;
    struct timeval start, end;


    //intialise MPI avec son MPI_COMM_WORLD
    int err = MPI_Init(&argc, &argv);

    if(err != MPI_SUCCESS){
      printf("Il y a eu une erreur lors de l'initialisation de MPI");
      return(3);
    }

    if (argc != 4) {
        printf("Veuillez entrer le nombre d'arguments comme suivant: #1 Choix du probleme a executer, #2 La valeur initiale des elements de la matrice, #3 le nombre d'iterations\n");
        MPI_Finalize();
        return 1;
    }
    problemeExecuter = atoi(argv[1]);
    valeurInitiale = atoi(argv[2]);
    nombreIteration = atoi(argv[3]);

    //Cette méthode est utilisée pour récupérer le numéro de processeur que nous

    gettimeofday(&start,0);
    MPI_Comm_rank(MPI_COMM_WORLD, &idProc);

    long *arrTemp;
    long *arrTempKMoins1;
    long arrGlobal[64];

    //Valeur qui sera changé à arrTemp[3] pour être envoyer à la matrice comme valeur de k-1 probleme 2
    long valeurColonne3;

    //initialise les matrice de pointeurs
    arrTemp = (long *)malloc(4 * sizeof(long) );
    arrTempKMoins1 = (long *)malloc(4 * sizeof(long) );

    //initialise chaque matrice de chaque processeur à la valeur passée en paramètre
    for(int i=0;i<4;i++){
        arrTemp[i] = valeurInitiale;
        arrTempKMoins1[i] = valeurInitiale;
    }


    //identifie si le numéro du processeur s'agit d'un numéro de processeur
    //pair / impair et initialise sa ligne avec la formule suivante
    //initialise la ligne et ses colonnes

    int matrixRow;
    int matrixStartCol;
    int matrixEndCol;

    if(idProc % 2 != 0){
        matrixRow = (idProc-1)/2;
        matrixStartCol = 4;
        matrixEndCol = 7;
    }else{
        matrixRow = (idProc)/2;
        matrixStartCol = 0;
        matrixEndCol = 3;
    }

    //pour chaque itération k effectue un calcul

    for(int k=1;k<=nombreIteration;k++){
        //Problème 2
        //si le processeur provient d'un nombre impair nous devons attendre les données du processeur pair
        if(idProc % 2 != 0 && problemeExecuter == 2){
            MPI_Recv(&valeurColonne3, 1, MPI_LONG, idProc-1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        for(int i=0;i<4;i++){
            //sleep 1 milliseconde par calcul
            usleep(1000);
            if(problemeExecuter == 1){
                arrTemp[i] = arrTempKMoins1[i] + (matrixRow + (matrixStartCol + i)) * k;
            }else{
                if(matrixStartCol == 0 && i == 0){
                    arrTemp[i] = arrTempKMoins1[i] + matrixRow * k;
                }else{
                    if(matrixStartCol == 4 && i == 0){
                        arrTemp[i] = arrTempKMoins1[i] + valeurColonne3 * (long)(k);
                    }else{
                        arrTemp[i] = arrTempKMoins1[i] + arrTemp[i-1] * k;
                    }
                }
            }
        }
        //Problème 2
        //si le processeur provient d'un numéro de processeur pair envoyer les
        //aux processeurs impairs pour qu'ils puissent effectué leur travail
        if(idProc % 2 == 0 && problemeExecuter == 2){
            MPI_Send(&arrTemp[3], 1, MPI_LONG, idProc+1, idProc, MPI_COMM_WORLD);
        }
        arrTempKMoins1 = arrTemp;
    }

    //Processeur 0 ici va récupérer tout les données de chaque processeur (incluant elle-même)
    //et les stocker dans la matrice global
    MPI_Gather(arrTemp, 4, MPI_LONG, arrGlobal, 4, MPI_LONG, 0, MPI_COMM_WORLD);


    //impression de la matrice finale
    if(idProc == 0){
        gettimeofday(&end,0);
        printf("Problème %d : impression de la matrice 'mat' :\n", problemeExecuter);
        for(int i=0;i<64;i++){
            if(i % 8 == 0 && i != 0){
                printf("\n");
            }
            printf("%ld " , arrGlobal[i]);
        }
        printf("\n");
        printf("Temps que l'application a pris pour terminer (en secondes): %f \n",(end.tv_sec + (end.tv_usec/1000000.0) -  start.tv_sec + (start.tv_usec/1000000.0)));
    }


    MPI_Finalize();
    return 0;
}
