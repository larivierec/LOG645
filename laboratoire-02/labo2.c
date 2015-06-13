#include <stdio.h>
#include <sys/time.h>
#include <unistd.h>
#include <stdlib.h>
#include <omp.h>

//declarations de nos methodes sequentielles
void printMatrix(long **mat);
void initMatrix(long **mat, int x);
void problemeNumero1(int nombreIteration, long **matKMoins1,long **mat);
void problemeNumero2(int nombreIteration, long **matKMoins1,long **mat);
void spinWait(int milliseconds);

int main(int argc, char ** argv) {

    //---------------------Initialisation des paramètres---------------------
    int problemeExecuter;
    int valeurInitiale;
    int nombreIteration;

    if (argc != 4) {
        printf("Veuillez entrer le nombre d'arguments comme suivant: #1 Choix du probleme a executer, #2 La valeur initiale des elements de la matrice, #3 le nombre d'iterations\n");
        return 1;
    }

    //assignons les variables passer en paramètres
    problemeExecuter = atoi(argv[1]);
    valeurInitiale = atoi(argv[2]);
    nombreIteration = atoi(argv[3]);
    //---------------------Initialisation des paramètres---------------------

    long **mat;
    mat = (long **)malloc(10*sizeof(long *));
    for(int i=0;i<10;i++){
        mat[i] = (long *)malloc(10*sizeof(long));
    }

    //intialise la matrice avec la valeur passer en paramètre
    initMatrix(mat, valeurInitiale);

    //assigne une nouvelle matrice de type double pointeur vers la matrice mat
    long **matkmoins1 = mat;

    //pour mesurer les millisecondes
    struct timeval startseq, endseq, startpar, endpar;

    //choisir un des problèmes à exécuter
    if(problemeExecuter == 1){

        //Start timer for sequential
        gettimeofday(&startseq,0);

        //Exécution séquentielle
        for(int k = 1; k <= nombreIteration;k++ ){
          for(int i=0; i<10; i++){
            for(int j=0; j<10; j++){
                mat[i][j] = matkmoins1[i][j] + i + j;
                spinWait(50);
            }
          }
          matkmoins1 = mat;
        }

        //End timer for sequential
        gettimeofday(&endseq,0);
        printf("Problème %d : impression de la matrice 'mat' pour l'execution sequentielle :\n", problemeExecuter);
        printMatrix(mat);
        float timeseq = (endseq.tv_sec + (endseq.tv_usec/1000000.0) -  startseq.tv_sec + (startseq.tv_usec/1000000.0));
        printf("Temps que l'application a pris pour terminer (en secondes): %f \n", timeseq);

        //-------------------------------------------------------
        initMatrix(mat, valeurInitiale);
        initMatrix(matkmoins1, valeurInitiale);

        //Start timer for parallel
        gettimeofday(&startpar,0);

        //omp_set_num_threads(10);
        #pragma omp parallel default(shared)
        {
            for(int k=1;k<=nombreIteration;k++){
                #pragma omp for
                for(int i=0;i<10;i++){
                    for(int j=0;j<10;j++){
                        mat[i][j] = matkmoins1[i][j] + (i + j);
                        matkmoins1[i][j] = mat[i][j];
                        spinWait(50);
                    }
                }
            }
        }

        //End timer for parallel
        gettimeofday(&endpar,0);
        printf("Problème %d : impression de la matrice 'mat' pour l'execution parallele :\n", problemeExecuter);
        printMatrix(mat);
        float timepar = (endpar.tv_sec + (endpar.tv_usec/1000000.0) -  startpar.tv_sec + (startpar.tv_usec/1000000.0));
        printf("Temps que l'application a pris pour terminer (en secondes): %f \n", timepar);
        printf("L'acceleration de l'application parallele vis-a-vis l'application sequentielle est %f \n", (timeseq/timepar));

    }else if(problemeExecuter == 2){

        //Start timer for sequential
        gettimeofday(&startseq,0);

        //Exécution séquentielle
        for(int k = 1; k <= nombreIteration;k++ ){
          for(int i=0; i<10; i++){
            for(int j=9; j>=0; j--){
                if(j==9){
                    mat[i][j] = matkmoins1[i][j] + i;
                }else{
                    mat[i][j] = matkmoins1[i][j] + mat[i][j+1];
                }
                spinWait(50);
            }
          }
          matkmoins1 = mat;
        }

        //End timer for sequential
        gettimeofday(&endseq,0);
        printf("Problème %d : impression de la matrice 'mat' pour l'execution sequentielle :\n", problemeExecuter);
        printMatrix(mat);
        float timeseq = (endseq.tv_sec + (endseq.tv_usec/1000000.0) -  startseq.tv_sec + (startseq.tv_usec/1000000.0));
        printf("Temps que l'application a pris pour terminer (en secondes): %f \n", timeseq);

        //-------------------------------------------------------
        initMatrix(mat, valeurInitiale);
        initMatrix(matkmoins1, valeurInitiale);

        //Start timer for parallel
        gettimeofday(&startpar,0);

        //omp_set_num_threads(10);
        #pragma omp parallel default(shared)
        {
            for(int k=1;k<=nombreIteration;k++){
                #pragma omp for
                for(int i=0;i<10;i++){
                    for(int j=9; j>=0; j--){
                        if(j==9){
                            mat[i][j] = matkmoins1[i][j] + i;
                        }else{
                            mat[i][j] = matkmoins1[i][j] + mat[i][j + 1];
                        }
                        matkmoins1[i][j] = mat[i][j];
                        spinWait(50);
                    }
                }
            }
        }

        //End timer for parallel
        gettimeofday(&endpar,0);
        printf("Problème %d : impression de la matrice 'mat' pour l'execution parallele :\n", problemeExecuter);
        printMatrix(mat);
        float timepar = (endpar.tv_sec + (endpar.tv_usec/1000000.0) -  startpar.tv_sec + (startpar.tv_usec/1000000.0));
        printf("Temps que l'application a pris pour terminer (en secondes): %f \n", timepar);
        printf("L'acceleration de l'application parallele vis-a-vis l'application sequentielle est %f \n", (timeseq/timepar));

    }else{
        return 2;
    }
    return 0;
}

//la méthode pour imprimer une matrice à l'écran
void printMatrix(long **mat){
    for(int i=0; i<10; i++){
        for(int j=0; j<10; j++){
            printf("%ld ", mat[i][j]);
        }
        printf("\n");
    }
}

//la méthode pour initialiser la matrice à valeur: x
void initMatrix(long **mat, int x){
    for(int i=0; i<10; i++){
        for(int j=0; j<10; j++){
            mat[i][j] = x;
        }
    }
}//fin initMatrix();

void spinWait(int milliseconds)
{
    struct timeval startTime;
    struct timeval endTime;
    gettimeofday(&startTime, NULL);
    do
    {
        gettimeofday(&endTime, NULL);
    } while ((endTime.tv_sec - startTime.tv_sec) * 1000000 + (endTime.tv_usec - startTime.tv_usec) < milliseconds * 1000);
}
