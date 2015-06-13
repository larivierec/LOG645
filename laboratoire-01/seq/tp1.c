#include <stdio.h>
#include <sys/time.h>
#include <unistd.h>
#include <stdlib.h>


//declarations de nos methodes sequentielles

void printMatrix(long **mat);
void initMatrix(long **mat, int x);
void problemeNumero1(int nombreIteration, long **matKMoins1,long **mat);
void problemeNumero2(int nombreIteration, long **matKMoins1,long **mat);

//main

int main(int argc, char *argv[])
{

    //pour mesurer les millisecondes
    struct timeval start, end;

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


    gettimeofday(&start,0);

    //matrix initiale
    long **mat;

    //initialise le mémoire de la matrice

    mat = (long **)malloc(8*sizeof(long *));
    for(int i=0;i<8;i++){
        mat[i] = (long *)malloc(8*sizeof(long));
    }

    //intialise la matrice avec la valeur passer en paramètre
    initMatrix(mat, valeurInitiale);

    //assigne une nouvelle matrice de type double pointeur vers la matrice mat
    long **matkmoins1 = mat;

    //choisir un des problèmes à exécuter
    if(problemeExecuter == 1){
      problemeNumero1(nombreIteration,matkmoins1, mat);
    }else if(problemeExecuter == 2){
      problemeNumero2(nombreIteration,matkmoins1, mat);
    }else{
        return 2;
    }
    gettimeofday(&end,0);


    printf("Problème %d : impression de la matrice 'mat' :\n", problemeExecuter);
    printMatrix(mat);
    printf("Temps que l'application a pris pour terminer (en secondes): %f \n",(end.tv_sec + (end.tv_usec/1000000.0) -  start.tv_sec + (start.tv_usec/1000000.0)));

    return 0;
}//fin main


//la méthode pour imprimer une matrice à l'écran
void printMatrix(long **mat){
    for(int i=0; i<sizeof(mat); i++){
        for(int j=0; j<sizeof(mat[0]); j++){
            printf("%ld ", mat[i][j]);
        }
        printf("\n");
    }
}


//la méthode pour initialiser la matrice à valeur: x
void initMatrix(long **mat, int x){
    for(int i=0; i<sizeof(mat); i++){
        for(int j=0; j<sizeof(mat[0]); j++){
            mat[i][j] = x;
        }
    }
}//fin initMatrix();


//l'algorithme qui calcule la matrice résultante du problème un
void problemeNumero1(int nombreIteration, long **matKMoins1,long **mat){
  for(int k = 1; k <= nombreIteration;k++ ){
    for(int i=0; i<sizeof(mat); i++){
      for(int j=0; j<sizeof(mat[0]); j++){
          usleep(1000);
          mat[i][j] = matKMoins1[i][j] + (i+j) * k;
      }
    }
    matKMoins1 = mat;
  }
}//fin problemeNumero1

//l'algorithme qui calcule la matrice résultante du problème du
void problemeNumero2(int nombreIteration, long **matKMoins1,long **mat){
  for(int k = 1; k <= nombreIteration;k++){
    for(int i=0; i<sizeof(mat); i++){
      for(int j=0; j<sizeof(mat[0]); j++){
        usleep(1000);
        if(j == 0){
          mat[i][j] = matKMoins1[i][j] + (i*k);
        }else{
          mat[i][j] = matKMoins1[i][j] + mat[i][j-1] * k;
        }
      }
    }
    matKMoins1 = mat;
  }
}//fin problemeNumero1
