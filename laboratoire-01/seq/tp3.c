#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#define TEMPS_ATTENTE 5

int main(int argc, char *argv[])
{
    //variables pour le programme
    int nombreLignes, nombreColonnes, nombrePasTemps, tempsDiscretise, tailleSubdivision, nbProc;

    int idProc;
    struct timeval start, end;


    if (argc != 7) {
        printf("Veuillez entrer le nombre d'arguments comme suivant: #1 Nombre de lignes de la matrice, #2 Nombre de colonnes de la matrice, #3 le nombre de pas de temps, #4 le temps discrétisé, #5 la taille d'un côté d'une subdivision, #6 le nombre de processus à utiliser.\n");
        return 1;
    }

    nombreLignes = atoi(argv[1]);
    nombreColonnes = atoi(argv[2]);
    nombrePasTemps = atoi(argv[3]);
    tempsDiscretise = atoi(argv[4]);
    tailleSubdivision = atoi(argv[5]);
    nbProc = atoi(argv[6]);


        //Sequential

        float **mat, **matkmoins1;

        //initialise le mémoire de la matrice
        mat = (float **)malloc(nombreLignes*sizeof(float *));
        for(int i=0;i<nombreLignes;i++){
            mat[i] = (float *)malloc(nombreColonnes*sizeof(float));
        }

        //Initialisation de la matrice
        for(int y=0; y<nombreLignes; y++){
            for(int x=0; x<nombreColonnes; x++){
                usleep(TEMPS_ATTENTE);
                mat[y][x] = (x*(nombreColonnes - x - 1)) * (y*(nombreLignes - y - 1));
            }
        }
        matkmoins1 = mat;

        //Imprimer la plaque
        printf("Impression de la plaque initiale (mat) :\n");
        for(int i=nombreLignes -1 ; i>=0; i--){
            for(int j=0; j<nombreColonnes; j++){
                printf("%.1f  ", mat[i][j]);
            }
            printf("\n");
        }

        float part1formula, part2formula, part3formula;


        for(int k=1; k<=nombrePasTemps; k++){
            for(int y=0; y<nombreLignes; y++){
                for(int x=0; x<nombreColonnes; x++){
                    usleep(TEMPS_ATTENTE);

                    part1formula = (1-(4*(tempsDiscretise/(pow(tailleSubdivision,2))))) * matkmoins1[y][x];
                    part2formula = (tempsDiscretise/(pow(tailleSubdivision,2)));
                    part3formula = 0;
                    printf("x : %d, y : %d\n", x, y);
                    if(y != 0){
                        printf("CRASH1.1 ?\n");
                        part3formula = part3formula + matkmoins1[y-1][x];
                    }
                    printf("CRASH1 ?\n");
                    if(y != nombreLignes - 1){
                        part3formula += matkmoins1[y+1][x];
                    }
                    printf("CRASH2 ?\n");
                    if(x != 0){
                        part3formula += matkmoins1[y][x-1];
                    }
                    printf("CRASH3 ?\n");
                    if(x != nombreColonnes - 1){
                        part3formula += matkmoins1[y][x+1];
                    }
                    printf("CRASH4 ?\n");
                    mat[y][x] = part1formula + (part2formula * part3formula);
                }
            }

            matkmoins1 = mat;

        }



        //We start the timer
        //gettimeofday(&start,0);




/*
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
*/


    return 0;
}
