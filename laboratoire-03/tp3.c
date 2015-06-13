#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include "mpi.h"
#include <sys/time.h>
#include <math.h>
#define TEMPS_ATTENTE 5

int main(int argc, char *argv[])
{
    //variables pour le programme
    int nombreLignes, nombreColonnes, nombrePasTemps, nbProc;
    float tempsDiscretise, tailleSubdivision;

    int idProc;
    struct timeval start, end;

    float vy, vx, hy, hx, v2y, v2x, h2y, h2x;
    int vnbProcUse, hnbProcUse, v2nbProcUse, h2nbProcUse, nbCommunications, vaire, haire, v2aire, h2aire;
    int ysizestart, ysizeend, xsizestart, xsizeend;
    MPI_Request req, req2;
    MPI_Status stat, stat2;

    if (argc != 7) {
        printf("Veuillez entrer le nombre d'arguments comme suivant: #1 Nombre de lignes de la matrice, #2 Nombre de colonnes de la matrice, #3 le nombre de pas de temps, #4 le temps discrétisé, #5 la taille d'un côté d'une subdivision, #6 le nombre de processus à utiliser.\n");
        return 1;
    }

    nombreLignes = atoi(argv[1]);\
    nombreColonnes = atoi(argv[2]);
    nombrePasTemps = atoi(argv[3]);
    tempsDiscretise = atof(argv[4]);
    tailleSubdivision = atof(argv[5]);
    nbProc = atoi(argv[6]);

    //intialise MPI avec son MPI_COMM_WORLD
    int err = MPI_Init(&argc, &argv);

    if(err != MPI_SUCCESS){
      printf("Il y a eu une erreur lors de l'initialisation de MPI");
      return(3);
    }

    //Cette méthode est utilisée pour récupérer le numéro de processeur que nous
    MPI_Comm_rank(MPI_COMM_WORLD, &idProc);
    MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

    float **matseq;
    matseq = (float **)malloc(nombreLignes*sizeof(float *));
    for(int i=0;i<nombreLignes;i++){
        matseq[i] = (float *)malloc(nombreColonnes*sizeof(float));
    }
    float **matseqkplus1;
    matseqkplus1 = (float **)malloc(nombreLignes*sizeof(float *));
    for(int i=0;i<nombreLignes;i++){
        matseqkplus1[i] = (float *)malloc(nombreColonnes*sizeof(float));
    }

    for(int y=0; y<nombreLignes; y++){
        for(int x=0; x<nombreColonnes; x++){
            usleep(TEMPS_ATTENTE);
            matseq[y][x] = (x*(nombreColonnes - x - 1)) * (y*(nombreLignes - y - 1));
        }
    }

    //intialisation du timer apres l'initialisation

    double timeStartSeq, timeEndSeq, TexecSeq;
    struct timeval ts;
    gettimeofday (&ts, NULL); // Début du chronomètre
    timeStartSeq = (double) (ts.tv_sec) + (double) (ts.tv_usec) / 1e6;

    if(idProc == 0){
        //Imprimer la plaque
        printf("Impression de la plaque initiale (matseq) :\n");
        for(int i=nombreLignes -1 ; i>=0; i--){
            for(int j=0; j<nombreColonnes; j++){
                printf("%.1f  ", matseq[i][j]);
            }
            printf("\n");
        }

        //CODE SEQUENTIEL
        float part1formula, part2formula, part3formula;

        for(int k=1; k<=nombrePasTemps; k++){
            for(int y=0; y<nombreLignes; y++){
                for(int x=0; x<nombreColonnes; x++){
                    usleep(TEMPS_ATTENTE);
                    if(x==0 || x==(nombreColonnes-1) || y==0 || y==(nombreLignes-1)){
                       matseqkplus1[y][x] = 0;
                    }else{

                        part1formula = (1-4*tempsDiscretise/tailleSubdivision * tailleSubdivision) * matseq[y][x];
                        part2formula = (tempsDiscretise/tailleSubdivision * tailleSubdivision);
                        part3formula = 0;

                        if(y != 0){
                            part3formula += matseq[y-1][x];
                        }
                        if(y != nombreLignes - 1){
                            part3formula += matseq[y+1][x];
                        }
                        if(x != 0){
                            part3formula += matseq[y][x-1];
                        }
                        if(x != nombreColonnes - 1){
                            part3formula += matseq[y][x+1];
                        }
                        matseqkplus1[y][x] = part1formula + (part2formula * part3formula);
                    }
                }
            }

            //copier la plaque

            for(int y=0; y<nombreLignes; y++){
                for(int x=0; x<nombreColonnes; x++){
                    matseq[y][x] = matseqkplus1[y][x];
                }
            }
        }

        //Imprimer la plaque apres n iterations
        printf("Impression de la plaque apres %d iteration (matseqkplus1) :\n", nombrePasTemps);
        for(int i=nombreLignes -1 ; i>=0; i--){
            for(int j=0; j<nombreColonnes; j++){
                printf("%.5f  ", matseqkplus1[i][j]);
            }
            printf("\n");
        }

        gettimeofday (&ts, NULL); // Fin du chronomètre
        timeEndSeq = (double) (ts.tv_sec) + (double) (ts.tv_usec) / 1e6;
        TexecSeq = timeEndSeq - timeStartSeq; //Temps d'exécution en secondes
        // Fin de l’exemple
        if(idProc == 0){
            printf("Temps d'execution Séquentiel: %f\n", TexecSeq);
        }
    }

    //CODE PARALLEL
    float **matpar;
    //initialise le mémoire de la matrice
    matpar = (float **)malloc(nombreLignes*sizeof(float *));
    for(int i=0;i<nombreLignes;i++){
        matpar[i] = (float *)malloc(nombreColonnes*sizeof(float));
    }
    //Initialisation de la matrice
    for(int y=0; y<nombreLignes; y++){
        for(int x=0; x<nombreColonnes; x++){
            usleep(TEMPS_ATTENTE);
            matpar[y][x] = (x*(nombreColonnes - x - 1)) * (y*(nombreLignes - y - 1));
        }
    }

    //Evaluation du matrice

    if(nombreColonnes >= nbProc){
        //Si on fait des groupes horizontales
        vy = (float)(nombreColonnes) / (float)(nbProc);
        vx = (float)(nombreLignes);
        vnbProcUse = nbProc;
    }else{
        vy = 1.0f;
        vx = (float)nombreLignes;
        vnbProcUse = nombreColonnes;
    }
    nbCommunications = vnbProcUse - 1;
    vaire = ((int)(vy) - 4) * (int)(vx);

    if(nombreLignes >= nbProc){
        //Si on fait des groupes verticales
        hy = (float)(nombreColonnes);
        hx = (float)(nombreLignes) / (float)(nbProc);
        hnbProcUse = nbProc;
    }else{
        hy = (float)nombreColonnes;
        hx = 1.0f;
        hnbProcUse = nombreLignes;
    }
    nbCommunications = hnbProcUse - 1;
    haire = (int)(hy) * ((int)(hx) - 4);

    //We start the timer
    double timeStart, timeEnd, Texec;
    struct timeval tp;
    gettimeofday (&tp, NULL); // Début du chronomètre
    timeStart = (double) (tp.tv_sec) + (double) (tp.tv_usec) / 1e6;

    //Coupure verticale
    if(nombreColonnes >= nombreLignes){ //Vertical cut
        if(idProc < vnbProcUse){

            ysizestart = (int)((int)vy * (int)(idProc)) + floorf(((float)(idProc) * (float)(vy - (int)vy)) * 10000000) / 10000000;
            ysizeend = ((int)((int)vy + (int)((idProc) * (int)vy)) + floorf(((float)(idProc + 1) * (float)(vy - (int)vy)) * 10000000) / 10000000) - 1;
            xsizestart = 0;
            xsizeend = (int)vx - 1;


            //--------------------Initialize variable for calculation---------------------------------------------
            float **mattocalculate;
            mattocalculate = (float **)malloc((int)((xsizeend - xsizestart) + 1)*sizeof(float *));
            for(int i=0;i<(int)((xsizeend - xsizestart) + 1);i++){
                mattocalculate[i] = (float *)malloc((int)((ysizeend - ysizestart) + 1)*sizeof(float));
            }
            float **mattocalculatekplus1;
            mattocalculatekplus1 = (float **)malloc((int)((xsizeend - xsizestart) + 1)*sizeof(float *));
            for(int i=0;i<(int)((xsizeend - xsizestart) + 1);i++){
                mattocalculatekplus1[i] = (float *)malloc((int)((ysizeend - ysizestart) + 1)*sizeof(float));
            }

            for(int y=0; y<=(ysizeend - ysizestart); y++){
                for(int x=(xsizeend - xsizestart) ; x>=0; x--){
                    mattocalculatekplus1[x][y] = 0;
                }
            }

            //Left column to send
            float *matleftcolumn;
            matleftcolumn = (float *)malloc((int)((xsizeend - xsizestart) + 1)*sizeof(float *));
            //Right column to send
            float *matrightcolumn;
            matrightcolumn = (float *)malloc((int)((xsizeend - xsizestart) + 1)*sizeof(float *));
            //Left column to send
            float *matleftcolumnOfProc;
            matleftcolumnOfProc = (float *)malloc((int)((xsizeend - xsizestart) + 1)*sizeof(float *));
            //Right column to send
            float *matrightcolumnOfProc;
            matrightcolumnOfProc = (float *)malloc((int)((xsizeend - xsizestart) + 1)*sizeof(float *));
            //--------------------Initialize variable for calculation---------------------------------------------


            int mat_y = 0;
            int mat_x = 0;
            float part1formula, part2formula, part3formula;
            int nbCaseCalculated = 0;
            int totalCaseToCalculate = ((ysizeend - ysizestart) + 1) * ((xsizeend - xsizestart) + 1);
            int nbColumns = (ysizeend - ysizestart);
            for(int k=1; k<=nombrePasTemps; k++){

                mat_y = 0;
                mat_x = 0;
                MPI_Request leftrequest, rightrequest;
                MPI_Status leftstatus, rightstatus;
                int leftflag = 0, rightflag = 0;
                int betweenindex = 1;
                nbCaseCalculated = 0;
                int initrecvleft = 0, initrecvright = 0;
                while(nbCaseCalculated < totalCaseToCalculate){

                    //--------------------First iteration-----------------------------------------------------------------
                    if(k == 1){ //If first iteration
                        //----------Left column START----------
                        mat_y = 0;
                        mat_x = 0; //Column equals zero
                        for(int y=xsizestart; y<=xsizeend; y++){
                            usleep(TEMPS_ATTENTE);
                            if( idProc == 0 || (idProc == (vnbProcUse - 1) && nbColumns == 0) || y==0 || y==(nombreLignes-1)){
                               mattocalculate[mat_y][mat_x] = 0.0;
                               matleftcolumn[mat_y] = 0.0;
                            }else{
                                part1formula = (1-4*tempsDiscretise/tailleSubdivision * tailleSubdivision) * matpar[y][ysizestart];
                                part2formula = (tempsDiscretise/tailleSubdivision * tailleSubdivision);
                                part3formula = 0;
                                if(y != 0){part3formula += matpar[y-1][ysizestart];}
                                if(y != nombreLignes - 1){part3formula += matpar[y+1][ysizestart];}
                                if(ysizestart != 0){part3formula += matpar[y][ysizestart-1];}
                                if(ysizestart != nombreColonnes - 1){part3formula += matpar[y][ysizestart+1];}
                                float tmpresult = part1formula + (part2formula * part3formula);
                                mattocalculate[mat_y][mat_x] = tmpresult;
                                matleftcolumn[mat_y] = tmpresult;
                            }
                            mat_y++;
                            nbCaseCalculated++;
                        }
                        if(idProc != 0 && nombrePasTemps > 1){
                            MPI_Isend(matleftcolumn, nombreLignes, MPI_FLOAT, idProc-1, k, MPI_COMM_WORLD, &leftrequest);
                        }
                        //----------Left column END----------

                        //----------Right column START----------
                        if(nbColumns > 0){ //0 equals 1 column, 1 equals 2 columns, etc...
                            mat_y = 0;
                            mat_x = nbColumns; //Column equals zero
                            for(int y=xsizestart; y<=xsizeend; y++){
                                usleep(TEMPS_ATTENTE);
                                if(idProc == (vnbProcUse - 1) || y==0 || y==(nombreLignes-1)){
                                   mattocalculate[mat_y][mat_x] = 0.0;
                                   matrightcolumn[mat_y] = 0.0;
                                }else{
                                    part1formula = (1-4*tempsDiscretise/tailleSubdivision * tailleSubdivision) * matpar[y][ysizeend];
                                    part2formula = (tempsDiscretise/tailleSubdivision * tailleSubdivision);
                                    part3formula = 0;
                                    if(y != 0){part3formula += matpar[y-1][ysizeend];}
                                    if(y != nombreLignes - 1){part3formula += matpar[y+1][ysizeend];}
                                    if(ysizeend != 0){part3formula += matpar[y][ysizeend-1];}
                                    if(ysizeend != nombreColonnes - 1){part3formula += matpar[y][ysizeend+1];}
                                    float tmpresult = part1formula + (part2formula * part3formula);
                                    mattocalculate[mat_y][mat_x] = tmpresult;
                                    matrightcolumn[mat_y] = tmpresult;
                                }
                                mat_y++;
                                nbCaseCalculated++;
                            }
                            if(idProc != (vnbProcUse-1)  && nombrePasTemps > 1){
                                MPI_Isend(matrightcolumn, nombreLignes, MPI_FLOAT, idProc+1, k, MPI_COMM_WORLD, &rightrequest);
                            }
                        }else{
                            if(idProc != (vnbProcUse-1)  && nombrePasTemps > 1){
                                MPI_Isend(matleftcolumn, nombreLignes, MPI_FLOAT, idProc+1, k, MPI_COMM_WORLD, &rightrequest);
                            }
                        }
                        //----------Right column END----------

                        //----------Between columns START----------
                        if(nbColumns > 1){
                            mat_y = 0;
                            mat_x = 1; //Second column
                            for(int y=xsizestart; y<=xsizeend; y++){
                                mat_x = 1;
                                for(int x=ysizestart+1; x<=ysizeend-1; x++){
                                    usleep(TEMPS_ATTENTE);
                                    if(y == 0 || y == (nombreLignes - 1)){
                                        mattocalculate[mat_y][mat_x] = 0.0;
                                    }else{
                                        part1formula = (1-4*tempsDiscretise/tailleSubdivision * tailleSubdivision) * matpar[y][x];
                                        part2formula = (tempsDiscretise/tailleSubdivision * tailleSubdivision);
                                        part3formula = 0;
                                        if(y != 0){part3formula += matpar[y-1][x];}
                                        if(y != nombreLignes - 1){part3formula += matpar[y+1][x];}
                                        if(x != 0){part3formula += matpar[y][x-1];}
                                        if(x != nombreColonnes - 1){part3formula += matpar[y][x+1];}
                                        mattocalculate[mat_y][mat_x] = part1formula + (part2formula * part3formula);
                                    }
                                    mat_x++;
                                    nbCaseCalculated++;
                                }
                                mat_y++;
                            }
                        }
                        //----------Between columns END----------

                    //--------------------Other iterations----------------------------------------------------------------
                  }else{ //Tout autre iteration

                        //Initialise le Ireceive
                        if(initrecvleft == 0 && idProc != 0){
                            MPI_Irecv(matleftcolumnOfProc, nombreLignes, MPI_FLOAT, idProc-1, MPI_ANY_TAG, MPI_COMM_WORLD, &leftrequest);
                            initrecvleft = 1;
                        }
                        if(initrecvright == 0 && idProc != vnbProcUse-1){
                            MPI_Irecv(matrightcolumnOfProc, nombreLignes, MPI_FLOAT, idProc+1, MPI_ANY_TAG, MPI_COMM_WORLD, &rightrequest);
                            initrecvright = 1;
                        }

                        if((nbColumns+1) == 1){

                            if(idProc != 0 && idProc != vnbProcUse-1){

                                if(leftflag == 0){
                                    MPI_Test( &leftrequest, &leftflag, &leftstatus );
                                    if(leftflag == 0 && nbCaseCalculated >= totalCaseToCalculate - ((xsizeend-xsizestart)+1)*2){
                                        MPI_Wait(&leftrequest,&leftstatus);
                                    }
                                }
                                if(rightflag == 0){
                                    MPI_Test( &rightrequest, &rightflag, &rightstatus );
                                    if(rightflag == 0 && nbCaseCalculated >= totalCaseToCalculate - ((xsizeend-xsizestart)+1)*2){
                                        MPI_Wait(&rightrequest,&rightstatus);
                                    }
                                }

                                if(leftflag == 1 && rightflag == 1){
                                    //----------Calculate column START----------
                                    mat_y = 0;
                                    mat_x = 0;
                                    for(int y=xsizestart; y<=xsizeend; y++){
                                        usleep(TEMPS_ATTENTE);
                                        if(y==0 || y==(nombreLignes-1)){
                                           mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                           matleftcolumn[mat_y] = 0.0;
                                        }else{
                                            part1formula = (1-4*tempsDiscretise/tailleSubdivision * tailleSubdivision) * mattocalculate[y][mat_x];
                                            part2formula = (tempsDiscretise/tailleSubdivision * tailleSubdivision);
                                            part3formula = 0;
                                            if(y != 0){part3formula += mattocalculate[y-1][mat_x];}
                                            if(y != nombreLignes - 1){part3formula += mattocalculate[y+1][mat_x];}
                                            if(ysizestart != 0){part3formula += matleftcolumnOfProc[(nombreLignes-1) - y];}
                                            if(ysizeend != nombreColonnes - 1){part3formula += matrightcolumnOfProc[(nombreLignes-1) - y];}
                                            float tmpresult = part1formula + (part2formula * part3formula);
                                            mattocalculatekplus1[mat_y][mat_x] = tmpresult;
                                            matleftcolumn[mat_y] = tmpresult;
                                        }
                                        mat_y++;
                                        nbCaseCalculated++;
                                    }
                                    leftflag = 2;
                                    rightflag = 2;
                                    if(k != nombrePasTemps){
                                        MPI_Isend(matleftcolumn, nombreLignes, MPI_FLOAT, idProc-1, k, MPI_COMM_WORLD, &leftrequest);
                                        MPI_Isend(matleftcolumn, nombreLignes, MPI_FLOAT, idProc+1, k, MPI_COMM_WORLD, &rightrequest);
                                    }
                                    //----------Calculate column END----------
                                }

                            }else if(idProc == 0){ //First processor

                                if(rightflag == 0){
                                    MPI_Test( &rightrequest, &rightflag, &rightstatus );
                                    if(rightflag == 0 && nbCaseCalculated >= totalCaseToCalculate - ((xsizeend-xsizestart)+1)){
                                        MPI_Wait(&rightrequest,&rightstatus);
                                    }
                                }

                                if(rightflag == 1){
                                    //----------Calculate last column START----------
                                    mat_y = 0;
                                    mat_x = 0;
                                    for(int y=xsizestart; y<=xsizeend; y++){
                                        usleep(TEMPS_ATTENTE);
                                        mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                        matrightcolumn[mat_y] = 0.0;
                                        mat_y++;
                                        nbCaseCalculated++;
                                    }
                                    rightflag = 2;
                                    if(k != nombrePasTemps){
                                        MPI_Isend(matrightcolumn, nombreLignes, MPI_FLOAT, idProc+1, k, MPI_COMM_WORLD, &rightrequest);
                                    }
                                    //----------Calculate last column END----------
                                }

                            }else{

                                if(leftflag == 0){
                                    MPI_Test( &leftrequest, &leftflag, &leftstatus );
                                    if(leftflag == 0 && nbCaseCalculated >= totalCaseToCalculate - ((xsizeend-xsizestart)+1)){
                                        MPI_Wait(&leftrequest,&leftstatus);
                                    }
                                }

                                if(leftflag == 1){
                                    //----------Left column START----------
                                    mat_y = 0;
                                    mat_x = 0;
                                    for(int y=xsizestart; y<=xsizeend; y++){
                                        mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                        matleftcolumn[mat_y] = 0.0;
                                        mat_y++;
                                        nbCaseCalculated++;
                                    }
                                    leftflag = 2;
                                    if(k != nombrePasTemps){
                                        MPI_Isend(matleftcolumn, nombreLignes, MPI_FLOAT, idProc-1, k, MPI_COMM_WORLD, &leftrequest);
                                    }
                                    //----------Left column END----------
                                }

                            }

                        }else if((nbColumns+1) == 2){

                            if(idProc != 0 && idProc != vnbProcUse-1){

                                if(leftflag == 0){
                                    MPI_Test( &leftrequest, &leftflag, &leftstatus );
                                    if(leftflag == 0 && nbCaseCalculated >= totalCaseToCalculate - ((xsizeend-xsizestart)+1)*2){
                                        MPI_Wait(&leftrequest,&leftstatus);
                                    }
                                }
                                if(rightflag == 0){
                                    MPI_Test( &rightrequest, &rightflag, &rightstatus );
                                    if(rightflag == 0 && nbCaseCalculated >= totalCaseToCalculate - ((xsizeend-xsizestart)+1)*2){
                                        MPI_Wait(&rightrequest,&rightstatus);
                                    }
                                }

                                if(leftflag == 1){
                                    //----------Calculate column START----------
                                    mat_y = 0;
                                    mat_x = 0;
                                    for(int y=xsizestart; y<=xsizeend; y++){
                                        usleep(TEMPS_ATTENTE);
                                        if(y==0 || y==(nombreLignes-1)){
                                           mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                           matleftcolumn[mat_y] = 0.0;
                                        }else{
                                            part1formula = (1-4*tempsDiscretise/tailleSubdivision * tailleSubdivision) * mattocalculate[y][mat_x];
                                            part2formula = (tempsDiscretise/tailleSubdivision * tailleSubdivision);
                                            part3formula = 0;
                                            if(y != 0){part3formula += mattocalculate[y-1][mat_x];}
                                            if(y != nombreLignes - 1){part3formula += mattocalculate[y+1][mat_x];}
                                            if(ysizestart != 0){part3formula += matleftcolumnOfProc[(nombreLignes-1) - y];}
                                            if(ysizestart != nombreColonnes - 1){part3formula += mattocalculate[y][mat_x+1];}
                                            float tmpresult = part1formula + (part2formula * part3formula);
                                            mattocalculatekplus1[mat_y][mat_x] = tmpresult;
                                            matleftcolumn[mat_y] = tmpresult;
                                        }
                                        mat_y++;
                                        nbCaseCalculated++;
                                    }
                                    leftflag = 2;
                                    if(k != nombrePasTemps){
                                        MPI_Isend(matleftcolumn, nombreLignes, MPI_FLOAT, idProc-1, k, MPI_COMM_WORLD, &leftrequest);
                                    }
                                    //----------Calculate column END----------
                                }
                                //We calculate the new values of the second column
                                if(rightflag == 1){
                                    //----------Calculate column START----------
                                    mat_y = 0;
                                    mat_x = 1; //The single column to calculate
                                    for(int y=xsizestart; y<=xsizeend; y++){
                                        usleep(TEMPS_ATTENTE);
                                        if(y==0 || y==(nombreLignes-1)){
                                           mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                           matleftcolumn[mat_y] = 0.0;
                                        }else{

                                            part1formula = (1-4*tempsDiscretise/tailleSubdivision * tailleSubdivision) * mattocalculate[y][mat_x];
                                            part2formula = (tempsDiscretise/tailleSubdivision * tailleSubdivision);
                                            part3formula = 0;
                                            if(y != 0){part3formula += mattocalculate[y-1][mat_x];}
                                            if(y != nombreLignes - 1){part3formula += mattocalculate[y+1][mat_x];}
                                            if(ysizeend != 0){part3formula += mattocalculate[y][mat_x-1];}
                                            if(ysizeend != nombreColonnes - 1){part3formula += matrightcolumnOfProc[(nombreLignes-1) - y];}
                                            float tmpresult = part1formula + (part2formula * part3formula);
                                            mattocalculatekplus1[mat_y][mat_x] = tmpresult;
                                            matleftcolumn[mat_y] = tmpresult;
                                        }
                                        mat_y++;
                                        nbCaseCalculated++;
                                    }
                                    rightflag = 2;
                                    if(k != nombrePasTemps){

                                        MPI_Isend(matleftcolumn, nombreLignes, MPI_FLOAT, idProc+1, k, MPI_COMM_WORLD, &rightrequest);
                                    }
                                    //----------Calculate column END----------
                                }

                            }else if(idProc == 0){ //First processor

                                if(rightflag == 0){
                                    MPI_Test( &rightrequest, &rightflag, &rightstatus );
                                    if(rightflag == 0 && nbCaseCalculated >= totalCaseToCalculate - ((xsizeend-xsizestart)+1)){
                                        MPI_Wait(&rightrequest,&rightstatus);
                                    }

                                }

                                if(leftflag != 2){
                                    //----------Calculate first column START----------
                                    mat_y = 0;
                                    mat_x = 0; //Column equals zero
                                    for(int y=xsizestart; y<=xsizeend; y++){
                                        mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                        matrightcolumn[mat_y] = 0.0;
                                        mat_y++;
                                        nbCaseCalculated++;
                                    }
                                    leftflag = 2;
                                    //----------Calculate first column END----------
                                }

                                //We calculate the new values of the second column
                                if(rightflag == 1){
                                    //----------Calculate column START----------
                                    mat_y = 0;
                                    mat_x = 1; //The single column to calculate
                                    for(int y=xsizestart; y<=xsizeend; y++){
                                        usleep(TEMPS_ATTENTE);
                                        if(y==0 || y==(nombreLignes-1)){
                                           mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                           matleftcolumn[mat_y] = 0.0;
                                        }else{


                                            part1formula = (1-4*tempsDiscretise/tailleSubdivision * tailleSubdivision) * mattocalculate[y][mat_x];
                                            part2formula = (tempsDiscretise/tailleSubdivision * tailleSubdivision);
                                            part3formula = 0;
                                            if(y != 0){part3formula += mattocalculate[y-1][mat_x];}
                                            if(y != nombreLignes - 1){part3formula += mattocalculate[y+1][mat_x];}
                                            if(ysizeend != 0){part3formula += mattocalculate[y][mat_x-1];}
                                            if(ysizeend != nombreColonnes - 1){part3formula += matrightcolumnOfProc[(nombreLignes-1) - y];}
                                            float tmpresult = part1formula + (part2formula * part3formula);

                                            mattocalculatekplus1[mat_y][mat_x] = tmpresult;

                                            matleftcolumn[mat_y] = tmpresult;
                                        }
                                        mat_y++;
                                        nbCaseCalculated++;
                                    }
                                    rightflag = 2;
                                    if(k != nombrePasTemps){

                                        MPI_Isend(matleftcolumn, nombreLignes, MPI_FLOAT, idProc+1, k, MPI_COMM_WORLD, &rightrequest);
                                    }
                                    //----------Calculate column END----------
                                }

                            }else{ //Last processor

                                if(leftflag == 0){
                                    //Check for the left processor
                                    MPI_Test( &leftrequest, &leftflag, &leftstatus );
                                    if(leftflag == 0 && nbCaseCalculated >= totalCaseToCalculate - ((xsizeend-xsizestart)+1)){
                                        MPI_Wait(&leftrequest,&leftstatus);
                                    }

                                }

                                if(rightflag != 2){
                                    //----------Calculate first column START----------
                                    mat_y = 0;
                                    mat_x = 1; //Last column
                                    for(int y=xsizestart; y<=xsizeend; y++){
                                        mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                        matrightcolumn[mat_y] = 0.0;
                                        mat_y++;
                                        nbCaseCalculated++;
                                    }
                                    rightflag = 2;
                                    //----------Calculate first column END----------
                                }

                                //We calculate the new values of the first column
                                if(leftflag == 1){
                                    //----------Calculate column START----------
                                    mat_y = 0;
                                    mat_x = 0; //The single column to calculate
                                    for(int y=xsizestart; y<=xsizeend; y++){
                                        usleep(TEMPS_ATTENTE);
                                        if(y==0 || y==(nombreLignes-1)){
                                           mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                           matleftcolumn[mat_y] = 0.0;
                                        }else{

                                            part1formula = (1-4*tempsDiscretise/tailleSubdivision * tailleSubdivision) * mattocalculate[y][mat_x];
                                            part2formula = (tempsDiscretise/tailleSubdivision * tailleSubdivision);
                                            part3formula = 0;
                                            if(y != 0){part3formula += mattocalculate[y-1][mat_x];}
                                            if(y != nombreLignes - 1){part3formula += mattocalculate[y+1][mat_x];}
                                            if(ysizestart != 0){part3formula += matleftcolumnOfProc[(nombreLignes-1) - y];}
                                            if(ysizestart != nombreColonnes - 1){part3formula += mattocalculate[y][mat_x+1];}
                                            float tmpresult = part1formula + (part2formula * part3formula);
                                            mattocalculatekplus1[mat_y][mat_x] = tmpresult;
                                            matleftcolumn[mat_y] = tmpresult;
                                        }
                                        mat_y++;
                                        nbCaseCalculated++;
                                    }
                                    leftflag = 2;
                                    if(k != nombrePasTemps){

                                        MPI_Isend(matleftcolumn, nombreLignes, MPI_FLOAT, idProc-1, k, MPI_COMM_WORLD, &leftrequest);
                                    }
                                    //----------Calculate column END----------
                                }
                            }

                        }else{ //If there are three columns or more

                            if(idProc != 0 && idProc != vnbProcUse-1){

                                if(leftflag == 0){
                                    //Check for the left processor
                                    MPI_Test( &leftrequest, &leftflag, &leftstatus );
                                    if(leftflag == 0 && nbCaseCalculated >= totalCaseToCalculate - ((xsizeend-xsizestart)+1)*2){
                                        MPI_Wait(&leftrequest,&leftstatus);
                                    }

                                }
                                if(rightflag == 0){
                                    MPI_Test( &rightrequest, &rightflag, &rightstatus );
                                    if(rightflag == 0 && nbCaseCalculated >= totalCaseToCalculate - ((xsizeend-xsizestart)+1)*2){
                                        MPI_Wait(&rightrequest,&rightstatus);
                                    }

                                }

                                //We calculate the new values of the first column
                                if(leftflag == 1){
                                    //----------Calculate column START----------
                                    mat_y = 0;
                                    mat_x = 0; //The single column to calculate
                                    for(int y=xsizestart; y<=xsizeend; y++){
                                        usleep(TEMPS_ATTENTE);
                                        if(y==0 || y==(nombreLignes-1)){
                                           mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                           matleftcolumn[mat_y] = 0.0;
                                        }else{
                                            part1formula = (1-4*tempsDiscretise/tailleSubdivision * tailleSubdivision) * mattocalculate[y][mat_x];
                                            part2formula = (tempsDiscretise/tailleSubdivision * tailleSubdivision);
                                            part3formula = 0;
                                            if(y != 0){part3formula += mattocalculate[y-1][mat_x];}
                                            if(y != nombreLignes - 1){part3formula += mattocalculate[y+1][mat_x];}
                                            if(ysizestart != 0){part3formula += matleftcolumnOfProc[(nombreLignes-1) - y];}
                                            if(ysizeend != nombreColonnes - 1){part3formula += mattocalculate[y][mat_x+1];}
                                            float tmpresult = part1formula + (part2formula * part3formula);
                                            mattocalculatekplus1[mat_y][mat_x] = tmpresult;
                                            matleftcolumn[mat_y] = tmpresult;
                                        }
                                        mat_y++;
                                        nbCaseCalculated++;
                                    }
                                    leftflag = 2;
                                    if(k != nombrePasTemps){

                                        MPI_Isend(matleftcolumn, nombreLignes, MPI_FLOAT, idProc-1, k, MPI_COMM_WORLD, &leftrequest);
                                    }
                                    //----------Calculate column END----------
                                }
                                //We calculate the new values of the last column
                                if(rightflag == 1){
                                    //----------Calculate column START----------
                                    mat_y = 0;
                                    mat_x = (ysizeend - ysizestart); //The last column to calculate
                                    for(int y=xsizestart; y<=xsizeend; y++){
                                        usleep(TEMPS_ATTENTE);
                                        if(y==0 || y==(nombreLignes-1)){
                                           mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                           matleftcolumn[mat_y] = 0.0;
                                        }else{

                                            part1formula = (1-4*tempsDiscretise/tailleSubdivision * tailleSubdivision) * mattocalculate[y][mat_x];
                                            part2formula = (tempsDiscretise/tailleSubdivision * tailleSubdivision);
                                            part3formula = 0;
                                            if(y != 0){part3formula += mattocalculate[y-1][mat_x];}
                                            if(y != nombreLignes - 1){part3formula += mattocalculate[y+1][mat_x];}
                                            if(ysizeend != 0){part3formula += mattocalculate[y][mat_x-1];}
                                            if(ysizeend != nombreColonnes - 1){part3formula += matrightcolumnOfProc[(nombreLignes-1) - y];}
                                            float tmpresult = part1formula + (part2formula * part3formula);
                                            mattocalculatekplus1[mat_y][mat_x] = tmpresult;
                                            matleftcolumn[mat_y] = tmpresult;
                                        }
                                        mat_y++;
                                        nbCaseCalculated++;
                                    }
                                    rightflag = 2;
                                    if(k != nombrePasTemps){

                                        MPI_Isend(matleftcolumn, nombreLignes, MPI_FLOAT, idProc+1, k, MPI_COMM_WORLD, &rightrequest);
                                    }
                                    //----------Calculate column END----------
                                }

                                if(betweenindex <= (ysizeend - ysizestart -1)){

                                    mat_y = 0;
                                    mat_x = betweenindex; //Second column
                                    for(int y=xsizestart; y<=xsizeend; y++){
                                        usleep(TEMPS_ATTENTE);

                                        if(y == 0 || y == (nombreLignes - 1)){
                                            mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                        }else{
                                            part1formula = (1-4*tempsDiscretise/tailleSubdivision * tailleSubdivision) * mattocalculate[y][mat_x];
                                            part2formula = (tempsDiscretise/tailleSubdivision * tailleSubdivision);
                                            part3formula = 0;
                                            if(y != 0){part3formula += mattocalculate[y-1][mat_x];}
                                            if(y != nombreLignes - 1){part3formula += mattocalculate[y+1][mat_x];}
                                            if(mat_x != 0){part3formula += mattocalculate[y][mat_x-1];}
                                            if(mat_x != nombreColonnes - 1){part3formula += mattocalculate[y][mat_x+1];}
                                            mattocalculatekplus1[mat_y][mat_x] = part1formula + (part2formula * part3formula);
                                        }
                                        nbCaseCalculated++;
                                        mat_y++;
                                    }

                                    betweenindex++;
                                }

                            }else if(idProc == 0){ //First processor

                                if(rightflag == 0){
                                    MPI_Test( &rightrequest, &rightflag, &rightstatus );
                                    if(rightflag == 0 && nbCaseCalculated >= totalCaseToCalculate - ((xsizeend-xsizestart)+1)){
                                        MPI_Wait(&rightrequest,&rightstatus);
                                    }

                                }

                                if(leftflag != 2){
                                    //----------Calculate first column START----------
                                    mat_y = 0;
                                    mat_x = 0; //Column equals zero
                                    for(int y=xsizestart; y<=xsizeend; y++){
                                        usleep(TEMPS_ATTENTE);
                                        mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                        mat_y++;
                                        nbCaseCalculated++;
                                    }
                                    leftflag = 2;
                                    //----------Calculate first column END----------
                                }

                                //We calculate the new values of the second column
                                if(rightflag == 1){
                                    //----------Calculate column START----------
                                    mat_y = 0;
                                    mat_x = (ysizeend - ysizestart); //The last column to calculate
                                    for(int y=xsizestart; y<=xsizeend; y++){
                                        usleep(TEMPS_ATTENTE);

                                        if(y==0 || y==(nombreLignes-1)){
                                           mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                           matrightcolumn[mat_y] = 0.0;
                                        }else{
                                            part1formula = (1-4*tempsDiscretise/tailleSubdivision * tailleSubdivision) * mattocalculate[y][mat_x];
                                            part2formula = (tempsDiscretise/tailleSubdivision * tailleSubdivision);
                                            part3formula = 0;

                                            if(y != 0){part3formula += mattocalculate[y-1][mat_x];}
                                            if(y != nombreLignes - 1){part3formula += mattocalculate[y+1][mat_x];}
                                            if(ysizeend != 0){part3formula += mattocalculate[y][mat_x-1];}
                                            if(ysizeend != nombreColonnes - 1){part3formula += matrightcolumnOfProc[(nombreLignes-1) - y];}
                                            float tmpresult = part1formula + (part2formula * part3formula);
                                            mattocalculatekplus1[mat_y][mat_x] = tmpresult;
                                            matrightcolumn[mat_y] = tmpresult;

                                        }

                                        mat_y++;
                                        nbCaseCalculated++;
                                    }
                                    rightflag = 2;
                                    if(k != nombrePasTemps){

                                        MPI_Isend(matrightcolumn, nombreLignes, MPI_FLOAT, idProc+1, k, MPI_COMM_WORLD, &rightrequest);
                                    }
                                    //----------Calculate column END----------
                                }

                                if(betweenindex <= (ysizeend - ysizestart - 1)){

                                    mat_y = 0;
                                    mat_x = betweenindex; //Second column
                                    for(int y=xsizestart; y<=xsizeend; y++){
                                        usleep(TEMPS_ATTENTE);

                                        if(y == 0 || y == (nombreLignes - 1)){
                                            mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                        }else{
                                            part1formula = (1-4*tempsDiscretise/tailleSubdivision * tailleSubdivision) * mattocalculate[y][mat_x];
                                            part2formula = (tempsDiscretise/tailleSubdivision * tailleSubdivision);
                                            part3formula = 0;
                                            if(y != 0){part3formula += mattocalculate[y-1][mat_x];}
                                            if(y != nombreLignes - 1){part3formula += mattocalculate[y+1][mat_x];}
                                            if(mat_x != 0){part3formula += mattocalculate[y][mat_x-1];}
                                            if(mat_x != nombreColonnes - 1){part3formula += mattocalculate[y][mat_x+1];}
                                            mattocalculatekplus1[mat_y][mat_x] = part1formula + (part2formula * part3formula);
                                        }
                                        nbCaseCalculated++;
                                        mat_y++;
                                    }

                                    betweenindex++;
                                }

                            }else{ //Last processor

                                if(leftflag == 0){
                                    //Check for the left processor
                                    MPI_Test( &leftrequest, &leftflag, &leftstatus );
                                    if(leftflag == 0 && nbCaseCalculated >= totalCaseToCalculate - ((xsizeend-xsizestart)+1)){
                                        MPI_Wait(&leftrequest,&leftstatus);
                                    }

                                }

                                if(rightflag != 2){
                                    //----------Calculate first column START----------
                                    mat_y = 0;
                                    mat_x = (ysizeend - ysizestart); //Last column
                                    for(int y=xsizestart; y<=xsizeend; y++){
                                        usleep(TEMPS_ATTENTE);
                                        mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                        matrightcolumn[mat_y] = 0.0;
                                        mat_y++;
                                        nbCaseCalculated++;
                                    }
                                    rightflag = 2;
                                    //----------Calculate first column END----------
                                }

                                //We calculate the new values of the first column
                                if(leftflag == 1){
                                    //----------Calculate column START----------
                                    mat_y = 0;
                                    mat_x = 0; //The single column to calculate
                                    for(int y=xsizestart; y<=xsizeend; y++){
                                        usleep(TEMPS_ATTENTE);
                                        if(y==0 || y==(nombreLignes-1)){
                                           mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                           matleftcolumn[mat_y] = 0.0;
                                        }else{
                                            part1formula = (1-4*tempsDiscretise/tailleSubdivision * tailleSubdivision) * mattocalculate[y][mat_x];
                                            part2formula = (tempsDiscretise/tailleSubdivision * tailleSubdivision);
                                            part3formula = 0;
                                            if(y != 0){part3formula += mattocalculate[y-1][mat_x];}
                                            if(y != nombreLignes - 1){part3formula += mattocalculate[y+1][mat_x];}
                                            if(ysizestart != 0){part3formula += matleftcolumnOfProc[(nombreLignes-1) - y];}
                                            if(ysizestart != nombreColonnes - 1){part3formula += mattocalculate[y][mat_x+1];}
                                            float tmpresult = part1formula + (part2formula * part3formula);
                                            mattocalculatekplus1[mat_y][mat_x] = tmpresult;
                                            matleftcolumn[mat_y] = tmpresult;
                                        }
                                        mat_y++;
                                        nbCaseCalculated++;
                                    }
                                    leftflag = 2;
                                    if(k != nombrePasTemps){

                                        MPI_Isend(matleftcolumn, nombreLignes, MPI_FLOAT, idProc-1, k, MPI_COMM_WORLD, &leftrequest);
                                    }
                                    //----------Calculate column END----------
                                }

                                if(betweenindex <= (ysizeend - ysizestart - 1)){

                                    mat_y = 0;
                                    mat_x = betweenindex;
                                    for(int y=xsizestart; y<=xsizeend; y++){
                                        usleep(TEMPS_ATTENTE);

                                        if(y == 0 || y == (nombreLignes - 1)){
                                            mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                        }else{
                                            part1formula = (1-4*tempsDiscretise/tailleSubdivision * tailleSubdivision) * mattocalculate[y][mat_x];
                                            part2formula = (tempsDiscretise/tailleSubdivision * tailleSubdivision);
                                            part3formula = 0;
                                            if(y != 0){part3formula += mattocalculate[y-1][mat_x];}
                                            if(y != nombreLignes - 1){part3formula += mattocalculate[y+1][mat_x];}
                                            if(mat_x != 0){part3formula += mattocalculate[y][mat_x-1];}
                                            if(mat_x != nombreColonnes - 1){part3formula += mattocalculate[y][mat_x+1];}
                                            mattocalculatekplus1[mat_y][mat_x] = part1formula + (part2formula * part3formula);
                                        }
                                        nbCaseCalculated++;
                                        mat_y++;
                                    }

                                    betweenindex++;
                                }
                            }

                        }

                    }

                }
                if(k > 1){
                    for(int y=0; y<=(ysizeend - ysizestart); y++){
                        for(int x=(xsizeend - xsizestart) ; x>=0; x--){
                            mattocalculate[x][y] = mattocalculatekplus1[x][y];
                        }
                    }
                }
            }



            //------------------------------------ENVOYE A PROCESSEUR 0------------------------------------
            float *partmattosend;
            //initialise le mémoire de la matrice
            partmattosend = (float *)malloc((totalCaseToCalculate+1)*sizeof(float *));
            int ind=0;
            for(int j=0; j<=(ysizeend - ysizestart); j++){
                for(int i=(xsizeend - xsizestart) ; i>=0; i--){
                    partmattosend[ind] = mattocalculate[i][j];
                    ind++;
                }
            }
            partmattosend[ind] = ysizestart;



            if(idProc != 0){


                MPI_Isend(partmattosend, totalCaseToCalculate+1, MPI_FLOAT, 0, idProc*(nombrePasTemps+1), MPI_COMM_WORLD, &req);

            }else{

                for(int j=0; j<=(ysizeend - ysizestart); j++){
                    for(int i=(xsizeend - xsizestart) ; i>=0; i--){
                        matpar[i][j] = mattocalculate[i][j];
                    }
                }
            }
            //------------------------------------ENVOYER A PROCESSEUR 0------------------------------------

            //Coupure verticale
            if(idProc == 0){
            for(int p=1; p<vnbProcUse; p++){

                int arr_length_recv;
                MPI_Status status;

                MPI_Probe(p, p*(nombrePasTemps+1), MPI_COMM_WORLD, &status);
                MPI_Get_count(&status, MPI_INT, &arr_length_recv);

                float *partmatrecv;
                partmatrecv = (float *)malloc(arr_length_recv*sizeof(float *));

                MPI_Recv(partmatrecv, arr_length_recv, MPI_FLOAT, p, p*(nombrePasTemps+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                int row = nombreLignes - 1;
                int col = partmatrecv[arr_length_recv-1];
                for(int i=0; i<arr_length_recv-1; i++){
                    matpar[row][col] = partmatrecv[i];
                    if(row == 0){
                        row = nombreLignes - 1;
                        col++;
                    }else{
                        row--;
                    }
                }
            }

            printf("Impression de la plaque apres %d iteration (marpar) :\n", nombrePasTemps);
            for(int i=nombreLignes -1; i>=0; i--){
                for(int j=0; j<=nombreColonnes -1; j++){
                    printf("%.5f  ", matpar[i][j]);
                }
                printf("\n");
            }
        }
    }
    }
    else{ //Coupure Horizontale
        if(idProc < hnbProcUse){
            ysizestart = 0;
            ysizeend = (int)hy - 1;
            xsizestart = (int)((int)hx * (int)(idProc)) + floorf(((float)(idProc) * (float)(hx - (int)hx)) * 10000000) / 10000000;
            xsizeend = ((int)((int)hx + (int)((idProc) * (int)hx)) + floorf(((float)(idProc + 1) * (float)(hx - (int)hx)) * 10000000) / 10000000) - 1;



            //--------------------Initialise variables pour le calcule---------------------------------------------
            float **mattocalculate;
            mattocalculate = (float **)malloc((int)((xsizeend - xsizestart) + 1)*sizeof(float *));
            for(int i=0;i<(int)((xsizeend - xsizestart) + 1);i++){
                mattocalculate[i] = (float *)malloc((int)((ysizeend - ysizestart) + 1)*sizeof(float));
            }
            float **mattocalculatekplus1;
            mattocalculatekplus1 = (float **)malloc((int)((xsizeend - xsizestart) + 1)*sizeof(float *));
            for(int i=0;i<(int)((xsizeend - xsizestart) + 1);i++){
                mattocalculatekplus1[i] = (float *)malloc((int)((ysizeend - ysizestart) + 1)*sizeof(float));
            }

            for(int y=0; y<=(ysizeend - ysizestart); y++){
                for(int x=(xsizeend - xsizestart) ; x>=0; x--){
                    mattocalculatekplus1[x][y] = 0;

                }
            }

            //Ligne du haut a envoyer
            float *mattopline;
            mattopline = (float *)malloc((int)((ysizeend - ysizestart) + 1)*sizeof(float *));
            //Colonne de droite a envoyer
            float *matbottomline;
            matbottomline = (float *)malloc((int)((ysizeend - ysizestart) + 1)*sizeof(float *));
            //Colonne de gauche a envoyer
            float *mattoplineOfProc;
            mattoplineOfProc = (float *)malloc((int)((ysizeend - ysizestart) + 1)*sizeof(float *));
            //Colonne de droite a envoyer
            float *matbottomlineOfProc;
            matbottomlineOfProc = (float *)malloc((int)((ysizeend - ysizestart) + 1)*sizeof(float *));

            int mat_y = 0;
            int mat_x = 0;
            float part1formula, part2formula, part3formula;
            int nbCaseCalculated = 0;
            int totalCaseToCalculate = ((ysizeend - ysizestart) + 1) * ((xsizeend - xsizestart) + 1);
            int nbLines = (xsizeend - xsizestart);

            for(int k=1; k<=nombrePasTemps; k++){

                mat_y = 0;
                mat_x = 0;
                MPI_Request toprequest, bottomrequest;
                MPI_Status topstatus, bottomstatus;
                int topflag = 0, bottomflag = 0;
                int betweenindex = 1;
                nbCaseCalculated = 0;
                int initrecvtop = 0, initrecvbottom = 0;
                while(nbCaseCalculated < totalCaseToCalculate){

                    //--------------------Premiere iteration-----------------------------------------------------------------
                    if(k == 1){

                        //----------Top line START----------
                        mat_y = 0;
                        mat_x = 0;
                        for(int x=ysizestart; x<=ysizeend; x++){
                            usleep(TEMPS_ATTENTE);
                            if( idProc == 0 || (idProc == (hnbProcUse - 1) && nbLines == 0) || x==0 || x==(nombreColonnes-1)){
                               mattocalculate[mat_y][mat_x] = 0.0;
                               mattopline[mat_x] = 0.0;
                            }else{
                                part1formula = (1-4*tempsDiscretise/tailleSubdivision * tailleSubdivision) * matpar[xsizestart][x];
                                part2formula = (tempsDiscretise/tailleSubdivision * tailleSubdivision);
                                part3formula = 0;
                                if(xsizestart != 0){part3formula += matpar[xsizestart-1][x];}
                                if(xsizestart != nombreLignes - 1){part3formula += matpar[xsizestart+1][x];}
                                if(x != 0){part3formula += matpar[xsizestart][x-1];}
                                if(x != nombreColonnes - 1){part3formula += matpar[xsizestart][x+1];}
                                float tmpresult = part1formula + (part2formula * part3formula);
                                mattocalculate[mat_y][mat_x] = tmpresult;
                                mattopline[mat_x] = tmpresult;

                            }
                            mat_x++;
                            nbCaseCalculated++;
                        }

                        if(idProc != 0 && nombrePasTemps > 1){

                            MPI_Isend(mattopline, nombreColonnes, MPI_FLOAT, idProc-1, k, MPI_COMM_WORLD, &toprequest);
                        }
                        //----------Top line END----------

                        //----------Bottom line START----------

                        if(nbLines > 0){
                            mat_y = nbLines;
                            mat_x = 0;
                            for(int x=ysizestart; x<=ysizeend; x++){

                                usleep(TEMPS_ATTENTE);
                                if(idProc == (hnbProcUse - 1) || x==0 || x==(nombreColonnes-1)){

                                   mattocalculate[mat_y][mat_x] = 0.0;
                                   matbottomline[mat_x] = 0.0;
                                }else{

                                    part1formula = (1-4*tempsDiscretise/tailleSubdivision * tailleSubdivision) * matpar[xsizeend][x];
                                    part2formula = (tempsDiscretise/tailleSubdivision * tailleSubdivision);
                                    part3formula = 0;

                                    if(xsizeend != 0){part3formula += matpar[xsizeend-1][x];}
                                    if(xsizeend != nombreLignes - 1){part3formula += matpar[xsizeend+1][x];}
                                    if(x != 0){part3formula += matpar[xsizeend][x-1];}
                                    if(x != nombreColonnes - 1){part3formula += matpar[xsizeend][x+1];}
                                    float tmpresult = part1formula + (part2formula * part3formula);
                                    mattocalculate[mat_y][mat_x] = tmpresult;
                                    matbottomline[mat_x] = tmpresult;
                                }
                                mat_x++;
                                nbCaseCalculated++;
                            }

                            if(idProc != (hnbProcUse-1)  && nombrePasTemps > 1){

                                MPI_Isend(matbottomline, nombreColonnes, MPI_FLOAT, idProc+1, k, MPI_COMM_WORLD, &bottomrequest);
                            }
                        }else{
                            if(idProc != (hnbProcUse-1)  && nombrePasTemps > 1){

                                MPI_Isend(mattopline, nombreColonnes, MPI_FLOAT, idProc+1, k, MPI_COMM_WORLD, &bottomrequest);
                            }
                        }
                        //----------Bottom line END----------

                        //----------Between columns START----------
                        if(nbLines > 1){
                            mat_y = 1;
                            mat_x = 0;
                            for(int y=xsizestart+1; y<=xsizeend-1; y++){
                                mat_x = 0;
                                for(int x=ysizestart; x<=ysizeend; x++){
                                    usleep(TEMPS_ATTENTE);
                                    if(x == 0 || x == (nombreColonnes - 1)){
                                        mattocalculate[mat_y][mat_x] = 0.0;
                                    }else{
                                        part1formula = (1-4*tempsDiscretise/tailleSubdivision * tailleSubdivision) * matpar[y][x];
                                        part2formula = (tempsDiscretise/tailleSubdivision * tailleSubdivision);
                                        part3formula = 0;
                                        if(y != 0){part3formula += matpar[y-1][x];}
                                        if(y != nombreLignes - 1){part3formula += matpar[y+1][x];}
                                        if(x != 0){part3formula += matpar[y][x-1];}
                                        if(x != nombreColonnes - 1){part3formula += matpar[y][x+1];}
                                        mattocalculate[mat_y][mat_x] = part1formula + (part2formula * part3formula);
                                    }
                                    mat_x++;
                                    nbCaseCalculated++;
                                }
                                mat_y++;
                            }
                        }
                        //----------Between columns END----------

                    //--------------------Tout autres iterations----------------------------------------------------------------
                    }else{
                        //Initialisation du Ireceive
                        if(initrecvtop == 0 && idProc != 0){

                            MPI_Irecv(mattoplineOfProc, nombreColonnes, MPI_FLOAT, idProc-1, k-1, MPI_COMM_WORLD, &toprequest);
                            initrecvtop = 1;
                        }
                        if(initrecvbottom == 0 && idProc != hnbProcUse-1){

                            MPI_Irecv(matbottomlineOfProc, nombreColonnes, MPI_FLOAT, idProc+1, k-1, MPI_COMM_WORLD, &bottomrequest);
                            initrecvbottom = 1;
                        }

                        if((nbLines+1) == 1){

                            if(idProc != 0 && idProc != hnbProcUse-1){
                                if(topflag == 0){

                                    MPI_Test( &toprequest, &topflag, &topstatus );
                                    if(topflag == 0 && nbCaseCalculated >= totalCaseToCalculate - ((ysizeend-ysizestart)+1)){
                                        MPI_Wait(&toprequest,&topstatus);
                                    }

                                }
                                if(bottomflag == 0){
                                    MPI_Test( &bottomrequest, &bottomflag, &bottomstatus );
                                    if(bottomflag == 0 && nbCaseCalculated >= totalCaseToCalculate - ((ysizeend-ysizestart)+1)){
                                        MPI_Wait(&bottomrequest,&bottomstatus);
                                    }

                                }

                                if(topflag == 1 && bottomflag == 1){
                                    //----------Calculate column START----------
                                    mat_y = 0;
                                    mat_x = 0;
                                    for(int x=ysizestart; x<=ysizeend; x++){
                                        usleep(TEMPS_ATTENTE);
                                        if(x==0 || x==(nombreColonnes-1)){
                                           mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                           mattopline[mat_x] = 0.0;
                                        }else{

                                            part1formula = (1-4*tempsDiscretise/tailleSubdivision * tailleSubdivision) * mattocalculate[mat_y][x];
                                            part2formula = (tempsDiscretise/tailleSubdivision * tailleSubdivision);
                                            part3formula = 0;
                                            if(xsizestart != 0){part3formula += mattoplineOfProc[x];}
                                            if(xsizeend != nombreLignes - 1){part3formula += matbottomlineOfProc[x];}
                                            if(x != 0){part3formula += mattocalculate[mat_y][x-1];}
                                            if(x != nombreColonnes - 1){part3formula += mattocalculate[mat_y][x+1];}
                                            float tmpresult = part1formula + (part2formula * part3formula);
                                            mattocalculatekplus1[mat_y][mat_x] = tmpresult;
                                            mattopline[mat_x] = tmpresult;
                                        }
                                        mat_x++;
                                        nbCaseCalculated++;
                                    }
                                    topflag = 2;
                                    bottomflag = 2;
                                    if(k != nombrePasTemps){

                                        MPI_Isend(mattopline, nombreColonnes, MPI_FLOAT, idProc-1, k, MPI_COMM_WORLD, &toprequest);
                                        MPI_Isend(mattopline, nombreColonnes, MPI_FLOAT, idProc+1, k, MPI_COMM_WORLD, &bottomrequest);
                                    }
                                    //----------Calculate column END----------
                                }

                            }else if(idProc == 0){
                                if(bottomflag == 0){
                                    MPI_Test( &bottomrequest, &bottomflag, &bottomstatus );
                                    if(bottomflag == 0 && nbCaseCalculated >= totalCaseToCalculate - ((ysizeend-ysizestart)+1)){
                                        MPI_Wait(&bottomrequest,&bottomstatus);
                                    }
                                }

                                if(bottomflag == 1){
                                    //----------Calculate last line START----------
                                    mat_y = 0;
                                    mat_x = 0;
                                    for(int x=ysizestart; x<=ysizeend; x++){
                                        usleep(TEMPS_ATTENTE);
                                        mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                        matbottomline[mat_x] = 0.0;
                                        mat_x++;
                                        nbCaseCalculated++;
                                    }
                                    bottomflag = 2;
                                    if(k != nombrePasTemps){

                                        MPI_Isend(matbottomline, nombreColonnes, MPI_FLOAT, idProc+1, k, MPI_COMM_WORLD, &bottomrequest);
                                    }
                                    //----------Calculate last line END----------
                                }

                            }else{

                                if(topflag == 0){
                                    MPI_Test( &toprequest, &topflag, &topstatus);
                                    if(topflag == 0 && nbCaseCalculated >= totalCaseToCalculate - ((ysizeend-ysizestart)+1)){
                                        MPI_Wait(&toprequest,&topstatus);
                                    }

                                }

                                if(topflag == 1){
                                    //----------Left column START----------
                                    mat_y = 0;
                                    mat_x = 0;
                                    for(int x=ysizestart; x<=ysizeend; x++){
                                        mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                        mattopline[mat_x] = 0.0;
                                        mat_x++;
                                        nbCaseCalculated++;
                                    }
                                    topflag = 2;
                                    if(k != nombrePasTemps){

                                        MPI_Isend(mattopline, nombreColonnes, MPI_FLOAT, idProc-1, k, MPI_COMM_WORLD, &toprequest);
                                    }
                                    //----------Left column END----------
                                }

                            }
                        }else if((nbLines+1) == 2){

                            if(idProc != 0 && idProc != hnbProcUse-1){

                                if(topflag == 0){
                                    MPI_Test( &toprequest, &topflag, &topstatus );
                                    if(topflag == 0 && nbCaseCalculated >= totalCaseToCalculate - ((ysizeend-ysizestart)+1)*2){
                                        MPI_Wait(&toprequest,&topstatus);
                                    }

                                }
                                if(bottomflag == 0){
                                    MPI_Test( &bottomrequest, &bottomflag, &bottomstatus );
                                    if(bottomflag == 0 && nbCaseCalculated >= totalCaseToCalculate - ((ysizeend-ysizestart)+1)*2){
                                        MPI_Wait(&bottomrequest,&bottomstatus);
                                    }

                                }


                                if(topflag == 1){
                                    //----------Calculate line START----------
                                    mat_y = 0;
                                    mat_x = 0;
                                    for(int x=ysizestart; x<=ysizeend; x++){
                                        usleep(TEMPS_ATTENTE);
                                        if(x==0 || x==(nombreColonnes-1)){
                                           mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                           mattopline[mat_x] = 0.0;
                                        }else{

                                            part1formula = (1-4*tempsDiscretise/tailleSubdivision * tailleSubdivision) * mattocalculate[mat_y][x];
                                            part2formula = (tempsDiscretise/tailleSubdivision * tailleSubdivision);
                                            part3formula = 0;
                                            if(xsizestart != 0){part3formula += mattoplineOfProc[x];}
                                            if(xsizestart != nombreLignes - 1){part3formula += mattocalculate[mat_y+1][x];}
                                            if(x != 0){part3formula += mattocalculate[mat_y][x-1];}
                                            if(x != nombreColonnes - 1){part3formula += mattocalculate[mat_y][x+1];}
                                            float tmpresult = part1formula + (part2formula * part3formula);
                                            mattocalculatekplus1[mat_y][mat_x] = tmpresult;
                                            mattopline[mat_x] = tmpresult;
                                        }
                                        mat_x++;
                                        nbCaseCalculated++;
                                    }
                                    topflag = 2;
                                    if(k != nombrePasTemps){

                                        MPI_Isend(mattopline, nombreColonnes, MPI_FLOAT, idProc-1, k, MPI_COMM_WORLD, &toprequest);
                                    }
                                    //----------Calculate line END----------
                                }

                                if(bottomflag == 1){
                                    //----------Calculate line START----------
                                    mat_y = 1;
                                    mat_x = 0;
                                    for(int x=ysizestart; x<=ysizeend; x++){
                                        usleep(TEMPS_ATTENTE);
                                        if(x==0 || x==(nombreColonnes-1)){
                                           mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                           matbottomline[mat_x] = 0.0;
                                        }else{

                                            part1formula = (1-4*tempsDiscretise/tailleSubdivision * tailleSubdivision) * mattocalculate[mat_y][x];
                                            part2formula = (tempsDiscretise/tailleSubdivision * tailleSubdivision);
                                            part3formula = 0;
                                            if(xsizeend != 0){part3formula += mattocalculate[mat_y-1][x];}
                                            if(xsizeend != nombreLignes - 1){part3formula += matbottomlineOfProc[x];}
                                            if(x != 0){part3formula += mattocalculate[mat_y][x-1];}
                                            if(x != nombreColonnes - 1){part3formula += mattocalculate[mat_y][x+1];}
                                            float tmpresult = part1formula + (part2formula * part3formula);
                                            mattocalculatekplus1[mat_y][mat_x] = tmpresult;
                                            matbottomline[mat_x] = tmpresult;
                                        }
                                        mat_x++;
                                        nbCaseCalculated++;
                                    }
                                    bottomflag = 2;
                                    if(k != nombrePasTemps){

                                        MPI_Isend(matbottomline, nombreColonnes, MPI_FLOAT, idProc+1, k, MPI_COMM_WORLD, &bottomrequest);
                                    }
                                    //----------Calculate line END----------
                                }

                            }else if(idProc == 0){

                                if(bottomflag == 0){

                                    MPI_Test( &bottomrequest, &bottomflag, &bottomstatus );
                                    if(bottomflag == 0 && nbCaseCalculated >= totalCaseToCalculate - ((ysizeend-ysizestart)+1)){
                                        MPI_Wait(&bottomrequest,&bottomstatus);
                                    }

                                }

                                if(topflag != 2){
                                    //----------Calculate first column START----------
                                    mat_y = 0;
                                    mat_x = 0;
                                    for(int x=ysizestart; x<=ysizeend; x++){
                                        mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                        mat_x++;
                                        nbCaseCalculated++;
                                    }
                                    topflag = 2;
                                    //----------Calculate first column END----------
                                }


                                if(bottomflag == 1){
                                    //----------Calculate column START----------
                                    mat_y = 1;
                                    mat_x = 0;
                                    for(int x=ysizestart; x<=ysizeend; x++){
                                        usleep(TEMPS_ATTENTE);
                                        if(x==0 || x==(nombreColonnes-1)){
                                           mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                           matbottomline[mat_x] = 0.0;
                                        }else{

                                            part1formula = (1-4*tempsDiscretise/tailleSubdivision * tailleSubdivision) * mattocalculate[mat_y][x];
                                            part2formula = (tempsDiscretise/tailleSubdivision * tailleSubdivision);
                                            part3formula = 0;
                                            if(xsizeend != 0){part3formula += mattocalculate[mat_y-1][x];}
                                            if(xsizeend != nombreLignes - 1){part3formula += matbottomlineOfProc[x];}
                                            if(x != 0){part3formula += mattocalculate[mat_y][x-1];}
                                            if(x != nombreColonnes - 1){part3formula += mattocalculate[mat_y][x+1];}
                                            float tmpresult = part1formula + (part2formula * part3formula);

                                            mattocalculatekplus1[mat_y][mat_x] = tmpresult;

                                            matbottomline[mat_x] = tmpresult;
                                        }
                                        mat_x++;
                                        nbCaseCalculated++;
                                    }
                                    bottomflag = 2;
                                    if(k != nombrePasTemps){

                                        MPI_Isend(matbottomline, nombreColonnes, MPI_FLOAT, idProc+1, k, MPI_COMM_WORLD, &bottomrequest);
                                    }
                                    //----------Calculate column END----------
                                }

                            }else{

                                if(topflag == 0){
                                    MPI_Test( &toprequest, &topflag, &topstatus );
                                    if(topflag == 0 && nbCaseCalculated >= totalCaseToCalculate - ((ysizeend-ysizestart)+1)){
                                        MPI_Wait(&toprequest,&topstatus);
                                    }

                                }

                                if(bottomflag != 2){
                                    //----------Calculate last line START----------
                                    mat_y = 1;
                                    mat_x = 0;
                                    for(int x=ysizestart; x<=ysizeend; x++){
                                        mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                        mat_x++;
                                        nbCaseCalculated++;
                                    }
                                    bottomflag = 2;
                                    //----------Calculate last line END----------
                                }

                                //We calculate the new values of the first column
                                if(topflag == 1){
                                    //----------Calculate column START----------
                                    mat_y = 0;
                                    mat_x = 0;
                                    for(int x=ysizestart; x<=ysizeend; x++){
                                        usleep(TEMPS_ATTENTE);
                                        if(x==0 || x==(nombreColonnes-1)){
                                           mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                           mattopline[mat_x] = 0.0;
                                        }else{

                                            part1formula = (1-4*tempsDiscretise/tailleSubdivision * tailleSubdivision) * mattocalculate[mat_y][x];
                                            part2formula = (tempsDiscretise/tailleSubdivision * tailleSubdivision);
                                            part3formula = 0;
                                            if(xsizestart != 0){part3formula += mattoplineOfProc[x];}
                                            if(xsizestart != nombreLignes - 1){part3formula += mattocalculate[mat_y+1][x];}
                                            if(x != 0){part3formula += mattocalculate[mat_y][x-1];}
                                            if(x != nombreColonnes - 1){part3formula += mattocalculate[mat_y][x+1];}
                                            float tmpresult = part1formula + (part2formula * part3formula);
                                            mattocalculatekplus1[mat_y][mat_x] = tmpresult;
                                            mattopline[mat_x] = tmpresult;
                                        }
                                        mat_x++;
                                        nbCaseCalculated++;
                                    }
                                    topflag = 2;
                                    if(k != nombrePasTemps){

                                        MPI_Isend(mattopline, nombreColonnes, MPI_FLOAT, idProc-1, k, MPI_COMM_WORLD, &toprequest);
                                    }
                                    //----------Calculate column END----------
                                }
                            }

                        }else{ //If there are three columns or more

                            if(idProc != 0 && idProc != hnbProcUse-1){

                                if(topflag == 0){
                                    MPI_Test( &toprequest, &topflag, &topstatus );
                                    if(topflag == 0 && nbCaseCalculated >= totalCaseToCalculate - ((ysizeend-ysizestart)+1)*2){
                                        MPI_Wait(&toprequest,&topstatus);
                                    }

                                }
                                if(bottomflag == 0){
                                    MPI_Test( &bottomrequest, &bottomflag, &bottomstatus );
                                    if(bottomflag == 0 && nbCaseCalculated >= totalCaseToCalculate - ((ysizeend-ysizestart)+1)*2){
                                        MPI_Wait(&bottomrequest,&bottomstatus);
                                    }

                                }


                                if(topflag == 1){
                                    //----------Calculate line START----------
                                    mat_y = 0;
                                    mat_x = 0;
                                    for(int x=ysizestart; x<=ysizeend; x++){
                                        usleep(TEMPS_ATTENTE);
                                        if(x==0 || x==(nombreColonnes-1)){
                                           mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                           mattopline[mat_x] = 0.0;
                                        }else{
                                            part1formula = (1-4*tempsDiscretise/tailleSubdivision * tailleSubdivision) * mattocalculate[mat_y][x];
                                            part2formula = (tempsDiscretise/tailleSubdivision * tailleSubdivision);
                                            part3formula = 0;
                                            if(xsizestart != 0){part3formula += mattoplineOfProc[x];}
                                            if(xsizeend != nombreLignes - 1){part3formula += mattocalculate[mat_y+1][x];}
                                            if(x != 0){part3formula += mattocalculate[mat_y][x-1];}
                                            if(x != nombreColonnes - 1){part3formula += mattocalculate[mat_y][x+1];}
                                            float tmpresult = part1formula + (part2formula * part3formula);
                                            mattocalculatekplus1[mat_y][mat_x] = tmpresult;
                                            mattopline[mat_x] = tmpresult;
                                        }
                                        mat_x++;
                                        nbCaseCalculated++;
                                    }
                                    topflag = 2;
                                    if(k != nombrePasTemps){

                                        MPI_Isend(mattopline, nombreColonnes, MPI_FLOAT, idProc-1, k, MPI_COMM_WORLD, &toprequest);
                                    }
                                    //----------Calculate line END----------
                                }

                                if(bottomflag == 1){
                                    //----------Calculate line START----------
                                    mat_y = (xsizeend - xsizestart);
                                    mat_x = 0;
                                    for(int x=ysizestart; x<=ysizeend; x++){
                                        usleep(TEMPS_ATTENTE);
                                        if(x==0 || x==(nombreColonnes-1)){
                                           mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                           matbottomline[mat_x] = 0.0;
                                        }else{

                                            part1formula = (1-4*tempsDiscretise/tailleSubdivision * tailleSubdivision) * mattocalculate[mat_y][x];
                                            part2formula = (tempsDiscretise/tailleSubdivision * tailleSubdivision);
                                            part3formula = 0;
                                            if(xsizeend != 0){part3formula += mattocalculate[mat_y-1][x];}
                                            if(xsizeend != nombreLignes - 1){part3formula += matbottomlineOfProc[x];}
                                            if(x != 0){part3formula += mattocalculate[mat_y][x-1];}
                                            if(x != nombreColonnes - 1){part3formula += mattocalculate[mat_y][x+1];}
                                            float tmpresult = part1formula + (part2formula * part3formula);
                                            mattocalculatekplus1[mat_y][mat_x] = tmpresult;
                                            matbottomline[mat_x] = tmpresult;
                                        }
                                        mat_x++;
                                        nbCaseCalculated++;
                                    }
                                    bottomflag = 2;
                                    if(k != nombrePasTemps){

                                        MPI_Isend(matbottomline, nombreColonnes, MPI_FLOAT, idProc+1, k, MPI_COMM_WORLD, &bottomrequest);
                                    }
                                    //----------Calculate column END----------
                                }

                                if(betweenindex <= (xsizeend - xsizestart -1)){

                                    mat_y = betweenindex;
                                    mat_x = 0;
                                    for(int x=ysizestart; x<=ysizeend; x++){
                                        usleep(TEMPS_ATTENTE);

                                        if(x == 0 || x == (nombreColonnes - 1)){
                                            mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                        }else{
                                            part1formula = (1-4*tempsDiscretise/tailleSubdivision * tailleSubdivision) * mattocalculate[mat_y][x];
                                            part2formula = (tempsDiscretise/tailleSubdivision * tailleSubdivision);
                                            part3formula = 0;
                                            if(mat_y != 0){part3formula += mattocalculate[mat_y-1][x];}
                                            if(mat_y != nombreLignes - 1){part3formula += mattocalculate[mat_y+1][x];}
                                            if(x != 0){part3formula += mattocalculate[mat_y][x-1];}
                                            if(x != nombreColonnes - 1){part3formula += mattocalculate[mat_y][x+1];}
                                            mattocalculatekplus1[mat_y][mat_x] = part1formula + (part2formula * part3formula);
                                        }
                                        nbCaseCalculated++;
                                        mat_x++;
                                    }

                                    betweenindex++;
                                }

                            }else if(idProc == 0){

                                if(bottomflag == 0){
                                    MPI_Test( &bottomrequest, &bottomflag, &bottomstatus );
                                    if(bottomflag == 0 && nbCaseCalculated >= totalCaseToCalculate - ((ysizeend-ysizestart)+1)){
                                        MPI_Wait(&bottomrequest,&bottomstatus);
                                    }

                                }

                                if(topflag != 2){
                                    //----------Calculate first line START----------
                                    mat_y = 0;
                                    mat_x = 0;
                                    for(int x=ysizestart; x<=ysizeend; x++){
                                        usleep(TEMPS_ATTENTE);
                                        mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                        mat_x++;
                                        nbCaseCalculated++;
                                    }
                                    topflag = 2;
                                    //----------Calculate first column END----------
                                }


                                if(bottomflag == 1){
                                    //----------Calculate column START----------
                                    mat_y = (xsizeend - xsizestart);
                                    mat_x = 0;
                                    for(int x=ysizestart; x<=ysizeend; x++){
                                        usleep(TEMPS_ATTENTE);

                                        if(x==0 || x==(nombreColonnes-1)){
                                           mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                           matbottomline[mat_x] = 0.0;
                                        }else{
                                            part1formula = (1-4*tempsDiscretise/tailleSubdivision * tailleSubdivision) * mattocalculate[mat_y][x];
                                            part2formula = (tempsDiscretise/tailleSubdivision * tailleSubdivision);
                                            part3formula = 0;

                                            if(xsizeend != 0){part3formula += mattocalculate[mat_y-1][x];}
                                            if(xsizeend != nombreLignes - 1){part3formula += matbottomlineOfProc[x];}
                                            if(x != 0){part3formula += mattocalculate[mat_y][x-1];}
                                            if(x != nombreColonnes - 1){part3formula += mattocalculate[mat_y][x+1];}
                                            float tmpresult = part1formula + (part2formula * part3formula);
                                            mattocalculatekplus1[mat_y][mat_x] = tmpresult;
                                            matbottomline[mat_x] = tmpresult;

                                        }

                                        mat_x++;
                                        nbCaseCalculated++;
                                    }
                                    bottomflag = 2;
                                    if(k != nombrePasTemps){
                                        MPI_Isend(matbottomline, nombreColonnes, MPI_FLOAT, idProc+1, k, MPI_COMM_WORLD, &bottomrequest);
                                    }
                                    //----------Calculate column END----------
                                }

                                if(betweenindex <= (xsizeend - xsizestart - 1)){

                                    mat_y = betweenindex;
                                    mat_x = 0;
                                    for(int x=ysizestart; x<=ysizeend; x++){
                                        usleep(TEMPS_ATTENTE);

                                        if(x == 0 || x == (nombreColonnes - 1)){
                                            mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                        }else{
                                            part1formula = (1-4*tempsDiscretise/tailleSubdivision * tailleSubdivision) * mattocalculate[mat_y][x];
                                            part2formula = (tempsDiscretise/tailleSubdivision * tailleSubdivision);
                                            part3formula = 0;
                                            if(mat_y != 0){part3formula += mattocalculate[mat_y-1][x];}
                                            if(mat_y != nombreLignes - 1){part3formula += mattocalculate[mat_y+1][x];}
                                            if(x != 0){part3formula += mattocalculate[mat_y][x-1];}
                                            if(x != nombreColonnes - 1){part3formula += mattocalculate[mat_y][x+1];}
                                            mattocalculatekplus1[mat_y][mat_x] = part1formula + (part2formula * part3formula);
                                        }
                                        nbCaseCalculated++;
                                        mat_x++;
                                    }

                                    betweenindex++;
                                }

                            }else{
                                if(topflag == 0){
                                    MPI_Test( &toprequest, &topflag, &topstatus);
                                    if(topflag == 0 && nbCaseCalculated >= totalCaseToCalculate - ((ysizeend-ysizestart)+1)){
                                        MPI_Wait(&toprequest,&topstatus);
                                    }

                                }

                                if(bottomflag != 2){
                                    //----------Calculate last line START----------
                                    mat_y = (xsizeend - xsizestart);
                                    mat_x = 0;
                                    for(int x=ysizestart; x<=ysizeend; x++){
                                        usleep(TEMPS_ATTENTE);
                                        mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                        mat_x++;
                                        nbCaseCalculated++;
                                    }
                                    bottomflag = 2;
                                    //----------Calculate last line END----------
                                }


                                if(topflag == 1){
                                    //----------Calculate column START----------
                                    mat_y = 0;
                                    mat_x = 0;
                                    for(int x=ysizestart; x<=ysizeend; x++){
                                        usleep(TEMPS_ATTENTE);
                                        if(x==0 || x==(nombreColonnes-1)){
                                           mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                           mattopline[mat_x] = 0.0;
                                        }else{
                                            part1formula = (1-4*tempsDiscretise/tailleSubdivision * tailleSubdivision) * mattocalculate[mat_y][x];
                                            part2formula = (tempsDiscretise/tailleSubdivision * tailleSubdivision);
                                            part3formula = 0;
                                            if(xsizestart != 0){part3formula += mattoplineOfProc[x];}
                                            if(xsizestart != nombreLignes - 1){part3formula += mattocalculate[mat_y+1][x];}
                                            if(x != 0){part3formula += mattocalculate[mat_y][x-1];}
                                            if(x != nombreColonnes - 1){part3formula += mattocalculate[mat_y][x+1];}
                                            float tmpresult = part1formula + (part2formula * part3formula);
                                            mattocalculatekplus1[mat_y][mat_x] = tmpresult;
                                            mattopline[mat_x] = tmpresult;
                                        }
                                        mat_x++;
                                        nbCaseCalculated++;
                                    }
                                    topflag = 2;
                                    if(k != nombrePasTemps){

                                        MPI_Isend(mattopline, nombreColonnes, MPI_FLOAT, idProc-1, k, MPI_COMM_WORLD, &toprequest);
                                    }
                                    //----------Calculate column END----------
                                }

                                if(betweenindex <= (xsizeend - xsizestart - 1)){

                                    mat_y = betweenindex;
                                    mat_x = 0;
                                    for(int x=ysizestart; x<=ysizeend; x++){
                                        usleep(TEMPS_ATTENTE);

                                        if(x == 0 || x == (nombreColonnes - 1)){
                                            mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                        }else{
                                            part1formula = (1-4*tempsDiscretise/tailleSubdivision * tailleSubdivision) * mattocalculate[mat_y][x];
                                            part2formula = (tempsDiscretise/tailleSubdivision * tailleSubdivision);
                                            part3formula = 0;
                                            if(mat_y != 0){part3formula += mattocalculate[mat_y-1][x];}
                                            if(mat_y != nombreLignes - 1){part3formula += mattocalculate[mat_y+1][x];}
                                            if(x != 0){part3formula += mattocalculate[mat_y][x-1];}
                                            if(x != nombreColonnes - 1){part3formula += mattocalculate[mat_y][x+1];}
                                            mattocalculatekplus1[mat_y][mat_x] = part1formula + (part2formula * part3formula);
                                        }
                                        nbCaseCalculated++;
                                        mat_x++;
                                    }
                                    betweenindex++;
                                }
                            }
                        }
                    }
                }

                //S'il y a plus d'une iteration
                if(k > 1){
                    for(int y=0; y<=(ysizeend - ysizestart); y++){
                        for(int x=(xsizeend - xsizestart) ; x>=0; x--){
                            mattocalculate[x][y] = mattocalculatekplus1[x][y];
                        }
                    }
                }
            }


            //------------------------------------ENVOYER A PROCESSEUR 0------------------------------------
            float *partmattosend;

            partmattosend = (float *)malloc((totalCaseToCalculate+1)*sizeof(float *));
            int ind=0;
            for(int i=0 ; i<=(xsizeend - xsizestart); i++){
                for(int j=0; j<=(ysizeend - ysizestart); j++){
                    partmattosend[ind] = mattocalculate[i][j];

                    ind++;
                }
            }
            partmattosend[ind] = xsizestart;



            if(idProc != 0){


                MPI_Isend(partmattosend, totalCaseToCalculate+1, MPI_FLOAT, 0, idProc*(nombrePasTemps+1), MPI_COMM_WORLD, &req);

            }else{

                for(int j=0; j<=(ysizeend - ysizestart); j++){
                    for(int i=(xsizeend - xsizestart) ; i>=0; i--){
                        matpar[i][j] = mattocalculate[i][j];
                    }
                }
            }
            //------------------------------------ENVOYER A PROCESSEUR 0------------------------------------

            //Coupure horizontale
            if(idProc == 0){
            for(int p=1; p<hnbProcUse; p++){

                int arr_length_recv;
                MPI_Status status;


                MPI_Probe(p, p*(nombrePasTemps+1), MPI_COMM_WORLD, &status);
                MPI_Get_count(&status, MPI_INT, &arr_length_recv);

                float *partmatrecv;
                partmatrecv = (float *)malloc(arr_length_recv*sizeof(float *));

                MPI_Recv(partmatrecv, arr_length_recv, MPI_FLOAT, p, p*(nombrePasTemps+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                int row = partmatrecv[arr_length_recv-1];
                int col = nombreColonnes - 1;
                for(int i=0; i<arr_length_recv-1; i++){
                    matpar[row][col] = partmatrecv[i];
                    if(col == 0){
                        row++;
                        col = nombreColonnes - 1;
                    }else{
                        col--;
                    }
                }
            }

            printf("Impression de la plaque apres %d iteration (marpar) :\n", nombrePasTemps);
            for(int i=nombreLignes -1; i>=0; i--){
                for(int j=0; j<=nombreColonnes -1; j++){
                    printf("%.5f  ", matpar[i][j]);
                }
                printf("\n");
            }
        }
        }
    }
    gettimeofday (&tp, NULL);
    timeEnd = (double) (tp.tv_sec) + (double) (tp.tv_usec) / 1e6;
    Texec = timeEnd - timeStart;

    if(idProc == 0){
        printf("Temps d'execution parallèle : %f\n", Texec);
        printf("L'accélération: %f\n",TexecSeq/Texec);
    }

    MPI_Finalize();
    return 0;
}
