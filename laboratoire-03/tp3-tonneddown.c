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

    float **matseq;
    //initialise le mémoire de la matrice
    matseq = (float **)malloc(nombreLignes*sizeof(float *));
    for(int i=0;i<nombreLignes;i++){
        matseq[i] = (float *)malloc(nombreColonnes*sizeof(float));
    }
    float **matseqkplus1;
    //initialise le mémoire de la matrice
    matseqkplus1 = (float **)malloc(nombreLignes*sizeof(float *));
    for(int i=0;i<nombreLignes;i++){
        matseqkplus1[i] = (float *)malloc(nombreColonnes*sizeof(float));
    }

    //Initialisation de la matrice
    for(int y=0; y<nombreLignes; y++){
        for(int x=0; x<nombreColonnes; x++){
            usleep(TEMPS_ATTENTE);
            matseq[y][x] = (x*(nombreColonnes - x - 1)) * (y*(nombreLignes - y - 1));
        }
    }
    //float **matseqkplus1 = matseq;

    if(idProc == 0){
        //Imprimer la plaque
        printf("Impression de la plaque initiale (matseq) :\n");
        for(int i=nombreLignes -1 ; i>=0; i--){
            for(int j=0; j<nombreColonnes; j++){
                printf("%.1f  ", matseq[i][j]);
            }
            printf("\n");
        }

        //Sequential
        float part1formula, part2formula, part3formula;

        for(int k=1; k<=nombrePasTemps; k++){
            for(int y=0; y<nombreLignes; y++){
                for(int x=0; x<nombreColonnes; x++){
                    usleep(TEMPS_ATTENTE);
                    if(x==0 || x==(nombreColonnes-1) || y==0 || y==(nombreLignes-1)){
                       matseqkplus1[y][x] = 0;
                    }else{

                        part1formula = (1-(4*(tempsDiscretise/(pow(tailleSubdivision,2))))) * matseq[y][x];
                        part2formula = (tempsDiscretise/(pow(tailleSubdivision,2)));
                        part3formula = 0;

                        if(y != 0){
                            part3formula += matseq[y-1][x];
                            //printf("matseq[%d][%d] : %.1f\n", y-1, x, matseq[y-1][x]);
                        }
                        if(y != nombreLignes - 1){
                            part3formula += matseq[y+1][x];
                            //printf("matseq[%d][%d] : %.1f\n", y+1, x, matseq[y+1][x]);
                        }
                        if(x != 0){
                            part3formula += matseq[y][x-1];
                            //printf("matseq[%d][%d] : %.1f\n", y, x-1, matseq[y][x-1]);
                        }
                        if(x != nombreColonnes - 1){
                            part3formula += matseq[y][x+1];
                            //printf("matseq[%d][%d] : %.1f\n", y, x+1, matseq[y][x+1]);
                        }
                        //printf("1matseq[%d][%d] : %.1f\n", y, x, matseq[y][x]);
                        //printf("1matseqkplus1[%d][%d] : %.1f\n", y, x, matseqkplus1[y][x]);
                        matseqkplus1[y][x] = part1formula + (part2formula * part3formula);
                        //printf("2matseq[%d][%d] : %.1f\n", y, x, matseq[y][x]);
                        //printf("2matseqkplus1[%d][%d] : %.1f\n", y, x, matseqkplus1[y][x]);
                        //printf("SEQ : P1 : %.1f, P2 : %.1f, P3 : %.1f, matseq[%d][%d] : %.1f : total : %.1f\n", part1formula, part2formula, part3formula, y, x, matseq[y][x], part1formula + (part2formula * part3formula));
                    }
                }
            }
            //We copy the value of each data to the current matrix
            for(int y=0; y<nombreLignes; y++){
                for(int x=0; x<nombreColonnes; x++){
                    matseq[y][x] = matseqkplus1[y][x];
                }
            }
            //Imprimer la plaque apres n iterations
            /*printf("Impression de la plaque apres %d iteration (matseqkplus1) :\n", k);
            for(int i=nombreLignes -1 ; i>=0; i--){
                for(int j=0; j<nombreColonnes; j++){
                    printf("%.3f  ", matseqkplus1[i][j]);
                }
                printf("\n");
            }*/
        }


        //Imprimer la plaque apres n iterations
        printf("Impression de la plaque apres %d iteration (matseqkplus1) :\n", nombrePasTemps);
        for(int i=nombreLignes -1 ; i>=0; i--){
            for(int j=0; j<nombreColonnes; j++){
                printf("%.3f  ", matseqkplus1[i][j]);
            }
            printf("\n");
        }
    }































    //PARALLEL
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

    if(nombreColonnes >= nbProc){
        //Si on fait des groupes horizontales
        vy = (float)(nombreColonnes) / (float)(nbProc);
        vx = (float)(nombreLignes);
        //if(((vy - (int)vy) < 0.25 && (vy - (int)vy) * nbProc >= vy && (vy - (int)vy) * nbProc < 2*vy) || (vy - (int)vy) == 0){
           vnbProcUse = nbProc;
        //}else{
           //vnbProcUse = nbProc - 1;
        //}
    }else{
        vy = 1.0f;
        vx = (float)nombreLignes;
        vnbProcUse = nombreColonnes;
    }
    nbCommunications = vnbProcUse - 1;
    vaire = ((int)(vy) - 4) * (int)(vx);
    //printf("vy:%.2f - vx:%.2f - nbConn:%d - aire:%d - vnbProcUse:%d\n", vy, vx, nbCommunications, vaire, vnbProcUse);

    if(nombreLignes >= nbProc){
        //Si on fait des groupes verticales
        hy = (float)(nombreColonnes);
        hx = (float)(nombreLignes) / (float)(nbProc);
        //if(((hx - (int)hx) < 0.25 && (hx - (int)hx) * nbProc >= hx && (hx - (int)hx) * nbProc < 2*hx) || (hx - (int)hx) == 0){
            hnbProcUse = nbProc;
        //}else{
            //hnbProcUse = nbProc - 1;
        //}
    }else{
        hy = (float)nombreColonnes;
        hx = 1.0f;
        hnbProcUse = nombreLignes;
    }
    nbCommunications = hnbProcUse - 1;
    haire = (int)(hy) * ((int)(hx) - 4);
    //printf("hy:%.2f - hx:%.2f - nbConn:%d - aire:%d - hnbProcUse:%d\n", hy, hx, nbCommunications, haire, hnbProcUse);

    //Si on fait des groupes horizontales divisés en deux
    if(nombreColonnes >= (nbProc / 2)){
        h2y = (float)(nombreColonnes) / ((float)(nbProc) / 2);
        h2x = (float)(nombreLignes) / 2;
        h2nbProcUse = (nbProc % 2 == 0 ? nbProc : nbProc - 1);
        nbCommunications = ((h2nbProcUse / 2) * 3) - 2;
        h2aire = abs(((int)(h2y) - 4) * ((int)(h2x) - 2));
        if(((int)(h2y) - 4) < 0 || ((int)(h2x) - 2) < 0){
            h2aire = h2aire * -1;
        }
        //printf("h2y:%.2f - h2x:%.2f - nbConn:%d - aire:%d - h2nbProcUse:%d\n", h2y, h2x, nbCommunications, h2aire, h2nbProcUse);
    }

    //Si on fait des groupes verticales divisés en deux
    if(nombreLignes >= (nbProc / 2)){
        v2y = (float)(nombreColonnes) / 2;
        v2x = (float)(nombreLignes) / ((float)(nbProc) / 2);
        v2nbProcUse = (nbProc % 2 == 0 ? nbProc : nbProc - 1);
        nbCommunications = ((v2nbProcUse / 2) * 3) - 2;
        v2aire = ((int)(v2y) - 2) * ((int)(v2x) - 4);
        if(((int)(v2y) - 2) < 0 || ((int)(v2x) - 4) < 0){
            v2aire = v2aire * - 1;
        }
        //printf("v2y:%.2f - v2x:%.2f - nbConn:%d - aire:%d - v2nbProcUse:%d\n", v2y, v2x, nbCommunications, v2aire, v2nbProcUse);
    }

    //We start the timer
    double timeStart, timeEnd, Texec;
    struct timeval tp;
    gettimeofday (&tp, NULL); // Début du chronomètre
    timeStart = (double) (tp.tv_sec) + (double) (tp.tv_usec) / 1e6;

    //Coupure verticale
    if(idProc < vnbProcUse){ //Eliminate all processors that we don't need
    //if(idProc == 4){ //Eliminate all processors that we don't need

        //Imprimer la plaque
        /*printf("Impression de la plaque initiale (matpar) :\n");
        for(int i=nombreLignes -1 ; i>=0; i--){
            for(int j=0; j<nombreColonnes; j++){
                printf("%.1f  ", matpar[i][j]);
            }
            printf("\n");
        }*/

        ysizestart = (int)((int)vy * (int)(idProc)) + floorf(((float)(idProc) * (float)(vy - (int)vy)) * 10000000) / 10000000;
        ysizeend = ((int)((int)vy + (int)((idProc) * (int)vy)) + floorf(((float)(idProc + 1) * (float)(vy - (int)vy)) * 10000000) / 10000000) - 1;
        xsizestart = 0;
        xsizeend = (int)vx - 1;

        //printf("IdProc : %d --> Largeur : %d\n", idProc, (ysizeend - ysizestart + 1));

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
                mattocalculatekplus1[y][x] = 0;
                //printf("IdProc:%d, Iteration:%d, [%d][%d]:%.3f=[%d][%d]:%.3f\n", idProc, k, y, x,mattocalculate[y][x],y, x, mattocalculatekplus1[y][x]);
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
            //printf("IdProc:%d, Current iteration: %d\n", idProc, k);
            mat_y = 0;
            mat_x = 0;
            MPI_Request leftrequest, rightrequest;
            MPI_Status leftstatus, rightstatus;
            int leftflag = 0, rightflag = 0;
            int betweenindex = 2;
            nbCaseCalculated = 0;
            int leftrecv = 0, rightrecv = 0;
            int firstrowcalculated = 0;
            int lastrowcalculated = 0;
            while(nbCaseCalculated < totalCaseToCalculate){
                //printf("Nbcasecalculated : %d/%d\n", nbCaseCalculated, totalCaseToCalculate);
                //--------------------First iteration-----------------------------------------------------------------
                if(k == 1){ //If first iteration
                    //First, we will calculate the extreme left and right column to send it quickly to the next processor
                    //----------Left column START----------
                    mat_y = 0;
                    mat_x = 0; //Column equals zero
                    for(int y=xsizestart; y<=xsizeend; y++){
                        if(idProc == 0 || y==0 || y==(nombreLignes-1)){
                           mattocalculate[mat_y][mat_x] = 0.0;
                           matleftcolumn[mat_y] = 0.0;
                        }else{
                            part1formula = (1-(4*(tempsDiscretise/(pow(tailleSubdivision,2))))) * matpar[y][ysizestart];
                            part2formula = (tempsDiscretise/(pow(tailleSubdivision,2)));
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
                    /*for(int x=0; x<=(xsizeend - xsizestart); x++){
                        printf("IdProc: %d - Value matleftcolumn[%d] : %.1f\n", idProc, x, matleftcolumn[x]);
                    }*/
                    if(idProc != 0 && nombrePasTemps > 1){
                        printf("idproc:%d send left column on iteration:%d\n", idProc, k);
                        MPI_Isend(matleftcolumn, nombreLignes, MPI_FLOAT, idProc-1, (k+1)*(idProc-1), MPI_COMM_WORLD, &leftrequest);
                    }
                    //----------Left column END----------

                    //----------Right column START----------
                    //On execute la colonne de droite seulement si la matrice à calculer par processeur est plus large que 1 colonne.
                    if(nbColumns > 0){ //0 equals 1 column, 1 equals 2 columns, etc...
                        mat_y = 0;
                        mat_x = nbColumns; //Column equals zero
                        for(int y=xsizestart; y<=xsizeend; y++){
                            //printf("%d, ysizeend:%d\n", iii, ysizeend);
                            if(idProc == nbProc - 1 || y==0 || y==(nombreLignes-1)){
                               //printf("iii : %d, idProc: %d, mattocalculate[%d][%d] : %.1f\n", iii, idProc, mat_y, mat_x, mattocalculate[mat_y][mat_x]);
                               mattocalculate[mat_y][mat_x] = 0.0;
                               matrightcolumn[mat_y] = 0.0;
                            }else{
                                part1formula = (1-(4*(tempsDiscretise/(pow(tailleSubdivision,2))))) * matpar[y][ysizeend];
                                part2formula = (tempsDiscretise/(pow(tailleSubdivision,2)));
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
                        /*for(int x=0; x<=(xsizeend - xsizestart); x++){
                            printf("IdProc: %d - Value matrightcolumn[%d] : %.1f\n", idProc, x, matrightcolumn[x]);
                        }*/
                        if(idProc != (vnbProcUse-1)  && nombrePasTemps > 1){
                            printf("idproc:%d send right column on iteration:%d\n", idProc,k);
                            MPI_Isend(matrightcolumn, nombreLignes, MPI_FLOAT, idProc+1, (k+1)*(idProc+1), MPI_COMM_WORLD, &rightrequest);
                        }
                    }else{
                        if(idProc != (vnbProcUse-1)  && nombrePasTemps > 1){
                            printf("idproc:%d send right column on iteration:%d", idProc,k);
                            MPI_Isend(matleftcolumn, nombreLignes, MPI_FLOAT, idProc+1, (k+1)*(idProc+1), MPI_COMM_WORLD, &rightrequest);
                        }
                    }
                    //----------Right column END----------

                    //----------Between columns START----------
                    if(nbColumns > 1){ //0 equals 1 column, 1 equals 2 columns, etc...
                        mat_y = 0;
                        mat_x = 1; //Second column
                        for(int y=xsizestart; y<=xsizeend; y++){
                            mat_x = 1;
                            for(int x=ysizestart+1; x<=ysizeend-1; x++){
                                if(y == 0 || y == (nombreLignes - 1)){
                                    mattocalculate[mat_y][mat_x] = 0.0;
                                }else{
                                    part1formula = (1-(4*(tempsDiscretise/(pow(tailleSubdivision,2))))) * matpar[y][x];
                                    part2formula = (tempsDiscretise/(pow(tailleSubdivision,2)));
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
                }else{ //All others iterations

                    //MAYBE use MPI_WAIT...
                    if(leftrecv == 0 && idProc != 0){
                        //printf("IdProc: %d, Allocation of the memory Irecv() left column, k=%d\n", idProc, k);
                        MPI_Irecv(matleftcolumnOfProc, nombreLignes, MPI_FLOAT, idProc-1, k*(idProc), MPI_COMM_WORLD, &leftrequest);
                        //MPI_Wait(&leftrequest,&leftstatus);
                        leftrecv = 1;
                    }
                    if(rightrecv == 0 && idProc != vnbProcUse-1){
                        //printf("IdProc: %d, Allocation of the memory Irecv() rightcolumn, k=%d\n", idProc, k);
                        MPI_Irecv(matrightcolumnOfProc, nombreLignes, MPI_FLOAT, idProc+1, k*(idProc), MPI_COMM_WORLD, &rightrequest);
                        //MPI_Wait(&rightrequest,&rightstatus);
                        rightrecv = 1;
                    }

                    if(leftflag == 0 && idProc != 0){
                        //Check for the left processor
                        //MPI_Iprobe( idProc-1, k, MPI_COMM_WORLD, &leftflag, &leftstatus );
                        MPI_Test( &leftrequest, &leftflag, &leftstatus );
                        //printf("IdProc: %d, Test du left...leftflag:%d, nbcasecalculated:%d, k=%d\n", idProc, leftflag, nbCaseCalculated,k);
                    }
                    if(rightflag == 0 && idProc != vnbProcUse-1){
                        //MPI_Iprobe( idProc+1, k, MPI_COMM_WORLD, &rightflag, &rightstatus );
                        MPI_Test( &rightrequest, &rightflag, &rightstatus );
                        //printf("IdProc: %d, Test du right...rightflag:%d, nbcasecalculated:%d, k=%d\n", idProc, rightflag, nbCaseCalculated,k);
                    }

                    //printf("IdProc: %d, leftflag:%d, rightflag:%d, nbcasecalculated:%d, k=%d\n", idProc,leftflag, rightflag, nbCaseCalculated,k);


                    if(nbColumns < 1){ //0 equals 1 column, 1 equals 2 columns, etc...
                        if(idProc == 0 && rightflag == 1){
                            //----------Right column START----------
                            mat_y = 0;
                            mat_x = nbColumns; //Column equals zero
                            for(int y=xsizestart; y<=xsizeend; y++){
                                mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                matrightcolumn[mat_y] = 0.0;
                                mat_y++;
                                nbCaseCalculated++;
                            }
                            if(idProc != (vnbProcUse-1) && k != nombrePasTemps){
                                printf("idproc:%d send right column on iteration:%d\n", idProc,k);
                                MPI_Isend(matrightcolumn, nombreLignes, MPI_FLOAT, idProc+1, (k+1)*(idProc+1), MPI_COMM_WORLD, &rightrequest);
                            }
                            //----------Right column END----------
                        }else if(idProc == vnbProcUse-1 && leftflag == 1){
                            //----------Left column START----------
                            mat_y = 0;
                            mat_x = 0; //Column equals zero
                            for(int y=xsizestart; y<=xsizeend; y++){
                                mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                matleftcolumn[mat_y] = 0.0;
                                mat_y++;
                                nbCaseCalculated++;
                            }
                            if(idProc != 0 && k != nombrePasTemps){
                                printf("idproc:%d send left column on iteration:%d\n", idProc,k);
                                MPI_Isend(matleftcolumn, nombreLignes, MPI_FLOAT, idProc-1, (k+1)*(idProc-1), MPI_COMM_WORLD, &leftrequest);
                            }
                            //----------Left column END----------
                        }else if(leftflag == 1 && rightflag == 1){
                            //----------Between column START----------
                            mat_y = 0;
                            mat_x = 0; //Column equals zero
                            for(int y=xsizestart; y<=xsizeend; y++){
                                if(y==0 || y==(nombreLignes-1)){
                                   mattocalculate[mat_y][mat_x] = 0.0;
                                   matleftcolumn[mat_y] = 0.0;
                                }else{
                                    part1formula = (1-(4*(tempsDiscretise/(pow(tailleSubdivision,2))))) * mattocalculate[y][ysizestart];
                                    part2formula = (tempsDiscretise/(pow(tailleSubdivision,2)));
                                    part3formula = 0;
                                    if(y != 0){part3formula += mattocalculate[y-1][ysizestart];}
                                    if(y != nombreLignes - 1){part3formula += mattocalculate[y+1][ysizestart];}
                                    //if(ysizestart != 0){part3formula += matpar[y][ysizestart-1];}
                                    if(ysizestart != 0){part3formula += matleftcolumnOfProc[(nombreLignes-1) - y];}
                                    //if(ysizeend != nombreColonnes - 1){part3formula += matpar[y][ysizeend+1];}
                                    if(ysizeend != nombreColonnes - 1){part3formula += matrightcolumnOfProc[(nombreLignes-1) - y];}
                                    float tmpresult = part1formula + (part2formula * part3formula);
                                    mattocalculatekplus1[mat_y][mat_x] = tmpresult;
                                    matleftcolumn[mat_y] = tmpresult;
                                }
                                mat_y++;
                                nbCaseCalculated++;
                            }
                            printf("idproc:%d send left and right column on iteration:%d\n", idProc,k);
                            MPI_Isend(matleftcolumn, nombreLignes, MPI_FLOAT, idProc-1, (k+1)*(idProc-1), MPI_COMM_WORLD, &leftrequest);
                            MPI_Isend(matleftcolumn, nombreLignes, MPI_FLOAT, idProc+1, (k+1)*(idProc+1), MPI_COMM_WORLD, &rightrequest);
                            //----------Between column END----------
                        }
                    }else{
                        //If the left processor send the columns, get it with a mpi_recv
                        if(leftflag == 1 && nbCaseCalculated < totalCaseToCalculate){
                            //printf("--IdProc:%d on iteration k=%d, Calculating values for first column leftflag=1 : \n", idProc, k);
                            //printf("IdProc: %d, Before receive...leftflag:%d, nbcasecalculated:%d\n", idProc, leftflag, nbCaseCalculated);
                            //MPI_Recv(matleftcolumnOfProc, nombreLignes, MPI_FLOAT, idProc-1, k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                            //printf("IdProc: %d, After receive...leftflag:%d, nbcasecalculated:%d\n", idProc, leftflag, nbCaseCalculated);
                            //----------Left column START----------
                            mat_y = 0;
                            mat_x = 0; //Column equals zero
                            //printf("\n");
                            for(int y=xsizestart; y<=xsizeend; y++){
                                if(idProc == 0 || y==0 || y==(nombreLignes-1)){
                                   mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                   matleftcolumn[mat_y] = 0.0;
                                   //printf("[%d][%d]:%.1f \n", y, mat_x, mattocalculate[mat_y][mat_x]);
                                }else{
                                    part1formula = (1-(4*(tempsDiscretise/(pow(tailleSubdivision,2))))) * mattocalculate[y][mat_x];
                                    part2formula = (tempsDiscretise/(pow(tailleSubdivision,2)));
                                    part3formula = 0;
                                    if(y != 0){part3formula += mattocalculate[y-1][mat_x];}
                                    if(y != nombreLignes - 1){part3formula += mattocalculate[y+1][mat_x];}
                                    //if(ysizestart != 0){part3formula += matpar[y][ysizestart-1];}
                                    if(ysizestart != 0){part3formula += matleftcolumnOfProc[(nombreLignes-1) - y];}
                                    if(ysizestart != nombreColonnes - 1){part3formula += mattocalculate[y][mat_x+1];}
                                    float tmpresult = part1formula + (part2formula * part3formula);
                                    mattocalculatekplus1[mat_y][mat_x] = tmpresult;
                                    matleftcolumn[mat_y] = tmpresult;

                                }
                                mat_y++;
                                nbCaseCalculated++;
                            }
                            //printf("idproc:%d RECEIVE and CALCULATE FIRST LEFT COLUMN nbCalculated:%d/%d iteration:%d\n", idProc, nbCaseCalculated, totalCaseToCalculate, k);
                            if(idProc != 0 && k != nombrePasTemps){
                                //printf("idproc:%d send left column on iteration:%d\n", idProc,k);
                                MPI_Isend(matleftcolumn, nombreLignes, MPI_FLOAT, idProc-1, (k+1)*(idProc-1), MPI_COMM_WORLD, &leftrequest);
                            }

                        }else if(idProc == 0 && firstrowcalculated != 1){
                            for(int y=xsizestart; y<=xsizeend; y++){
                                mattocalculatekplus1[y][0] = 0.0;
                                nbCaseCalculated++;
                            }
                            firstrowcalculated = 1;
                            leftflag = 1;
                            //printf("Idproc:%d didnt receive and calculate the first row.. nbCalculated:%d/%d iteration:%d\n", idProc, nbCaseCalculated, totalCaseToCalculate, k);
                        }
                        //----------Left column END----------

                        if(rightflag == 1 && nbCaseCalculated < totalCaseToCalculate){
                            //printf("--IdProc:%d on iteration k=%d, Calculating values for last column rightflag=1 : \n", idProc, k);
                            //printf("IdProc: %d, Before receive...rightflag:%d, nbcasecalculated:%d\n", idProc, rightflag, nbCaseCalculated);
                            //MPI_Recv(matrightcolumnOfProc, nombreLignes, MPI_FLOAT, idProc+1, k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                            //printf("IdProc: %d, After receive...rightflag:%d, nbcasecalculated:%d\n", idProc, rightflag, nbCaseCalculated);
                            mat_y = 0;
                            mat_x = nbColumns; //Column equals zero
                            //printf("\n");
                            for(int y=xsizestart; y<=xsizeend; y++){
                                if(idProc == nbProc - 1 || y==0 || y==(nombreLignes-1)){
                                   //printf("iii : %d, idProc: %d, mattocalculate[%d][%d] : %.1f\n", iii, idProc, mat_y, mat_x, mattocalculate[mat_y][mat_x]);
                                   mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                   matrightcolumn[mat_y] = 0.0;
                                   //printf("[%d][%d]:%.1f \n", y, mat_x, mattocalculate[mat_y][mat_x]);
                                }else{
                                    part1formula = (1-(4*(tempsDiscretise/(pow(tailleSubdivision,2))))) * mattocalculate[y][mat_x];
                                    part2formula = (tempsDiscretise/(pow(tailleSubdivision,2)));
                                    part3formula = 0;
                                    if(y != 0){part3formula += mattocalculate[y-1][mat_x];}
                                    if(y != nombreLignes - 1){part3formula += mattocalculate[y+1][mat_x];}
                                    //printf("IdProc:%d, iteration:%d, ysizeend:%d\n", idProc, k, ysizeend);
                                    if(ysizeend != 0){part3formula += mattocalculate[y][mat_x-1];}
                                    //if(ysizeend != nombreColonnes - 1){part3formula += matpar[y][ysizeend+1];}
                                    if(ysizeend != nombreColonnes - 1){part3formula += matrightcolumnOfProc[(nombreLignes-1) - y];}
                                    float tmpresult = part1formula + (part2formula * part3formula);
                                    mattocalculatekplus1[mat_y][mat_x] = tmpresult;
                                    matrightcolumn[mat_y] = tmpresult;
                                }
                                mat_y++;
                                nbCaseCalculated++;
                            }
                            //printf("idproc:%d RECEIVE and CALCULATE LAST RIGHT COLUMN nbCalculated:%d/%d iteration:%d\n", idProc, nbCaseCalculated, totalCaseToCalculate, k);
                            //printf("\n");
                            if(idProc != (vnbProcUse-1) && k != nombrePasTemps){
                                //printf("idproc:%d send right column on iteration:%d\n", idProc,k);
                                MPI_Isend(matrightcolumn, nombreLignes, MPI_FLOAT, idProc+1, (k+1)*(idProc+1), MPI_COMM_WORLD, &rightrequest);
                            }
                        }else if(idProc == vnbProcUse-1 && lastrowcalculated != 1){
                            for(int y=xsizestart; y<=xsizeend; y++){
                                mattocalculatekplus1[y][(ysizeend - ysizestart)] = 0.0;
                                nbCaseCalculated++;
                            }
                            lastrowcalculated = 1;
                            rightflag = 1;
                            //printf("Idproc:%d didnt receive and calculate the last row.. nbCalculated:%d/%d iteration:%d\n", idProc, nbCaseCalculated, totalCaseToCalculate, k);
                        }

                        //----------Between columns START----------
                        if(leftflag == 1 && nbCaseCalculated < totalCaseToCalculate){
                            //printf("-IdProc:%d on iteration k=%d, Calculating values for leftflag TRUE or FIRST IDPROC : ", idProc, k);
                            //We can calculate the first column in between
                            mat_y = 0;
                            mat_x = 1; //Second column
                            for(int y=xsizestart; y<=xsizeend; y++){
                                //printf("Calcul de la case [%d][%d] \n", y, mat_x);
                                if(y == 0 || y == (nombreLignes - 1)){
                                    mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                    //printf("[%d][%d]:%.3f \n", mat_y, mat_x, mattocalculate[mat_y][mat_x]);
                                }else{
                                    part1formula = (1-(4*(tempsDiscretise/(pow(tailleSubdivision,2))))) * mattocalculate[y][mat_x];
                                    part2formula = (tempsDiscretise/(pow(tailleSubdivision,2)));
                                    part3formula = 0;
                                    if(y != 0){part3formula += mattocalculate[y-1][mat_x];}
                                    if(y != nombreLignes - 1){part3formula += mattocalculate[y+1][mat_x];}
                                    if(leftflag==1){
                                        if(idProc != 0){part3formula += matleftcolumnOfProc[(nombreLignes-1) - y];}
                                    }else{
                                        if(mat_x != 0){part3formula += mattocalculate[y][mat_x-1];}
                                    }
                                    if(mat_x != nombreColonnes - 1){part3formula += mattocalculate[y][mat_x+1];}
                                    mattocalculatekplus1[mat_y][mat_x] = part1formula + (part2formula * part3formula);
                                }
                                nbCaseCalculated++;
                                mat_y++;
                            }
                            leftflag = 2;
                            //printf("\n");
                            //printf("idproc:%d CALCULATE SECOND ROW nbCalculated:%d/%d iteration:%d\n", idProc, nbCaseCalculated, totalCaseToCalculate, k);
                        }

                        if(rightflag == 1 && nbCaseCalculated < totalCaseToCalculate){
                            //printf("-IdProc:%d on iteration k=%d, Calculating values for rightflag TRUE or LAST IDPROC : ", idProc, k);
                            //We can calculate the last column in between
                            mat_y = 0;
                            mat_x = (ysizeend - ysizestart-1); //Before last column
                            for(int y=xsizestart; y<=xsizeend; y++){
                                //printf("Calcul de la case [%d][%d] \n", y, mat_x);
                                if(y == 0 || y == (nombreLignes - 1)){
                                    mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                    //printf("[%d][%d]:%.3f \n", mat_y, mat_x, mattocalculate[mat_y][mat_x]);
                                }else{
                                    part1formula = (1-(4*(tempsDiscretise/(pow(tailleSubdivision,2))))) * mattocalculate[y][mat_x];
                                    part2formula = (tempsDiscretise/(pow(tailleSubdivision,2)));
                                    part3formula = 0;
                                    if(y != 0){part3formula += mattocalculate[y-1][mat_x];}
                                    if(y != nombreLignes - 1){part3formula += mattocalculate[y+1][mat_x];}
                                    if(mat_x != 0){part3formula += mattocalculate[y][mat_x-1];}
                                    if(rightflag == 1){
                                        if(idProc != vnbProcUse-1){part3formula += matrightcolumnOfProc[(nombreLignes-1) - y];}
                                    }else{
                                        if(mat_x != nombreColonnes - 1){part3formula += mattocalculate[y][mat_x+1];}
                                    }
                                    //if(mat_x != nombreColonnes - 1){part3formula += mattocalculate[y][mat_x+1];}
                                    mattocalculatekplus1[mat_y][mat_x] = part1formula + (part2formula * part3formula);
                                }
                                nbCaseCalculated++;
                                mat_y++;
                            }
                            rightflag = 2;
                            //printf("\n");
                            //printf("idproc:%d CALCULATE AVANT DERNIERE COLUMN nbCalculated:%d/%d iteration:%d\n", idProc, nbCaseCalculated, totalCaseToCalculate, k);
                        }

                        //if(betweenindex < (ysizeend - ysizestart + 1 - 2)){
                        if(betweenindex <= (ysizeend - ysizestart - 2) && nbCaseCalculated < totalCaseToCalculate){
                            //printf("IdProc:%d on iteration k=%d (%d<=%d), Calculating values : ", idProc, k, betweenindex, (ysizeend - ysizestart - 2));
                            mat_y = 0;
                            mat_x = betweenindex; //Second column
                            for(int y=xsizestart; y<=xsizeend; y++){
                                //printf("[%d][%d] ", y, mat_x);
                                if(y == 0 || y == (nombreLignes - 1)){
                                    mattocalculatekplus1[mat_y][mat_x] = 0.0;
                                }else{
                                    part1formula = (1-(4*(tempsDiscretise/(pow(tailleSubdivision,2))))) * mattocalculate[y][mat_x];
                                    part2formula = (tempsDiscretise/(pow(tailleSubdivision,2)));
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
                            //printf("idproc:%d CALCULATE COLUMN : %d nbCalculated:%d/%d iteration:%d\n", idProc, betweenindex, nbCaseCalculated, totalCaseToCalculate, k);
                            betweenindex++;
                        }
                        //----------Between columns END----------
                    }
                }
                //printf("NbCaseToCalculate : %d, TotalToCalculate %d\n", nbCaseToCalculate, ((ysizeend - ysizestart) + 1) * ((xsizeend - xsizestart) + 1));
            }

            //If it's more than the first iteration
            if(k > 1){
                //We copy the value of each data to the current matrix
                for(int y=0; y<=(ysizeend - ysizestart); y++){
                    for(int x=(xsizeend - xsizestart) ; x>=0; x--){
                        mattocalculate[y][x] = mattocalculatekplus1[y][x];
                        //printf("IdProc:%d, Iteration:%d, [%d][%d]:%.3f=[%d][%d]:%.3f\n", idProc, k, y, x,mattocalculate[y][x],y, x, mattocalculatekplus1[y][x]);
                    }
                }
            }

            /*printf("Impression de la plaque du PROC = %d a l'iteration = %d :\n", idProc, k);
            for(int i=(xsizeend - xsizestart) ; i>=0; i--){
                for(int j=0; j<=(ysizeend - ysizestart); j++){
                    printf("%.1f  ", mattocalculate[i][j]);
                }
                printf("\n");
            }*/
        }
        //printf("idProc:%d iteration's' done \n", idProc);


        /*int ind=0;
        for(int j=0; j<=(ysizeend - ysizestart); j++){
            for(int i=(xsizeend - xsizestart) ; i>=0; i--){
                matparfinal[ind] = mattocalculate[i][j];
                //printf("%.1f [%d][%d]\n", matparfinal[xsizestart + i][ysizestart + j], xsizestart + i, ysizestart + j);
                //printf("IdProc : %d -> [%d][%d]  value [%d][%d] : %.1f\n", idProc, xsizestart + i, ysizestart + j, i, j, mattocalculate[i][j]);
                ind++;
            }
            //printf("\n");
        }*/
        //printf("TEST1 : %.1f [%d]\n", matparfinal[8], 8);






        //------------------------------------SEND TO PROCESSOR 0------------------------------------
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
        //printf("Come on dude : %.1f and ysizestart: %d\n", partmattosend[ind], ysizestart);

        printf("END : IdProc:%d\n", idProc);
        if(idProc != 0){
            //printf("2IdProc:%d\n", idProc);
            //printf("IdPROC: %d - Send data to idproc = 0 FINAL\n", idProc);
            MPI_Isend(partmattosend, totalCaseToCalculate+1, MPI_FLOAT, 0, idProc*(nombrePasTemps+1), MPI_COMM_WORLD, &req);
            //printf("3IdProc:%d ----- sent :)\n", idProc);
        }else{
            //printf("4IdProc:%d\n", idProc);
            for(int j=0; j<=(ysizeend - ysizestart); j++){
                for(int i=(xsizeend - xsizestart) ; i>=0; i--){
                    matpar[i][j] = mattocalculate[i][j];
                }
            }
        }
        //------------------------------------SEND TO PROCESSOR 0------------------------------------
    }

    //MPI_Barrier(MPI_COMM_WORLD);

    if(idProc == 0){
        for(int p=1; p<vnbProcUse; p++){
            //int startcolumn;
            //MPI_Recv(&startcolumn, 1, MPI_INT, p, p, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //printf("IdProc receive : %d, startcolumn : %d\n", p,startcolumn);

            /*
            int arr_length_recv;
            MPI_Status status;

            // Probe du message qui s'en vient du processeur
            MPI_Probe(p, p, MPI_COMM_WORLD, &status);
            // Quand le probe a termine, le status a la grosseur du message, donc on va chercher la grosseur avec MPI_Get_count
            MPI_Get_count(&status, MPI_INT, &arr_length_recv);

            float *partmatrecv;
            partmatrecv = (float *)malloc(arr_length_recv*sizeof(float *));

            //MPI_Recv(matcalculated, nombreLignes * (nbColumnsOfProc+1), MPI_FLOAT, p, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(partmatrecv, arr_length_recv, MPI_FLOAT, p, p, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            */

            int arr_length_recv;
            MPI_Status status;

            // Probe du message qui s'en vient du processeur
            MPI_Probe(p, p*(nombrePasTemps+1), MPI_COMM_WORLD, &status);
            // Quand le probe a termine, le status a la grosseur du message, donc on va chercher la grosseur avec MPI_Get_count
            MPI_Get_count(&status, MPI_INT, &arr_length_recv);

            float *partmatrecv;
            partmatrecv = (float *)malloc(arr_length_recv*sizeof(float *));

            //printf("ready to receive data from p:%d",p);
            //MPI_Recv(matcalculated, nombreLignes * (nbColumnsOfProc+1), MPI_FLOAT, p, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(partmatrecv, arr_length_recv, MPI_FLOAT, p, p*(nombrePasTemps+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE);




            int row = nombreLignes - 1;
            //int col = startcolumn;
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
                printf("%.3f  ", matpar[i][j]);
            }
            printf("\n");
        }
    }


    gettimeofday (&tp, NULL); // Fin du chronomètre
    timeEnd = (double) (tp.tv_sec) + (double) (tp.tv_usec) / 1e6;
    Texec = timeEnd - timeStart; //Temps d'exécution en secondes
    // Fin de l’exemple
    if(idProc == 0){
        printf("Temps d'execution : %f\n", Texec);
    }

    //NE PAS SUPPRIMER ---------------------------------------------------------
    //initialise le mémoire de la matrice
    /*matcolleft = (float **)malloc(((int)vx)*sizeof(float *));
    matcollcenter = (long **)malloc(vy-2*sizeof(long *));
    for(int i=0;i<8;i++){
        mat[i] = (long *)malloc(8*sizeof(long));
    }
    matcollright = (float **)malloc(((int)vx)*sizeof(float *));*/

    //float **matcolleftkmoins1 = matcolleft;
    //float **matcollrightkmoins1 = matcollright;
    //NE PAS SUPPRIMER ---------------------------------------------------------


    MPI_Finalize();
    return 0;
}
