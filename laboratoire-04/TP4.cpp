#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <CL/cl.h>
#include <ctime>
#include <cmath>

//trouver sur internet pour la tailled d'un fichier
#define MEM_SIZE (128)
#define MAX_SOURCE_SIZE (0x100000)

using namespace std;

int main(int argc, char* argv[])
{
	//N
	float nombreLignes;
	//M
	float nombreColonnes;
	float nombrePasTemps;
	float tempsDiscretise;
	float tailleSubdivision;

	//exit si nous avons pas 6 parametres
	if(argc != 6){
		cout << "Il vous manque des paramètres, les parametres par defaut seront misent";
		nombreLignes = 10.0f;
		nombreColonnes = 10.0f;
		nombrePasTemps = 5000.0f;
		tempsDiscretise = 0.005f;
		tailleSubdivision = 0.1f;
	}
	else{
		nombreLignes = atof(argv[1]);
		nombreColonnes = atof(argv[2]);
		nombrePasTemps = atof(argv[3]);
		tempsDiscretise = atof(argv[4]);
		tailleSubdivision = atof(argv[5]);
	}



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
            matseq[y][x] = (x*(nombreColonnes - x - 1)) * (y*(nombreLignes - y - 1));
        }
    }

    //intialisation du timer apres l'initialisation
	clock_t wStartSeq, wEndSeq;
	wStartSeq = clock();

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
    cout << "Impression de la plaque apres " << nombrePasTemps << "iteration (matseqkplus1) :\n";
    for(int i=nombreLignes -1 ; i>=0; i--){
        for(int j=0; j<nombreColonnes; j++){
            printf("%.5f  ", matseqkplus1[i][j]);
        }
        printf("\n");
    }

    // Fin du chronomètre
	wEndSeq = clock();
	double wTimeElapsed = (double)(wEndSeq - wStartSeq) / CLOCKS_PER_SEC * 1000.0;

	cout << "Temps ecoule sequentiel: \t" << wTimeElapsed << "\n";
	// Fin de l’exemple


	//partie parallele (OpenCL)

	float *matriceInitiale;
	float *matriceFinale;

	clock_t wStartPar, wEndPar;

	float longueurVecteur = nombreLignes*nombreColonnes;
	size_t bytes = longueurVecteur*sizeof(float);

	matriceInitiale = (float*)malloc(bytes);
	matriceFinale = (float*)malloc(bytes);

	cl_mem d_matriceInitiale;
	cl_mem d_matriceFinale;

	//cree une matrice 1xLONGUEUR avec la logique d'une matrice MxN
	//matrice initiale

	int secondIndex = 0;
	bool FirstCol = true;
	for(int colNum = 0; colNum < nombreColonnes; colNum++){
	  if(FirstCol){
		for(int rowNumber = 0; rowNumber < nombreLignes; rowNumber++){
		  matriceInitiale[rowNumber] = abs((colNum * (nombreColonnes - colNum - 1.0)) * ( rowNumber* (nombreLignes - (rowNumber) - 1.0)));
		  if(rowNumber == nombreLignes-1){
			secondIndex = rowNumber;
			FirstCol = false;
		  }
		}
	  }
	  else{
		for(int rowNumber = 1; rowNumber <= nombreLignes; rowNumber++){
		  matriceInitiale[secondIndex+rowNumber] = abs((colNum * (nombreColonnes - colNum - 1.0)) * ((rowNumber-1) * (nombreLignes - (rowNumber-1) - 1.0)));
		  if(rowNumber == nombreLignes){
			secondIndex = secondIndex + rowNumber;
		  }
		}
	  }
	}

	cout << "Affichage de la matrice initale parallele : \n";

	for(int i = 0; i < nombreLignes; i++){
		for(int cols = 0; cols < nombreColonnes ;cols++){
			cout << matriceInitiale[(int)(cols*nombreLignes)+i] << "\t";
		}
		cout << "\n";
	}

	size_t globalThreadMap, localThreadMap;
	localThreadMap = longueurVecteur;
	globalThreadMap = (ceil(longueurVecteur/(float)localThreadMap)*localThreadMap)*(int)nombrePasTemps;


	cl_device_id *device_id = NULL;
	cl_context context = NULL;
	cl_command_queue command_queue = NULL;
	cl_program program = NULL;
	cl_kernel kernel = NULL;
	cl_platform_id platform_id = NULL;
	cl_uint ret_num_devices;
	cl_uint ret_num_platforms;
	cl_int status;

	FILE *fp;
	char fileName[] = "./TP4.cl";
	char *source_str;
	size_t source_size;

	//load le fichier source .cl

	fp = fopen(fileName, "r");
	if (!fp) {
		fprintf(stderr, "Failed to load kernel.\n");
	exit(1);
	}
	source_str = (char*)malloc(MAX_SOURCE_SIZE);
	source_size = fread(source_str, 1, MAX_SOURCE_SIZE, fp);
	fclose(fp);

	status = clGetPlatformIDs(1, NULL, &ret_num_platforms);

	if(status != CL_SUCCESS){
		cout << "Erreur: Platformes non-trouvees!" << endl;
		return 56;
	}

	cl_platform_id* platforms = (cl_platform_id* )malloc(ret_num_platforms* sizeof(cl_platform_id));
	status = clGetPlatformIDs(ret_num_platforms, platforms, NULL);

	int index = 0;
	bool gpuSet = false;

	//Assigne le GPU pour faire le traitement, sinon exit.

	while(index < ret_num_platforms && gpuSet == false){
		platform_id = platforms[index];
		status = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, 0, NULL, &ret_num_devices);
		if(status == CL_SUCCESS){
			gpuSet = true;
			free(platforms);
		}
		index++;
	}

	if(status != CL_SUCCESS){
		free(platforms);
		cout << "Pas de carte graphique disponible." << endl;
		cout << "Exit." << endl;
		return 57;
	}
	else
	{
		device_id = (cl_device_id*)malloc(ret_num_devices * sizeof(cl_device_id));
		status = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, ret_num_devices, device_id, NULL);
	}


	context = clCreateContext(NULL, 1, device_id, NULL, NULL, &status);

	command_queue = clCreateCommandQueue(context, *device_id, 0, &status);

	program = clCreateProgramWithSource(context, 1, (const char **)&source_str,
	(const size_t *)&source_size, &status);

	status = clBuildProgram(program, 1, device_id, NULL, NULL, NULL);

	kernel = clCreateKernel(program, "calculateMatrix", &status);

	d_matriceInitiale = clCreateBuffer(context, CL_MEM_READ_WRITE,bytes,NULL,&status);
	d_matriceFinale = clCreateBuffer(context,CL_MEM_READ_WRITE,bytes,NULL,&status);

	

	//parametres OpenCL
	status = clSetKernelArg(kernel, 0, sizeof(cl_mem), &d_matriceInitiale);
	status = clSetKernelArg(kernel, 1, sizeof(cl_mem), &d_matriceFinale);
	status = clSetKernelArg(kernel, 2, sizeof(float), &nombreLignes);
	status = clSetKernelArg(kernel, 3, sizeof(float), &nombreColonnes);
	status = clSetKernelArg(kernel, 4, sizeof(float), &nombrePasTemps);
	status = clSetKernelArg(kernel, 5, sizeof(float), &tempsDiscretise);
	status = clSetKernelArg(kernel, 6, sizeof(float), &tailleSubdivision);

	status = clEnqueueWriteBuffer(command_queue,d_matriceInitiale,CL_TRUE,0,
		bytes ,matriceInitiale,0,NULL,NULL);


	wStartPar = clock();
	bool toRead = false;
	for(int k = 1; k <= (int)nombrePasTemps; k++){

		status = clEnqueueNDRangeKernel(command_queue,kernel,1,NULL,
			&globalThreadMap,NULL,0,NULL,NULL);

		if(k%2 > 0){
			status = clEnqueueReadBuffer(command_queue, d_matriceFinale, CL_TRUE, 0,bytes,matriceFinale, 0, NULL, NULL);
			status = clEnqueueWriteBuffer(command_queue,d_matriceInitiale,CL_TRUE,0,bytes,matriceFinale,0,NULL,NULL);
		}else{
			status = clEnqueueReadBuffer(command_queue, d_matriceFinale, CL_TRUE, 0,bytes,matriceInitiale, 0, NULL, NULL);
			status = clEnqueueWriteBuffer(command_queue,d_matriceInitiale,CL_TRUE,0,bytes,matriceInitiale,0,NULL,NULL);
		}

		if(k == (int)nombrePasTemps){
			if(k%2 > 0){
				toRead = false;
			}
			else{
				toRead = true;
			}
		}
	}

	status = clFinish(command_queue);
	wEndPar = clock();

	cout << wStartPar << "\n";
	cout << wEndPar << "\n";



	double wTimeElapsedParallele = (double)(wEndPar - wStartPar) / CLOCKS_PER_SEC * 1000.0;

	cout << "Temps ecoule parallel: \t" << wTimeElapsedParallele << "\n";
	cout << "Acceleration par rapport a l'execution sequentielle:\t\n" << wTimeElapsed/wTimeElapsedParallele << "\n";


	cout << "Affichage de la matrice finale\n";
	for(int i = 0; i < nombreLignes; i++){
		for(int cols = 0; cols < nombreColonnes ;cols++){
			if(toRead == false){
				cout << matriceFinale[(int)(cols*nombreLignes)+i] << "\t";
			}
			else{
				cout << matriceInitiale[(int)(cols*nombreLignes)+i] << "\t";
			}
		}
		cout << "\n";
	}


	//Libere la memoire utilisee
	status = clFlush(command_queue);
	status = clReleaseKernel(kernel);
	status = clReleaseProgram(program);
	status = clReleaseMemObject(d_matriceInitiale);
	status = clReleaseMemObject(d_matriceFinale);
	status = clReleaseCommandQueue(command_queue);
	status = clReleaseContext(context);

	free(source_str);
	free(matriceFinale);
	free(matriceInitiale);
	free(matseq);
	free(matseqkplus1);


	int waitForExit = 0;
	cout << "Veuillez entrez un chiffre pour quitter" << "\n";
	cin >> waitForExit;
	return 0;
}
