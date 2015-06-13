__kernel void calculateMatrix(__global float* matInit, __global float *matFinale, const float nombreLignes,
							  const float nombreColonnes, const float nombrePasTemps,
							  const float tempsDiscretise, const float tailleSubdivision){

	int globalID = get_global_id(0);
	int nombreLines = (int) nombreLignes;
	int nombreCol = (int) nombreColonnes;

	float formuleP1 = 0;
	float formuleP2 = 0;
	float formuleP3 = 0;

	if(globalID < nombreLines * nombreCol){

		if((globalID < nombreLignes) || (globalID % (int)nombreLignes) == 0 || (globalID >= ((nombreLines * nombreCol) - nombreLines)) || (((globalID-nombreLines)+1) % (int)nombreLignes) == 0){
			matFinale[globalID] = 0.0;
		}
		else{
			formuleP1 = (1-(4*tempsDiscretise/tailleSubdivision*tailleSubdivision)) * matInit[globalID];
			formuleP2 = (tempsDiscretise/tailleSubdivision*tailleSubdivision);
			formuleP3 = matInit[globalID-nombreLines] + matInit[globalID+nombreLines] + matInit[globalID-1] + matInit[globalID+1];
			matFinale[globalID] = formuleP1 + (formuleP2 * formuleP3);
		}
	}
}
