##################################################################
# Fichier “Makefile”
# make all: mettre a jour le tout
# make clean: enlever les version executable
################ #################################################
CFLAGS = -g -std=c99 -D _BSD_SOURCE
COMP = mpicc $(CFLAGS)
ALL = tp1_par

all: $(ALL)
clean:
	rm $(ALL)
tp1_par: tp1_par.c
	$(COMP) tp1_par.c -o tp1_par
