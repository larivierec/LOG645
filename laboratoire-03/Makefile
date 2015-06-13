##################################################################
# Fichier “Makefile”
# make all: mettre a jour le tout
# make clean: enlever les version executable
################ #################################################
CFLAGS = -g -std=c99 -D _BSD_SOURCE
COMP = mpicc $(CFLAGS)
ALL = tp3

all: $(ALL)
clean:
	rm $(ALL)
tp3: tp3.c
	$(COMP) tp3.c -o tp3
