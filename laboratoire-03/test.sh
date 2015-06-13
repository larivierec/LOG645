chmod a+x run.sh
make
./run.sh $1 $2 $3 $4 $5 $6
make clean
# Le repertoire en devrait pas contenir d'executable
ls
