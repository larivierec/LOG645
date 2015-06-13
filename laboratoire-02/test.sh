chmod a+x run.sh
make
./run.sh $1 $2 $3
make clean
# Le repertoire en devrait pas contenir d'executable
ls
