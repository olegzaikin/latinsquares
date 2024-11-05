n=$1

for f in ./dls_n${n}_first_row_ascending_*_inverse*.cnf
do
 echo "Processing $f"
 base_name=$(basename -- "$f" .c)
 #echo $base_cnfname
 clasp --models 0 --quiet --enum-mode=bt --configuration=crafty $f
done
