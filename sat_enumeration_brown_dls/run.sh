n=$1
./sat_enumeration_brown_dls $n
find . -name "brown_dls_n${n}_cover*" -exec cat {} \; > brown_dls_n${n}
sort -u brown_dls_n${n} > brown_dls_n${n}_nodupl
wc -l brown_dls_n${n}_nodupl
rm brown_dls_n${n}
rm brown_dls_n${n}_cover*
