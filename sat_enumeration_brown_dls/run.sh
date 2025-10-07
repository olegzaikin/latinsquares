n=$1
cms=$2
./sat_enumeration_brown_dls $n $cms > out_enum_brown_dls_n${n} -cpunum=12
# Find all horizontal symmetric DLS and copy them to one file:
find . -name "brown_dls_n${n}_cover*_horiz" -exec cat {} \; > brown_dls_n${n}_horiz
rm -f brown_dls_n${n}_cover*_horiz
# Find all vertically symmetric DLS and copy them to one file:
find . -name "brown_dls_n${n}_cover*_vertic" -exec cat {} \; > brown_dls_n${n}_vertic
rm -f brown_dls_n${n}_cover*_vertic
# Sort all horizontal symmetric DLS:
sort -u brown_dls_n${n}_horiz > brown_dls_n${n}_horiz_nodupl
rm -f brown_dls_n${n}_horiz
wc -l brown_dls_n${n}_horiz_nodupl
# Sort all vertical symmetric DLS:
sort -u brown_dls_n${n}_vertic > brown_dls_n${n}_vertic_nodupl
rm -f brown_dls_n${n}_vertic
wc -l brown_dls_n${n}_vertic_nodupl
# Sort all DLS:
cat brown_dls_n${n}_horiz_nodupl brown_dls_n${n}_vertic_nodupl > brown_dls_n${n}
sort -u brown_dls_n${n} > brown_dls_n${n}_nodupl
rm -f brown_dls_n${n}
wc -l brown_dls_n${n}_nodupl
# Find intersections between horizontal and vertical:
comm -12 brown_dls_n${n}_horiz_nodupl brown_dls_n${n}_vertic_nodupl > brown_dls_n${n}_horiz-and-vertic_nodupl
wc -l brown_dls_n${n}_horiz-and-vertic_nodupl
