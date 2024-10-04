for i in $(seq 0 31);
do
root -b -l <<EOF
.L TestEfficiency.cxx+
.L RawCorrelation.cxx+
TestEfficiency($i)
RawCorrelation($i)
.q
EOF
done