for i in $(seq 1 10);
do
root -b -l <<EOF
.L TestEfficiency.cxx+
.L RawCorrelation.cxx+
TestEfficiency($i)
RawCorrelation($i)
.q
EOF
done