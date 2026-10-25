# Build T0-restricted sub-instances (rows <= T0) for selected D with the line's own tools (filt_T.py, mkinst.py),
# into audit/instT<T0>/.  Usage: sh subinst.sh T0 D1,D2,...
T=$1; DS=$2; PY=${PY:-python3}
mkdir -p instT$T
$PY ../filt_T.py $T ../c3_1e10.txt ../chm200_pool.txt > recT$T.txt
rm -f instT$T/index.tsv
(cd .. && $PY mkinst.py audit/instT$T $DS audit/recT$T.txt) && wc -l instT$T/index.tsv
