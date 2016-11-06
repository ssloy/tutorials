>analyzer.log
for i in `seq 1 191`; do cat logic-1-$i >> analyzer.log ; done
rm logic-1-*
rm version
rm metadata
