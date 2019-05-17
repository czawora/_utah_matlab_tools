
f=construct_spikeInfoMS
mcc -R -singleCompThread -I ../  -mv ../$f.m \
&& num_lines=$(cat "run_"$f".sh" | wc -l) \
&& num_comment=$(grep -x "^#.*$" "run_"$f".sh" | wc -l) \
&& f0="run_"$f".sh" \
&& f2="run_"$f"_swarm.sh" \
&& (head -n $num_comment $f0 > $f2;
cat ../matlab_header.txt >> $f2;
tail -n $((num_lines-num_comment)) $f0 >> $f2;) \
&& echo "compiled!"

chmod 777 *
