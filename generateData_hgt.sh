#!/bin/bash
set -e
set -x

cwd=`pwd`

# rep is number species trees, genes number of locus trees per species trees, species number of species, b is distribution of birth rates, t is distribution of maximum tree lengths, tr is distribution of transfer rates 

rep=10
genes=500
species=10000
for b in 0.0000005; do
	for t in 100000000; do
		#for tr in 0 0.0000002 0.0000005 0.00000002 0.000000005 0.000000002; do 
		for tr in norm-lognormal; do 
		#for tr in  0; do 
			#simphy -rs $rep -rl U:$genes,$genes -rg 1 -st U:$t,$t -si U:1,1 -sl U:$species,$species -sb U:$b,$b -cp U:200000,200000  -hs L:1.5,1 -hl L:1.2,1 -hg l:1.4,1 -cu E:10000000 -so U:1,1 -od 1 -or 0 -v 3  -cs 293745 -o model.$species.$t.$b.$tr -lt U:$tr,$tr -lk 1| tee log.$species.$t.$b.$tr;
            simphy -rs "$rep" -rl f:"$genes" -rg 1 -sb f:"$b" -sd f:`python -c "print($b/1.2)"` -st f:$t -sl f:"$species"  -si f:1 -sp f:500000 -su f:4e-08 -hh f:1 -hs ln:1.5,1 -hl ln:1.3692114,0.6931472 -hg ln:1.5,1 -cs 14907 -gt n:-18,0.4 -lt ln:gt,0.75 -lk 1 -lb f:0 -ld f:0 -v 3 -o model.$species.$t.$b.$tr -ot 0 -op 1 -od 1 > log.txt
			for r in `ls -d model.$species.$t.$b.$tr/*/`; do 
			        sed -i -e "s/_0_0//g" $r/g_trees*.trees; 
			done
			perl ./post_stidsim.pl `pwd`/model.$species.$t.$b.$tr 1 # You can comment this line in your initial tests
			for r in `ls -d model.$species.$t.$b.$tr/*/`; do
				cat $r/g_trees*.trees > $r/truegenetrees;
				rm  $r/g_trees*.trees;
				# ( cd $r && ~/bin/indelible )
			done
			echo look at model.$species.$t.$b.$tr/*/l_trees.trees for figuring out the number of transfer events
		done
	done
done

cd $cwd
