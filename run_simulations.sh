for seed in {41801..41810..1}
do
  for idf in Id2
  do
    for k0 in 3 4
    do
      for ke in 1 2
      do
        for src in 1 3
        do
	        for infl in 0.2 1 3
	        do
	          for spc in 0 50
	          do
	            for mod in baseT sparseT
	            do
		            echo "Rscript --no-restore 3_simulation_study.R $k0 $ke $src $infl $spc $seed $mod $idf &> progress/sim_$mod\_$idf\_K0$k0\_Ke$ke\_$src\present_infl$infl\_$spc\sparse_seed$seed\.txt" || continue
	            done
              done
          	done
		done
	  done
	done
  done
done | parallel --jobs 10

wait

# jobs can't be more than 1/3 of total cores available because each runs 3 chains
