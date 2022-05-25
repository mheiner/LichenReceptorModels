for seed in {41901..41910..1}
do
  for idf in Id2
  do
    for k0 in 4 5
    do
      for ke in 0 1 2
      do
        for infl in 0.5 1 2
        do
          for mod in baseT sparseT
          do
            echo "Rscript --no-restore  1_run_models.R $mod $idf $k0 $ke $infl $seed &> progress/sim_$mod\_$idf\_K0$k0\_Ke$ke\_infl$infl\_seed$seed\.txt" || continue
          done
		    done
		  done
		done
	done
done | parallel --jobs 10

wait

# jobs can't be more than 1/3 of total cores available because each runs 3 chains