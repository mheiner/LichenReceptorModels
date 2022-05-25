for seed in {42701..42704..1}
do
  for idf in Id2
  do
    for k0 in 5
    do
      for ke in 2
      do
        for infl in 1
        do
          for mod in baseT sparseT # base sparse nullT null
          do
            echo "Rscript --no_restore  1_run_models.R $mod $idf $k0 $ke $infl $seed &> progress/sim_$mod\_$idf\_K0$k0\_Ke$ke\_infl$infl\_seed$seed\.txt" || continue
          done
		    done
		  done
		done
	done
done | parallel --jobs 10

wait

# jobs can't be more than 1/4 of total cores available because each runs 3 or 4 chains