foreach joe (`ls ALL_SAMPLES/ | grep Sample`)
  echo $joe
  qsub -q baliga run_starproc.sh "${joe}"
end

