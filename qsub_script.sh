foreach joe (`ls discovery2014/ALL_SAMPLES/ | grep Sample`)
  echo $joe
  qsub -q baliga script.py --sample $joe
end

