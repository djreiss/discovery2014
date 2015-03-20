#!/bin/csh -f

## email: abort (a), beginning (b), end (e)
##$ -m abe
#$ -m a
#$ -M dreiss@systemsbiology.org
#$ -P Baliga
#$ -cwd
##$ -l h_rt=10000:00:00
##$ -l cores=8
##$ -t 701-999
##### Tell it you're using 20 threads, so we don't fill the whole node with 40 threads (20 cores!)
#$ -pe serial 8
##$ -o output/$TASK_ID.log 
#$ -j y

#uname -a
#pwd

setenv PATH ".:/tools/bin:${PATH}"
rehash
echo "HERE"
echo ${1}
pwd

#cd /proj/omics4tb/dreiss/discovery2014/TEST
/users/dreiss/miniconda/bin/python starproc.py --sample ${1} >>&${1}.out


