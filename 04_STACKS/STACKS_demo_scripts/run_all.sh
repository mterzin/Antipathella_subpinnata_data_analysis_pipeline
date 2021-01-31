#!/bin/bash

cd $(dirname $0)

# e.g.
# cd /path/to/demo_scripts
# qsub -N top.xx -joe -o $PWD -q catchenlab -l walltime=168:00:00 <<< "cd $PWD; ./run_all.sh"

export PATH="$STACKS147/bin":"$PATH"

scripts='
0.prepare_env.sh
1.survey_lanes.sh
2.clean_lanes.sh
3.genome_db.sh
4.tests_denovo.sh
5.tests_ref.sh
6.full_denovo.sh
7.full_ref.sh
8.genepop_analyses.sh
'

for script in $scripts ;do
	echo
	echo "$script..."
	/usr/bin/time ./$script
done
