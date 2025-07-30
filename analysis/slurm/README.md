## Folder Contents

Scripts in this file were used to analyse the whole set of sls-translatable paths using slurm.

### analyse_output.py

Main script to analyse the results from obj_tester.sh. Does a few statistics and plots some figures. 

### check_drt_output.py 

Script which takes DrTrafo output and checks if a given aCFP is satisified. 

### make_que.sh 

Creates a txt file containing numbers 1-n for obj_tester.sh

### obj_tester.sh

Main bash script which gets send to the cluster for analysis. Takes an index from the que file and then processes the corresponding aCFP and simulates it on the domain level and generates a nucleotide sequence and checks wether or not the design was succesfull. Repeats the process until it found a sucessfull design or it reaches the try limit. 

### slurm_obj_tester.py

Corresponding python script which is used by obj_tester.sh to do the whole analysis. 

