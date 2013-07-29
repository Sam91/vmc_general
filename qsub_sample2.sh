# This is a sample PBS script. It will request 1 processor on 1 node
# for 4 hours.
#   
#   Request 1 processors on 1 node 
#   
# PPP -l nodes=1:ppn=1
#
#   Request 4 hours of walltime
#
#PBS -l walltime=28:00:00
#
#   Request 1 gigabyte of memory per process
#
#PBS -l pmem=400mb
#
#   Request that regular output and terminal output go to the same file
#
#PBS -j oe
#
#   The following is the body of the script. By default,
#   PBS scripts execute in your home directory, not the
#   directory from which they were submitted. The following
#   line places you in the directory from which the job
#   was submitted.
#
cd $PBS_O_WORKDIR
cd /home/sbieri/vmc_general/q_scripts
#
#   Now we want to run the program "hello".  "hello" is in
#   the directory that this script is being submitted from,
#   $PBS_O_WORKDIR.
#

export LD_LIBRARY_PATH=/opt/openmpi/lib:/home/sbieri/vmc_general/lib

echo " "
echo " "
echo "Job started on `hostname` at `date`"
echo " "

#./main_u1hyb 12 48 100 20 50 300 40 20
/home/sbieri/vmc_general/prod/executable AA AB AC AD AE AF AG AH AI AJ AK AL AM AN AO AP

echo " "
echo " "
echo "Job Ended at `date`"
echo " "

