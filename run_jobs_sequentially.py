import os
import datetime
import glob
import subprocess
import time
from argparse import ArgumentParser
from contextlib import contextmanager

# run like: python run_jobs_sequentially.py --device opencl:0.0
parser = ArgumentParser()
parser.add_argument("--device", help="set opencl device")
args = parser.parse_args()

@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)

pwd = os.getcwd()
jobs_waiting_file = "jobs_waiting.txt"
jobs_running_file = "jobs_running_%s.txt"%(args.device)
jobs_finished_file = "jobs_finished.txt"
while True:
    with open(jobs_waiting_file, "r") as f:
        directory = f.readline().strip('\n')
        os.system('echo "$(tail -n +2 %s)" > %s'%(jobs_waiting_file,jobs_waiting_file))
    if not directory:
        print("INFO +++ No more jobs left to run ... DONE")
        break

    # do some stuff
    if os.path.isdir(directory):
        print('INFO +++ starting simulation in %s'%(directory))
        with open(jobs_running_file, "w") as output:
            output.write(directory + '\n')

        # do the actual simulation
        with cd(directory):
            os.system("bash run_jobfiles.sh --device %s"%(args.device))
            #time.sleep(10)

        with open(jobs_running_file, 'w') as f:
            f.write('')
   
    else:
        print('WARNING +++ directory does not exist: "%s"'%(directory))    

    time.sleep(10)
       
    with open(jobs_finished_file, "a") as output:
        output.write('%s    %s\n'%(str(datetime.datetime.now())[:-7], directory))
