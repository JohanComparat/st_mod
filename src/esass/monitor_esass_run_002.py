import os
import sys
import subprocess
from select import select

# Start two subprocesses
procs = [
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_SimBKG_CLUseed002_processing_0720.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_SimBKG_CLUseed002_processing_0750.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_SimBKG_CLUseed002_processing_0780.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_SimBKG_CLUseed002_processing_0810.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_SimBKG_CLUseed002_processing_0840.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_SimBKG_CLUseed002_processing_0870.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_SimBKG_CLUseed002_processing_0900.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_SimBKG_CLUseed002_processing_0930.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_SimBKG_CLUseed002_processing_0960.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_SimBKG_CLUseed002_processing_0990.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_SimBKG_CLUseed002_processing_1020.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_SimBKG_CLUseed002_processing_1050.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_SimBKG_CLUseed002_processing_1080.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_SimBKG_CLUseed002_processing_1110.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_SimBKG_CLUseed002_processing_1140.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_SimBKG_CLUseed002_processing_1170.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_SimBKG_CLUseed002_processing_1200.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_SimBKG_CLUseed002_processing_1230.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_SimBKG_CLUseed002_processing_1260.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_SimBKG_CLUseed002_processing_1290.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_SimBKG_CLUseed002_processing_1320.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_SimBKG_CLUseed002_processing_1350.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_SimBKG_CLUseed002_processing_1380.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_SimBKG_CLUseed002_processing_1410.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_SimBKG_CLUseed002_processing_1440.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_SimBKG_CLUseed002_processing_1470.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    # subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    # subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    # subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    # subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    # subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    # subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    # subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    # subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    # subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    # subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    # subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    # subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    # subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    # subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    # subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
]
# Continuously monitor output
while procs:
    # Remove processes that have finished
    for p in procs[:]:
        if p.poll() is not None:
            # Drain remaining output
            for line in p.stdout:
                print(f"[PID {p.pid}]: {line}", end='')
            p.stdout.close()
            procs.remove(p)

    if not procs:
        break

    # Wait until any process has output
    rlist, _, _ = select([p.stdout for p in procs], [], [], 0.1)
    for pipe in rlist:
        line = pipe.readline()
        if line:
            # Identify which process produced this
            pid = next(p.pid for p in procs if p.stdout is pipe)
            print(f"[PID {pid}]: {line}", end='')
