import os
import sys
import subprocess
from select import select

# Start two subprocesses
procs = [
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_1080.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_1110.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_1140.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_1170.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_1200.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_1230.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_1260.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_1290.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_1320.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_1350.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_1380.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_1410.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_1440.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_1470.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_1500.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_1530.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_1560.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_1590.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
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