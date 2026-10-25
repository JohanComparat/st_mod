import os
import sys
import subprocess
from select import select

# Start two subprocesses
procs = [
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed009_SimBKG_CLUseed081_processing_0360.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed009_SimBKG_CLUseed081_processing_0540.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed009_SimBKG_CLUseed081_processing_0720.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed082_processing_0000.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed082_processing_0180.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed082_processing_0360.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed082_processing_0540.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed002_SimBKG_CLUseed083_processing_0000.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed002_SimBKG_CLUseed083_processing_0180.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed002_SimBKG_CLUseed083_processing_0360.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed002_SimBKG_CLUseed083_processing_0540.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed002_SimBKG_CLUseed083_processing_0720.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed002_SimBKG_CLUseed083_processing_0900.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed002_SimBKG_CLUseed083_processing_1080.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed002_SimBKG_CLUseed083_processing_1260.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed003_SimBKG_CLUseed084_processing_0000.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed003_SimBKG_CLUseed084_processing_0180.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed003_SimBKG_CLUseed084_processing_0360.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed003_SimBKG_CLUseed084_processing_0540.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed003_SimBKG_CLUseed084_processing_0720.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed003_SimBKG_CLUseed084_processing_0900.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed003_SimBKG_CLUseed084_processing_1080.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed003_SimBKG_CLUseed084_processing_1260.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed085_processing_0000.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed085_processing_0180.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed085_processing_0360.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed085_processing_0540.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed085_processing_0720.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed085_processing_0900.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed085_processing_1080.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed085_processing_1260.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed005_SimBKG_CLUseed086_processing_0000.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
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

