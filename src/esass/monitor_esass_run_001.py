import os
import sys
import subprocess
from select import select

# Start two subprocesses
procs = [
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_AGNseed001_SimBKG_processing_0000.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_AGNseed001_SimBKG_processing_0050.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_AGNseed001_SimBKG_processing_0100.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_AGNseed001_SimBKG_processing_0150.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_AGNseed001_SimBKG_processing_0200.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_AGNseed001_SimBKG_processing_0250.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_AGNseed001_SimBKG_processing_0300.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_AGNseed001_SimBKG_processing_0350.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_AGNseed001_SimBKG_processing_0400.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_AGNseed001_SimBKG_processing_0450.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_AGNseed001_SimBKG_processing_0500.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_AGNseed001_SimBKG_processing_0550.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_AGNseed001_SimBKG_processing_0600.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_0000.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_0050.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_0100.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_0150.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_0200.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_0250.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_0300.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_0350.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_0400.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_0450.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_0500.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_0550.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_0600.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_0650.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_0700.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_0750.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_0800.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_SimBKG_CLUseed001_processing_0000.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_SimBKG_CLUseed001_processing_0050.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_SimBKG_CLUseed001_processing_0100.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_SimBKG_CLUseed001_processing_0150.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_SimBKG_CLUseed001_processing_0200.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e4_merge_SimBKG_CLUseed001_processing_0250.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
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
