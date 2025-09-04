import os
import sys
import subprocess
from select import select




# Start two subprocesses
procs = [
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed003_SimBKG_CLUseed003_processing_0000.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed004_SimBKG_CLUseed004_processing_0000.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed004_SimBKG_CLUseed004_processing_0030.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed004_SimBKG_CLUseed004_processing_0060.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed004_SimBKG_CLUseed004_processing_0090.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed004_SimBKG_CLUseed004_processing_0120.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed004_SimBKG_CLUseed004_processing_0150.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed004_SimBKG_CLUseed004_processing_0180.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed005_SimBKG_CLUseed005_processing_0000.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed005_SimBKG_CLUseed005_processing_0030.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed005_SimBKG_CLUseed005_processing_0060.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed005_SimBKG_CLUseed005_processing_0090.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed005_SimBKG_CLUseed005_processing_0120.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed005_SimBKG_CLUseed005_processing_0150.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed005_SimBKG_CLUseed005_processing_0180.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed005_SimBKG_CLUseed005_processing_0210.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed005_SimBKG_CLUseed005_processing_0240.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed005_SimBKG_CLUseed005_processing_0270.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed005_SimBKG_CLUseed005_processing_0300.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed005_SimBKG_CLUseed005_processing_0330.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed005_SimBKG_CLUseed005_processing_0360.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed005_SimBKG_CLUseed005_processing_0390.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed005_SimBKG_CLUseed005_processing_0420.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed005_SimBKG_CLUseed005_processing_0450.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed005_SimBKG_CLUseed005_processing_0480.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed005_SimBKG_CLUseed005_processing_0510.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed005_SimBKG_CLUseed005_processing_0540.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed005_SimBKG_CLUseed005_processing_0570.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed005_SimBKG_CLUseed005_processing_0600.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed005_SimBKG_CLUseed005_processing_0630.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed005_SimBKG_CLUseed005_processing_0660.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed006_SimBKG_CLUseed006_processing_0000.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed006_SimBKG_CLUseed006_processing_0030.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed006_SimBKG_CLUseed006_processing_0060.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed006_SimBKG_CLUseed006_processing_0090.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed006_SimBKG_CLUseed006_processing_0120.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed006_SimBKG_CLUseed006_processing_0150.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed006_SimBKG_CLUseed006_processing_0180.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed006_SimBKG_CLUseed006_processing_0210.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/GE_e5_merge_AGNseed006_SimBKG_CLUseed006_processing_0240.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
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
