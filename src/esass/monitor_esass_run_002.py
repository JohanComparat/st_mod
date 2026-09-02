import os
import sys
import subprocess
from select import select

# Start two subprocesses
procs = [
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed009_SimBKG_CLUseed054_processing_0000.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed009_SimBKG_CLUseed054_processing_0050.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed009_SimBKG_CLUseed054_processing_0100.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed009_SimBKG_CLUseed054_processing_0150.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed009_SimBKG_CLUseed054_processing_0200.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed009_SimBKG_CLUseed054_processing_0250.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed009_SimBKG_CLUseed054_processing_0300.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed009_SimBKG_CLUseed054_processing_0350.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed009_SimBKG_CLUseed054_processing_0400.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed009_SimBKG_CLUseed054_processing_0450.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed009_SimBKG_CLUseed054_processing_0500.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed009_SimBKG_CLUseed054_processing_0550.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed009_SimBKG_CLUseed054_processing_0600.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed009_SimBKG_CLUseed054_processing_0650.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed009_SimBKG_CLUseed054_processing_0700.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed009_SimBKG_CLUseed054_processing_0750.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed009_SimBKG_CLUseed054_processing_0800.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed009_SimBKG_CLUseed054_processing_0850.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed009_SimBKG_CLUseed054_processing_0900.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed009_SimBKG_CLUseed054_processing_0950.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed009_SimBKG_CLUseed054_processing_1000.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed009_SimBKG_CLUseed054_processing_1050.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed009_SimBKG_CLUseed054_processing_1100.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed009_SimBKG_CLUseed054_processing_1150.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed009_SimBKG_CLUseed054_processing_1200.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed009_SimBKG_CLUseed054_processing_1250.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed009_SimBKG_CLUseed054_processing_1300.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed009_SimBKG_CLUseed054_processing_1350.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed009_SimBKG_CLUseed054_processing_1400.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed009_SimBKG_CLUseed054_processing_1450.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed009_SimBKG_CLUseed054_processing_1500.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed009_SimBKG_CLUseed054_processing_1550.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
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

/home/idies/workspace/erosim/runs/GE_e5_merge_AGNseed005_SimBKG_CLUseed023_processing_0000.sh
/home/idies/workspace/erosim/runs/GE_e5_merge_AGNseed005_SimBKG_CLUseed023_processing_0005.sh
/home/idies/workspace/erosim/runs/GE_e5_merge_AGNseed005_SimBKG_CLUseed023_processing_0010.sh
/home/idies/workspace/erosim/runs/GE_e5_merge_AGNseed005_SimBKG_CLUseed023_processing_0015.sh
/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed003_SimBKG_CLUseed039_processing_1440.sh
/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed003_SimBKG_CLUseed039_processing_1485.sh
/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed003_SimBKG_CLUseed039_processing_1530.sh
/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed003_SimBKG_CLUseed039_processing_1575.sh
/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed005_SimBKG_CLUseed032_processing_1620.sh
/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed005_SimBKG_CLUseed032_processing_1665.sh
/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed005_SimBKG_CLUseed032_processing_1710.sh
/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed005_SimBKG_CLUseed032_processing_1755.sh
/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed005_SimBKG_CLUseed032_processing_1800.sh
/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed005_SimBKG_CLUseed032_processing_1845.sh
/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed005_SimBKG_CLUseed032_processing_1890.sh
/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed005_SimBKG_CLUseed032_processing_1935.sh
/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed005_SimBKG_CLUseed032_processing_1980.sh
/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed005_SimBKG_CLUseed032_processing_2025.sh
/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed005_SimBKG_CLUseed032_processing_2070.sh
/home/idies/workspace/erosim/runs/GE_e5_merge_AGNseed002_SimBKG_CLUseed038_processing_0000.sh
/home/idies/workspace/erosim/runs/GE_e5_merge_AGNseed002_SimBKG_CLUseed038_processing_0045.sh
/home/idies/workspace/erosim/runs/GE_e5_merge_AGNseed002_SimBKG_CLUseed038_processing_0090.sh
/home/idies/workspace/erosim/runs/GE_e5_merge_AGNseed002_SimBKG_CLUseed038_processing_0135.sh
/home/idies/workspace/erosim/runs/GE_e5_merge_AGNseed002_SimBKG_CLUseed038_processing_0180.sh
/home/idies/workspace/erosim/runs/GE_e5_merge_AGNseed002_SimBKG_CLUseed038_processing_0225.sh
/home/idies/workspace/erosim/runs/GE_e5_merge_AGNseed002_SimBKG_CLUseed038_processing_0270.sh
/home/idies/workspace/erosim/runs/GE_e5_merge_AGNseed002_SimBKG_CLUseed038_processing_0315.sh
/home/idies/workspace/erosim/runs/GE_e5_merge_AGNseed002_SimBKG_CLUseed038_processing_0360.sh
/home/idies/workspace/erosim/runs/GE_e5_merge_AGNseed002_SimBKG_CLUseed038_processing_0405.sh
/home/idies/workspace/erosim/runs/GE_e5_merge_AGNseed002_SimBKG_CLUseed038_processing_0450.sh
/home/idies/workspace/erosim/runs/GE_e5_merge_AGNseed002_SimBKG_CLUseed038_processing_0495.sh
/home/idies/workspace/erosim/runs/GE_e5_merge_AGNseed002_SimBKG_CLUseed038_processing_0540.sh