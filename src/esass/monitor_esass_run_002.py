import os
import sys
import subprocess
from select import select

# Start two subprocesses
procs = [
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed049_processing_1600.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed049_processing_1625.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed049_processing_1650.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed049_processing_1675.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed049_processing_1700.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed049_processing_1725.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed049_processing_1750.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed049_processing_1775.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed049_processing_1800.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed049_processing_1825.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed049_processing_1850.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed049_processing_1875.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed049_processing_1900.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed049_processing_1925.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed049_processing_1950.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed049_processing_1975.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed049_processing_2000.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed049_processing_2025.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed049_processing_2050.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed049_processing_2075.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed049_processing_2100.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed049_processing_2125.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed049_processing_2150.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed049_processing_2175.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed049_processing_2200.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed049_processing_2225.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed049_processing_2250.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed049_processing_2275.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed049_processing_2300.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed049_processing_2325.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed049_processing_2350.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["bash", "/home/idies/workspace/erosim/runs/GE_e4_merge_AGNseed004_SimBKG_CLUseed049_processing_2375.sh"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
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