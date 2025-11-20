import os
import sys
import subprocess
from select import select

# Start two subprocesses
procs = [subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/make_summary_skymap.py", "GE_e4_merge_AGNseed00{0}_SimBKG".format(seed)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True) for seed in range(10)] + [subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/make_summary_skymap.py", "GE_e4_merge_AGNseed00{0}_SimBKG_CLUseed00{0}".format(seed)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True) for seed in range(10)] + [subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/make_summary_skymap.py", "GE_e4_merge_SimBKG_CLUseed00{0}".format(seed)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True) for seed in range(10)] + [subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/make_summary_skymap.py", "GE_e4_merge_SimBKG"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True)] + [subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/make_summary_skymap.py", "GE_e5_merge_AGNseed00{0}_SimBKG_CLUseed00{0}".format(seed)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True) for seed in range(10)] + [subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/make_summary_skymap.py", "GE_e5_mergefagn0p15_AGNseed00{0}_SimBKG_CLUseed00{0}".format(seed)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True) for seed in range(19)]

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
