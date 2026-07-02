import os
import sys
import subprocess
from select import select

#Define seed list
clu_seed = list(range(1,118))
agn_seed = list(range(1,10))*13
t_tot_seed = list(zip(agn_seed, clu_seed))
tot_seed = [(str(seed[0]).zfill(3), str(seed[1]).zfill(3)) for seed in t_tot_seed]

# Start two subprocesses
procs = [subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/make_summary_skymap.py", "GE_e4_merge_AGNseed{0}_SimBKG_CLUseed{1}".format(seed[0], seed[1])], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True) for seed in tot_seed] + [subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/make_summary_skymap.py", "GE_e5_merge_AGNseed{0}_SimBKG_CLUseed{1}".format(seed[0], seed[1])], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True) for seed in tot_seed]

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