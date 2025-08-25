import os
import sys
import subprocess
from select import select

# Start two subprocesses
procs = [
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e4_merge_AGNseed002_SimBKG"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e4_merge_AGNseed002_SimBKG_CLUseed002"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e4_merge_SimBKG_CLUseed002"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e4_merge_AGNseed003_SimBKG"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e4_merge_AGNseed003_SimBKG_CLUseed003"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e4_merge_SimBKG_CLUseed003"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e4_merge_AGNseed004_SimBKG"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e4_merge_AGNseed004_SimBKG_CLUseed004"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e4_merge_SimBKG_CLUseed004"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e4_merge_AGNseed005_SimBKG"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e4_merge_AGNseed005_SimBKG_CLUseed005"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e4_merge_SimBKG_CLUseed005"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e4_merge_AGNseed006_SimBKG"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e4_merge_AGNseed006_SimBKG_CLUseed006"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e4_merge_SimBKG_CLUseed006"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e4_merge_AGNseed007_SimBKG"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e4_merge_AGNseed007_SimBKG_CLUseed007"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e4_merge_SimBKG_CLUseed007"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e4_merge_AGNseed008_SimBKG"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e4_merge_AGNseed008_SimBKG_CLUseed008"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e4_merge_SimBKG_CLUseed008"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e4_merge_SimBKG"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e4_merge_AGNseed001_SimBKG"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e4_merge_AGNseed001_SimBKG_CLUseed001"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e4_merge_SimBKG_CLUseed001"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e5_merge_AGNseed001_SimBKG_CLUseed001"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e5_merge_AGNseed002_SimBKG_CLUseed002"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e5_merge_AGNseed003_SimBKG_CLUseed003"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e5_merge_AGNseed004_SimBKG_CLUseed004"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e5_merge_AGNseed005_SimBKG_CLUseed005"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e5_merge_AGNseed006_SimBKG_CLUseed006"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e5_merge_AGNseed007_SimBKG_CLUseed007"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e5_merge_AGNseed008_SimBKG_CLUseed008"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e5_mergefagn0p15_AGNseed001_SimBKG_CLUseed001"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e5_mergefagn0p15_AGNseed002_SimBKG_CLUseed002"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e5_mergefagn0p15_AGNseed003_SimBKG_CLUseed003"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e5_mergefagn0p15_AGNseed004_SimBKG_CLUseed004"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e5_mergefagn0p15_AGNseed005_SimBKG_CLUseed005"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e5_mergefagn0p15_AGNseed006_SimBKG_CLUseed006"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e5_mergefagn0p15_AGNseed007_SimBKG_CLUseed007"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", "GE_e5_mergefagn0p15_AGNseed008_SimBKG_CLUseed008"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
    # subprocess.Popen(["python", "/home/idies/workspace/erosim/software/st_mod/src/esass/runs/make_summary_skymap.py", ""], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, text=True, close_fds=True),
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
