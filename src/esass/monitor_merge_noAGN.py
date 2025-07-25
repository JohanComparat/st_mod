import os
import sys
import subprocess
from select import select

# Start two subprocesses
procs = [
    subprocess.Popen(["python", "merge_events_noAGN.py", "1"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                     bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "merge_events_noAGN.py", "2"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                     bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "merge_events_noAGN.py", "3"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                     bufsize=1, text=True, close_fds=True),
    subprocess.Popen(["python", "merge_events_noAGN.py", "4"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                     bufsize=1, text=True, close_fds=True),
    # subprocess.Popen(["python", "merge_events_noAGN.py", "5"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    #                  bufsize=1, text=True, close_fds=True),
    # subprocess.Popen(["python", "merge_events_noAGN.py", "6"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    #                  bufsize=1, text=True, close_fds=True),
    # subprocess.Popen(["python", "merge_events_noAGN.py", "7"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    #                  bufsize=1, text=True, close_fds=True),
    # subprocess.Popen(["python", "merge_events_noAGN.py", "8"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    #                  bufsize=1, text=True, close_fds=True),
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
