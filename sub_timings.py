from sage.all import *
proof.all(False)

from sub_ftf import sub_FTF
from params.params import sub_Params

import time
import logging
# usage: logging.getLogger('sub_ftf').setLevel(logging.DEBUG/INFO)
logging.getLogger('sub_ftf').setLevel(logging.INFO)
logging.getLogger('sub_orientation').setLevel(logging.INFO)

def run_times(sec):

    params = sub_Params(sec)
    alice = sub_FTF(params)
    B = params.B

    n_runs = 10
    action_alice = []

    for run in range(n_runs):
        print(f'Run {run+1}')
        alice = sub_FTF(params)

        t0 = time.time()

        # Do the action
        sigma = alice.action()
        t1 = time.time()

        print(f'Action: {t1-t0:.3f}s')
        print(f'-'*50)

        action_alice.append(t1-t0)


    avg_action = sum(action_alice)/n_runs
    print(f'='*50)
    print(f'{sec = }')
    print(f'Result of {n_runs} runs')
    print(f'Avg action: {avg_action:.3f}s')
    with open('out.txt', 'a') as fh:
        fh.write(f'='*50 + '\n')
        fh.write(f'{sec = }\n')
        fh.write(f'Result of {n_runs} runs' + '\n')
        fh.write(f'Avg action: {avg_action:.3f}s' + '\n')

# run_times(333)
run_times(1184)
