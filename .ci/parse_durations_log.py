#!/usr/bin/env python3

from collections import defaultdict
import os
import json
import time

ci_folder = os.path.dirname(__file__)


def read_log():
    start_token = '= slowest test durations ='
    start_token_seen = False
    for line in open(os.path.join(ci_folder, 'durations.log')):
        if start_token_seen:
            try:
                dur, kind, test_id = line.split()
            except:
                return
            else:
                if dur[0] not in '0123456789':
                    return
            if kind != 'call':
                continue
            if dur[-1] != 's':
                raise NotImplementedError("expected seconds")
            yield test_id, float(dur[:-1])
        elif start_token in line:
            start_token_seen = True


def main(ref_timing, limits=(10, .1)):
    """
    parses durations.log (made by generate_durations_log.sh)
    """
    groupings = [defaultdict(list) for _ in range(len(limits))]
    accumul_n = [0 for _ in range(len(limits))]
    accumul_t = [0.0 for _ in range(len(limits))]
    for test_id, dur in read_log():
        if test_id.startswith('sympy/utilities/tests/test_code_quality.py'):
            continue # white-listed (worth running since it catches many errors)
        for idx, lim in enumerate(limits):
            if dur/ref_timing >= lim:
                fname, tname = test_id.split('::')
                groupings[idx][fname].append(tname)
                accumul_t[idx] += dur
                accumul_n[idx] += 1
                break
    json_data = json.dumps([{k: sorted(v) for k, v in gr.items()}
                            for gr in groupings], indent=4, sort_keys=True)
    open(os.path.join(ci_folder, 'durations.json'), 'wt').write(json_data)
    print('number in group, accumulated_time: %s' %
          str(list(zip(accumul_n, accumul_t))))


def slow_function():
    t = time.time()
    a = 0
    for i in range(5):
        a += sum(x**.3 - x**i for x in range(1000000) if x % 3 == 0)
    return time.time() - t


if __name__ == '__main__':
    ref_time = slow_function()
    main(ref_time)
