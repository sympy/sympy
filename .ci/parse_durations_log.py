#!/usr/bin/env python3

import os
import json

ci_folder = os.path.dirname(__file__)

def read_log():
    stop_token = '= FAILURES ='
    start_token = '= slowest test durations ='
    start_token_seen = False
    for line in open(os.path.join(ci_folder, 'durations.log')):
        if stop_token in line:
            return
        elif start_token_seen:
            time, kind, test_id = line.split()
            if kind != 'call':
                continue
            if time[-1] != 's':
                raise NotImplementedError("expected seconds")
            yield test_id, float(time[:-1])
        elif start_token in line:
            start_token_seen = True


def main(limits=(10, .1)):
    """
    parses durations.log (see generate_durations_log.sh for how to generate that file)
    """
    groupings = [[] for _ in range(len(limits))]
    for test_id, time in read_log():
        for idx, lim in enumerate(limits):
            if time >= lim:
                groupings[idx].append(test_id)
    open(os.path.join(ci_folder, 'durations.json'), 'wt').write(json.dumps([sorted(gr) for gr in groupings]))


if __name__ == '__main__':
    main()
