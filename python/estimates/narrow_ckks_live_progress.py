#!/usr/bin/env python3
"""Observe existing Linux CKKS jobs without reading keys or controlling jobs.

Process activity and file updates are observations, not accuracy, completion,
or remaining-time estimates. PID/start-time pairs protect against PID reuse.
"""
import argparse
import datetime
import json
import os
from pathlib import Path
import time


def observe_process(pid, expected_start):
    try:
        fields = (Path('/proc') / str(pid) / 'stat').read_text().rsplit(')', 1)[1].split()
        start = int(fields[19])
        if start != expected_start:
            return dict(pid=pid, expected_start_ticks=expected_start,
                        state='identity_mismatch', observed_start_ticks=start)
        return dict(pid=pid, start_ticks=start, state=fields[0],
                    cpu_ticks=int(fields[11]) + int(fields[12]),
                    rss_bytes=int(fields[21]) * os.sysconf('SC_PAGE_SIZE'))
    except FileNotFoundError:
        return dict(pid=pid, expected_start_ticks=expected_start, state='absent')
    except (OSError, ValueError, IndexError) as error:
        return dict(pid=pid, expected_start_ticks=expected_start,
                    state='observation_error', error=str(error))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--job', action='append', required=True, metavar='NAME:PID:START_TICKS')
    parser.add_argument('--artifact', type=Path, action='append', default=[])
    parser.add_argument('--output', type=Path, required=True, help='New JSONL history file')
    parser.add_argument('--summary', type=Path, required=True, help='New latest-observation JSON file')
    parser.add_argument('--interval-seconds', type=float, default=30)
    parser.add_argument('--duration-seconds', type=float, default=1800)
    args = parser.parse_args()
    if not (0 < args.interval_seconds <= 60 and 0 <= args.duration_seconds <= 3600):
        parser.error('interval must be in (0,60], duration in [0,3600]')
    if (args.output.exists() or args.summary.exists() or
            args.output.resolve() == args.summary.resolve()):
        parser.error('history and summary must be distinct new files')
    jobs = {}
    try:
        for entry in args.job:
            name, pid, start = entry.rsplit(':', 2)
            if not name or name in jobs or int(pid) <= 0 or int(start) < 0:
                raise ValueError('invalid or duplicate job')
            jobs[name] = int(pid), int(start)
    except ValueError as error:
        parser.error(str(error))
    os.umask(0o077)
    ticks_per_second = os.sysconf('SC_CLK_TCK')
    previous = {}
    start_time = time.monotonic()
    previous_time = start_time
    with args.output.open('x') as history:
        while True:
            now = time.monotonic()
            current = {}
            for name, (pid, start) in jobs.items():
                row = observe_process(pid, start)
                old = previous.get(name, {})
                if 'cpu_ticks' in row and 'cpu_ticks' in old:
                    advanced = row['cpu_ticks'] - old['cpu_ticks']
                    if advanced >= 0:
                        row['cpu_seconds_since_previous'] = advanced / ticks_per_second
                current[name] = row
            artifacts = {}
            for path in args.artifact:
                try:
                    info = path.stat()
                    artifacts[str(path)] = dict(bytes=info.st_size, inode=info.st_ino,
                                                mtime_ns=info.st_mtime_ns)
                except FileNotFoundError:
                    artifacts[str(path)] = dict(state='absent')
                except OSError as error:
                    artifacts[str(path)] = dict(state='observation_error', error=str(error))
            observation = dict(observed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                               elapsed_seconds=now-start_time,
                               interval_seconds=now-previous_time, jobs=current,
                               artifacts=artifacts,
                               scope='Read-only activity observation; does not establish test accuracy or completion')
            payload = json.dumps(observation, allow_nan=False)
            history.write(payload+'\n')
            history.flush()
            pending = args.summary.with_suffix(args.summary.suffix+'.pending')
            with pending.open('x') as file:
                file.write(payload+'\n')
            pending.replace(args.summary)
            previous = current
            previous_time = now
            remaining = args.duration_seconds - (time.monotonic()-start_time)
            if remaining <= 0:
                break
            time.sleep(min(args.interval_seconds, remaining))


if __name__ == '__main__':
    main()
