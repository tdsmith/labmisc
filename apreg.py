#!/usr/bin/env python
"Regularizes .csv exports from Anton-Parr RheoPlus software."

import unicodecsv as csv
import codecs
import re
import pandas as pd

def cast(field):
    try:
        field = int(field)
        return field
    except ValueError:
        pass
    try:
        field = float(field)
        return field
    except ValueError:
        return field

def parse(fp):
    # guess delimiter
    c = fp.read(1)
    while c and c not in ",\t":
        c = fp.read(1)
    fp.seek(0)

    reader = csv.reader(fp, delimiter=c or ",", encoding='latin-1')

    in_data_table = False
    frames = []

    # record will be a dict of lists {'column_name': [row0, row1...], ...}
    record = None
    record_name = None
    header = None
    units_next = False
    in_data_table = False
    pt_duration_next = False
    pt_duration = 0
    pt_duration_unit = ""

    def closeout(frame, record_name):
        realname, rep = record_name.rsplit(None, 1)
        frame['Experiment'] = realname
        frame['Rep'] = rep
        time_col = "Time [%s]" % pt_duration_unit
        interval_col = "Interval [%s]" % pt_duration_unit
        if time_col not in frame.columns:
            frame[time_col] = frame["Meas. Pts."].astype(int) * pt_duration
            frame[interval_col] = pt_duration

    for line in reader:
        if not line or line == [u'']:
            continue

        if line[0] == 'Data Series Information' or line[0].startswith(u'Interval:'):
            # new record; terminate last record
            if record:
                frame = pd.DataFrame(record)
                closeout(frame, record_name)
                frames.append(frame)
            # reset state
            record = {}
            in_data_table, units_next = False, False
            header = None

        if units_next:
            units_next = False
            in_data_table = True
            for i, unit in enumerate(line):
                if unit:
                    header[i] = '%s %s' % (header[i], unit)
            continue

        if in_data_table:
            row = {}
            # protect against repeated headers
            seen_headers = set()
            for i, cell in enumerate(line):
                if header[i] in seen_headers:
                    continue
                seen_headers.add(header[i])
                if cell == "******":
                    cell = ""
                record.setdefault(header[i], []).append(cell)

        if pt_duration_next:
            pt_duration_next = False
            for cell in line:
                if not cell.startswith("Meas. Pt."):
                    continue
                break
            m = re.match(r"Meas\. Pt\. Duration ([\d.]+) (\w+)", cell)
            pt_duration, pt_duration_unit = m.groups()
            pt_duration = float(pt_duration)

        if line[0] == 'Name:':
            record_name = line[-1]

        if line[0] == 'Meas. Pts.':
            units_next = True
            header = line

        if line[0] == 'Time Setting:':
            pt_duration_next = True

    # off the end of the file; close out
    frame = pd.DataFrame(record)
    closeout(frame, record_name)
    frames.append(frame)

    return pd.tools.merge.concat(frames, ignore_index=True)

if __name__ == '__main__':
    import sys
    import os.path
    for filename in sys.argv[1:]:
        with open(filename, mode='rb') as fp:
            table = parse(fp)
        table.to_csv('regularized-' + os.path.basename(filename),
                na_rep = 'NA',
                index = False,
                encoding = 'utf-8')
