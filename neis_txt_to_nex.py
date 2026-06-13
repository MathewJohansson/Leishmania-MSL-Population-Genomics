#!/usr/bin/env python3
import sys
import re

def convert(infile, outfile):
    with open(infile, 'r') as f:
        lines = [l.rstrip('\n') for l in f.readlines()]
    ntax = int(lines[0].strip().split()[0])
    samples = []
    rows = []
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        parts = re.split(r'\s+', line)
        samples.append(parts[0])
        rows.append(parts[1:])
    with open(outfile, 'w') as f:
        f.write('#NEXUS\n')
        f.write('BEGIN Taxa;\n')
        f.write(f'DIMENSIONS ntax={ntax};\n')
        f.write('TAXLABELS\n')
        for i, s in enumerate(samples, 1):
            f.write(f"[{i}] '{s}'\n")
        f.write(';\n')
        f.write('END;\n')
        f.write('BEGIN Distances;\n')
        f.write(f'DIMENSIONS ntax={ntax};\n')
        f.write('FORMAT labels=left diagonal=yes missing=? triangle=lower;\n')
        f.write('MATRIX\n')
        for i, (s, row) in enumerate(zip(samples, rows)):
            lower = row[:i + 1]
            vals = '  '.join(lower)
            f.write(f"'{s}' {vals}\n")
        f.write(';\n')
        f.write('END;\n')
    print(f"Done: {outfile} ({ntax} taxa)")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python3 neis_txt_to_nex.py input.txt output.nex")
        sys.exit(1)
    convert(sys.argv[1], sys.argv[2])
