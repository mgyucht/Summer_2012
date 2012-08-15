#!/usr/bin/env python2

for p0 in range(23, 40):
  p = round(p0 / float(50), 2)
  out = '/scratch/gpfs/myucht/%.2g/%.4g/' % (p, 0.0001259)
  print('%sstress_data_%d.txt' % (out, 1))
