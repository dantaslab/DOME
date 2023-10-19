#!/usr/bin/env python3

"""Filter out the duplicates (A B vs B A) from the joint BLASTn output file


Usage: python3 s22_nodupBLASTn.py <BLASTn file> <filtered output>

Args:
        BLASTn file: joint BLASTn output file; .txt file
  filtered output: outputs with duplicates filtered out; .txt file
"""

import subprocess
import os
import sys

class DupFilter(object):
    def __init__(self, infile, outfile):
        self.infile = infile
        self.outfile = outfile
        self.set1 = '{}.set1'.format(outfile)
        self.set2 = '{}.set2'.format(outfile)
        self.set1sorted = '{}.sorted'.format(self.set1)
        self.set2sorted = '{}.sorted'.format(self.set2)
    
    def get_key(self, line, reverse=False):
        if line == b'':
            return None
        if not reverse:
            return line.split(b'\t',3)[0:2]
        else:
            return line.split(b'\t',3)[1::-1]

    def step1(self):
        with open(self.set1, "wb") as t1, open(self.set2, 'wb') as t2, open(self.outfile, 'wb') as f, open('single.debug', 'wb') as u:
            for line in open(self.infile, "rb"):
              pair = line.split(b'\t',3)[:2]
              if pair[0] > pair[1]:
                  t1.write(line)
              elif pair[0] < pair[1]:
                  t2.write(line)
              else:
                  u.write(line)
                  f.write(line)

    def step2(self):
        try:
            p1 = subprocess.Popen(['/usr/bin/sort', '-k', '1,1', '-k', '2,2', '-o', self.set1sorted, self.set1])
            p2 = subprocess.Popen(['/usr/bin/sort', '-k', '2,2', '-k', '1,1', '-o', self.set2sorted, self.set2])
            p1.wait()
            p2.wait()
        except subprocess.CalledProcessError as e:
            raise(e)

    def step3(self):
        with open(self.set1sorted, "rb") as t1, open(self.set2sorted, 'rb') as t2, open(self.outfile, 'r+b') as f, open('unique.debug', 'wb') as u, open('paired.debug', 'wb') as p:
            f.seek(0,2)
            l1 = t1.readline()
            l2 = t2.readline()
            while True:
                k1,k2 = self.get_key(l1), self.get_key(l2, True)
                if k1 == k2 == None:
                    break
                elif k1 == k2:
                    f.write(l1)
                    p.write(l1)
                    l1 = t1.readline()
                    l2 = t2.readline()
                elif k2 is None or k1 < k2:
                    f.write(l1)
                    u.write(l1)
                    l1 = t1.readline()
                elif k1 is None or k2 < k1:
                    f.write(l2)
                    u.write(l2)
                    l2 = t2.readline()
    def step4(self):
        os.unlink(self.set1)
        os.unlink(self.set2)
        os.unlink(self.set1sorted)
        os.unlink(self.set2sorted)

if __name__ == '__main__':
    import time

    #Print out the doc string and exit if the number of input parameters is not correct
    if len(sys.argv) != 3:
      print(__doc__)
      exit(1)
    
    infile = sys.argv[1]
    outfile = sys.argv[2]
    d = DupFilter(infile, outfile)
    t1 = time.time()
    d.step1()
    print("Step1: ", time.time() - t1)
    t2 = time.time()
    d.step2()
    print("Step2: ", time.time() - t2)
    t3 = time.time()
    d.step3()
    print("Step3: ", time.time() - t3)
    t4 = time.time()
    d.step4()
    print("Total: ", time.time() - t1)