#!/usr/bin/env python

import sys

for i in range(256):
    
    sys.stdout.write("/* 0x%02x */" %(i))

    offset = 0
    sys.stdout.write(" %d," %(offset))
    for j in range(4):
        offset += ((i >> (6 - j * 2)) & 0x03) + 1
        sys.stdout.write(" %d," %(offset))

    sys.stdout.write('\n')
