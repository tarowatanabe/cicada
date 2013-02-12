#!/usr/bin/env python

import sys
import re

opt = re.compile(r"^(\s+)(--\S+)\s+(arg \(=.+?\)|\[=arg\(=\S+\)\]|arg|)(.*?)$")

for line in sys.stdin:
    result = opt.search(line)

    if not result:
        sys.stdout.write(line)
        continue
    
    sys.stdout.write('\n')
    if result.group(3):
        sys.stdout.write("%s**%s** `%s` %s\n" %(result.group(1), result.group(2), result.group(3), result.group(4)))
    else:
        sys.stdout.write("%s**%s** %s\n" %(result.group(1), result.group(2), result.group(4)))
    
    

