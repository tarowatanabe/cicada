#!/usr/bin/env python

import sys
import re

entities = {}

pattern = re.compile(r"^\s*\"(.+?)\",\s+\"(.+?)\"")


for line in sys.stdin:
    
    result = pattern.search(line)
    if not result: continue
    
    entity = result.group(1)
    mapped = result.group(2)

    if len(mapped) <= 1 or mapped[0] != '.' or mapped[-1] != '.': continue
    if entity[0] != '&' or entity[-1] != ';': continue

    mapped = mapped.replace('\\', '')
    
    entities[entity] = mapped

names = entities.keys()

def compare(x, y):
    if len(x) != len(y):
        return cmp(len(y), len(x))
    else:
        return cmp(x, y)

names.sort(compare)

for name in names:
    sys.stdout.write("\"'%s' <> '%s';\",\n" %(name, entities[name]))
    
