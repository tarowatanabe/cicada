#!/usr/bin/env python

import sys
import re

pattern = re.compile(r"^\s*\"(.+?)\s+<>\s+(.+?)\s*;\s*\",\s*$")

entity_map = {}
uspatent_map = {}

for line in open('sgml_entity_table.hpp'):
    result = pattern.search(line)

    if not result: continue

    entity = result.group(1)
    ustr   = result.group(2)
    
    entity_map[entity] = ustr


for line in open('sgml_uspatent_table.hpp'):
    result = pattern.search(line)

    if not result: continue

    entity = result.group(1)
    ustr   = result.group(2)
    
    uspatent_map[ustr] = entity

uspatents = uspatent_map.keys()

def compare(x, y):
    if len(x) != len(y):
        return cmp(len(y), len(x))
    else:
        return cmp(x, y)

uspatents.sort(compare)

for uspatent in uspatents:

    entity = uspatent_map[uspatent]
    
    if not entity_map.has_key(entity):
        sys.stderr.write("no key? " + uspatent+' '+entity+'\n')
        continue
    
    sys.stdout.write("\"%s <> %s;\",\n" %(uspatent, entity_map[entity]))
