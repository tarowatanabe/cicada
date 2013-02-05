#!/usr/bin/env python

import sys
import re

entities = {}

pattern = re.compile(r"^\s*<!ENTITY\s+(\S+)\s+(\S+)\s*>")

pattern_amp = re.compile(r"&#38;")
pattern_plane = re.compile(r"%plane([a-fA-F0-9]+);")

for line in sys.stdin:
    
    result = pattern.search(line)
    if not result: continue

    name = result.group(1)
    mapped = result.group(2)
    
    if mapped and mapped[0] == '"':
        mapped = mapped[1:]
    if mapped and mapped[-1] == '"':
        mapped = mapped[:-1]

    if name and mapped:
        
        mapped = pattern_amp.sub(r"&", mapped)
        mapped = pattern_plane.sub(r"&#x0\1", mapped)
        
        entities[name] = mapped

pattern_entity = re.compile(r"&#([a-zA-Z0-9]+);")

inverse = {}
for name in entities.keys():
    umapped = ""
    for mapped in pattern_entity.findall(entities[name]):
        value = 0
        if mapped[0] == 'x':
            value = int(mapped[1:], 16)
        else:
            value = int(mapped)
        if value > 0xffff:
            umapped += "\\\\U%08x" %(value)
        else:
            umapped += "\\\\u%04x" %(value)

    inverse[umapped] = name

maps = inverse.keys()
maps.sort()
maps.reverse()

for mapped in maps:
    sys.stdout.write("\"%s <> '&%s;';\",\n" %(mapped, inverse[mapped]))
            
        


    
    
