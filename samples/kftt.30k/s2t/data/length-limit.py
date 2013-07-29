#!/usr/bin/env python

from optparse import OptionParser

parser = OptionParser()

parser.add_option("--source",   action="store", type="string", default="")
parser.add_option("--target",   action="store", type="string", default="")
parser.add_option("--output-source",   action="store", type="string", default="")
parser.add_option("--output-target",   action="store", type="string", default="")

(options, args) = parser.parse_args()

in_src = open(options.source)
in_trg = open(options.target)

out_src = open(options.output_source, 'w')
out_trg = open(options.output_target, 'w')

id = 0
while 1:

    source = in_src.readline()
    target = in_trg.readline()

    if not source or not target: break
    
    if len(source.split()) > 20: continue

    out_src.write(source)
    out_trg.write(str(id) + ' ||| ' + target)
    
    id += 1
