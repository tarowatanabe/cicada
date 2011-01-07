#!/usr/bin/env python
#
#  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
#

import sys
import os.path
import subprocess

from optparse import OptionParser, make_option

opt_parser = OptionParser(
    option_list=[
        make_option("--root-count", action="store", type="string", default="-",
                    help="root count file", metavar='FILE')
        ])

(options, args) = opt_parser.parse_args()

def istream(name):
    if name == '-':
        return sys.stdin
    
    (prefix, ext) = os.path.splitext(name)
    
    if ext == '.gz':
        return subprocess.Popen(['gunzip', '-c', name], stdout=subprocess.PIPE).stdout
    elif ext == '.bz2':
        return subprocess.Popen(['bunzip2', '-c', name], stdout=subprocess.PIPE).stdout
    else:
        return open(name, 'r')

def is_non_terminal(non_terminal):
    return len(non_terminal) > 2 and non_terminal[0] == '[' and non_terminal[-1] == ']'

non_terminals = set()
for line in istream(options.root_count):
    tokens = line.split()
    if not tokens: continue

    if not is_non_terminal(tokens[0]): continue

    non_terminals.add(tokens[0])

for non_terminal in non_terminals:
    print non_terminal
    
