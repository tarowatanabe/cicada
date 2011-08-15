#!/usr/bin/env python
#
#  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
#

###
### a cicada.config generator
###

###
### We will generate common cicada.config from indexed grammar
### (SCFG and Tree-to-string)
### 

import UserList
import os
import os.path

from optparse import OptionParser, make_option

opt_parser = OptionParser(
    option_list=[
    
    ### grammars
    make_option("--grammar",      default=[], action="append", type="string", help="indexed grammar"),
    make_option("--tree-grammar", default=[], action="append", type="string", help="indexed tree-grammar"),
    
    ### only affect grammars
    make_option("--max-span", default=15, action="store", type="int",
                metavar="LENGTH", help="maximum span size (default: 15)"),
    
    ### glues
    make_option("--goal", default="[s]", action="store", type="string", help="goal non-terminal (default: [s])"),
    make_option("--glue", default="[x]", action="store", type="string", help="non-terminal for glue rules (default: [x])"),
    
    make_option("--straight", default=None, action="store_true", help="straight gulue rule"),
    make_option("--invert",   default=None, action="store_true", help="invert gulue rule"),
    
    make_option("--insertion", default=None, action="store_true", help="insertion grammar"),
    make_option("--deletion",  default=None, action="store_true", help="deletion grammar"),
    make_option("--fallback",  default=None, action="store_true", help="fallback tree-grammar"),
    
    ### feature functions

    ## operations...
    
    # cicada composition
    make_option("--phrase",   default=None, action="store_true", help="phrase-based grammar"),
    make_option("--scfg",     default=None, action="store_true", help="SCFG"),
    make_option("--tree",     default=None, action="store_true", help="tree-to-string"),
    make_option("--tree-cky", default=None, action="store_true", help="string-to-{string,tree}"),
    
    ### outputs
    make_option("--file",      default=None, action="store_true", help="output in a file"),
    make_option("--directory", default=None, action="store_true", help="output in a directory"),
    
    ## debug messages
    make_option("--debug", default=0, action="store", type="int"),
    ])


class Grammar(UserList.UserList):
    
    def __init__(self, grammar_dir="", max_span=15):
        
        UserList.UserList.__init__(self)
        
        path_files = os.path.join(grammar_dir, 'files');
        
        if not os.path.exists(path_files):
            raise ValueError, "no path to files: %s" %(path_files)
        
        self.grammar_dir = grammar_dir
        
        for line in open(path_files):
            
            name = line.strip()
            if not name: continue
            
            path = os.path.join(grammar_dir, name)
            if not os.path.exists(path):
                raise ValueError, "no path to scores: %s" %(path)
            
            self.append("grammar = " + path + ":max-span=%d" %(max_span))

class TreeGrammar(UserList.UserList):
    
    def __init__(self, grammar_dir="", max_span=15):
        
        UserList.UserList.__init__(self)
        
        path_files = os.path.join(grammar_dir, 'files');
        
        if not os.path.exists(path_files):
            raise ValueError, "no path to files: %s" %(path_files)
        
        self.grammar_dir = grammar_dir
        
        for line in open(path_files):
            
            name = line.strip()
            if not name: continue
            
            path = os.path.join(grammar_dir, name)
            if not os.path.exists(path):
                raise ValueError, "no path to scores: %s" %(path)
            
            self.append("tree-grammar = " + path)

(options, args) = opt_parser.parse_args()

### grammars

def non_terminal(x, index=0):
    if not x:
        raise ValueError, "no non-terminal? %s"  %(x)
    if x[0] != '[' or x[-1] != ']':
        raise ValueError, "invalid non-terminal? %s" %(x)
    if index == 0:
        return x
    else:
        return '[' + x[1:-1] + ',' + str(index) + ']'


options.goal = non_terminal(options.goal)
options.glue = non_terminal(options.glue)

print "goal = %s" %(options.goal)
print

if options.grammar:
    print "# grammar"
    for indexed in options.grammar:
        grammar = Grammar(grammar_dir=indexed)
        
        for transducer in grammar:
            print transducer
    print

if options.straight or options.invert:
    
    straight = "false"
    invert = "false"

    if options.straight:
        straight = "true"
    if options.invert:
        invert = "true"
    print "# straight glue rule: %s ||| %s %s ||| %s %s" %(options.goal,
                                                           non_terminal(options.goal, 1), non_terminal(options.glue, 2),
                                                           non_terminal(options.goal, 1), non_terminal(options.glue, 2))
    print "# inverted glue rule: %s ||| %s %s ||| %s %s" %(options.goal,
                                                           non_terminal(options.goal, 1), non_terminal(options.glue, 2),
                                                           non_terminal(options.glue, 2), non_terminal(options.goal, 1))
    print "grammar = glue:goal=%s,non-terminal=%s,straight=%s,invert=%s" %(options.goal, options.glue, straight, invert)
    print

if options.insertion or options.deletion:
    if options.insertion:
        print "# insertion grammar %s ||| terminal ||| terminal" %(options.glue)
        print "grammar = insertion:non-terminal=%s" %(options.glue)
    if options.deletion:
        print "# deletion grammar %s ||| terminal ||| <epsilon>" %(options.glue)
        print "grammar = deletion:non-terminal=%s" %(options.glue)
    print

if options.tree_grammar:
    print "# tree-grammar"
    for indexed in options.tree_grammar:
        grammar = TreeGrammar(grammar_dir=indexed)
    
        for transducer in grammar:
            print transducer
    print

if options.fallback:
    print "tree-grammar = fallback:non-terminal=%s" %(options.glue)
    print

### feature-functions



### operations


print "# operations"
if options.scfg:
    print "operation = compose-cky"
elif options.phrase:
    print "operation = compose-phrase"
elif options.tree:
    print "operation = compose-tree"
elif options.tree_cky:
    print "operation = compose-tree-cky"
else:
    raise ValueError, "no operations? --{scfg,phrase,tree,tree-cky}"
print "operation = push-bos-eos"
print "operation = apply:prune=true,size=200,${weights}"
print "operation = remove-bos-eos:forest=true"

if options.file:
    print "operation = output:file=${file},kbest=${kbest},${weights}"
elif options.directory:
    print "operation = output:directory=${directory},kbest=${kbest},${weights}"
else:
    raise ValueError, "no output? --{file,directory}"
