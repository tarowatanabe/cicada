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
from optparse import OptionParser, make_option

opt_parser = OptionParser(
    option_list=[
            
    ### grammars
    make_option("--grammar",      default="", action="store", type="string", help="indexed grammar"),
    make_option("--tree-grammar", default="", action="store", type="string", help="indexed tree-grammar"),
    
    ### only affect grammars
    make_option("--max-span", default=15, action="store", type="int",
                metavar="LENGTH", help="maximum span size (default: 15)"),
    
    ### glues
    make_option("--goal", default="[s]", action="store", type="string", help="goal non-terminal (default: [s])"),
    make_option("--glue", default="[x]", action="store", type="string", help="non-terminal for glue rules (default: [x])"),

    ## operations...
    
    # cicada composition
    make_option("--phrase",   default=None, action="store_true", help="phrase-based grammar"),
    make_option("--scfg",     default=None, action="store_true", help="SCFG"),
    make_option("--tree",     default=None, action="store_true", help="tree-to-string"),
    make_option("--tree-cky", default=None, action="store_true", help="string-to-{string,tree}"),
    
    ### feature application strategy
    make_option("--apply-cube", default=None, action="store_true", help="feature application strategy (cube_pruning)"),
    make_option("--apply-grow", default=None, action="store_true", help="feature application strategy (cube_growing)"),
    
    ### some recommended configuration
    make_option("--push-bos-eos", default=None, action="store_true", help="sentence bondary marker pushing"),
    
    ### outputs
    make_option("--output-file",      default=None, action="store_true", help="output in a file"),
    make_option("--output-directory", default=None, action="store_true", help="output in a directory"),
    
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
        
        for line in open(paht_file):
            
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
        
        for line in open(paht_file):
            
            name = line.strip()
            if not name: continue
            
            path = os.path.join(grammar_dir, name)
            if not os.path.exists(path):
                raise ValueError, "no path to scores: %s" %(path)
            
            self.append("tree-grammar = " + path)

(options, args) = opt_parser.parse_args()
