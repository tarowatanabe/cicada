#!/usr/bin/env python
#
#  Copyright(C) 2011-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
#

###
### a cicada.config generator
###

###
### We will generate common cicada.config from indexed grammar
### (SCFG and Tree-to-string)
### 

import sys
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
                metavar="LENGTH", help="maximum span size (default: %default)"),
    
    ### glues
    make_option("--goal", default="[s]", action="store", type="string", help="goal non-terminal (default: %default)"),
    make_option("--glue", default="[x]", action="store", type="string", help="non-terminal for glue rules (default: %default)"),
    
    make_option("--straight", default=None, action="store_true", help="straight gulue rule for string-to-string"),
    make_option("--invert",   default=None, action="store_true", help="invert gulue rule for string-to-string"),
    make_option("--tree-straight", default=None, action="store_true", help="straight gulue rule for string-to-tree"),
    make_option("--tree-invert",   default=None, action="store_true", help="invert gulue rule for string-to-tree"),
    
    make_option("--insertion", default=None, action="store_true", help="insertion grammar for string-to-{string,tree}"),
    make_option("--deletion",  default=None, action="store_true", help="deletion grammar for string-to-{string,tree}"),
    make_option("--fallback-glue", default=None, action="store_true", help="fallback tree-grammar with glue non-terminal"),
    make_option("--fallback",      default=None, action="store_true", help="fallback tree-grammar"),
    
    ### feature functions
    make_option("--feature-ngram", default=[], action="append", type="string", help="ngram feature"),
    make_option("--feature-kenlm", default=[], action="append", type="string", help="kenlm feature"),
    make_option("--feature-lexicon", default="", action="store", type="string", help="lexicon feature"),
    
    ## operations...
    
    # cicada composition
    make_option("--phrase",   default=None, action="store_true", help="phrase-based grammar"),
    make_option("--scfg",     default=None, action="store_true", help="SCFG"),
    make_option("--ghkm",     default=None, action="store_true", help="GHKM (or tree-to-{string,tree})"),
    make_option("--tree",     default=None, action="store_true", help="tree-to-{string,tree}"),
    make_option("--tree-cky", default=None, action="store_true", help="string-to-{string,tree}"),

    # beam size
    make_option("--beam", default=200, type=int, action="store", help="beam size (default: %default)"),
    
    ## debug messages
    make_option("--debug", default=0, action="store", type="int"),
    ])

def escape_path(path):
    if ':' in path:
        return '"' + path + '"'
    else:
        return path

class Grammar(UserList.UserList):
    
    def __init__(self, grammar_dir="", max_span=15):

        self.max_span = max_span
        
        UserList.UserList.__init__(self)
        
        path_files = os.path.join(grammar_dir, 'files');
        
        if not os.path.exists(path_files):
            self.append(grammar_dir)
        else:
            for line in open(path_files):
                name = line.strip()
                if not name: continue
                
                self.append(os.path.join(grammar_dir, name))

    def append(self, path):
        if not os.path.exists(path):
            raise ValueError, "no path to grammar: %s" %(path)
        
        UserList.UserList.append(self, "grammar = " + escape_path(path) + ":max-span=%s" %(self.max_span))

class TreeGrammar(UserList.UserList):
    
    def __init__(self, grammar_dir="", max_span=15):

        self.max_span = max_span
        
        UserList.UserList.__init__(self)
        
        path_files = os.path.join(grammar_dir, 'files');
        
        if not os.path.exists(path_files):
            self.append(grammar_dir)
        else:
            for line in open(path_files):
                name = line.strip()
                if not name: continue
            
                self.append(os.path.join(grammar_dir, name))

    def append(self, path):
        if not os.path.exists(path):
            raise ValueError, "no path to grammar: %s" %(path)
        
        UserList.UserList.append(self, "tree-grammar = " + escape_path(path) + ":max-span=%s" %(self.max_span))

def non_terminal(x, index=0):
    if not x:
        raise ValueError, "no non-terminal? %s"  %(x)
    if x[0] != '[' or x[-1] != ']':
        raise ValueError, "invalid non-terminal? %s" %(x)
    if index == 0:
        return x
    else:
        return '[' + x[1:-1] + ',' + str(index) + ']'

if __name__ == '__main__':
    
    opt_parser.disable_interspersed_args()
    
    (options, args) = opt_parser.parse_args()

    # check compositions...
    num_composition = 0
    if options.phrase:
        num_composition += 1
    if options.scfg:
        num_composition += 1
    if options.ghkm:
        num_composition += 1
    if options.tree:
        num_composition += 1
    if options.tree_cky:
        num_composition += 1

    if num_composition > 1:
        raise ValueError, "only one of --{phrase,scfg,ghkm,tree,tree-cky} can be specified"

    ### grammars
    options.goal = non_terminal(options.goal)
    options.glue = non_terminal(options.glue)
    
    print "#"
    print "# config file for cicada generated by %s" %(' '.join(sys.argv))
    print "#"
    print

    print "# goal for parsing/composition"
    print "goal = %s" %(options.goal)
    print

    if options.grammar:
        print "#"
        print "# grammar. For details, see \"cicada --grammar-list\""
        print "#"
        for indexed in options.grammar:
            grammar = Grammar(grammar_dir=indexed, max_span=options.max_span)
            
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
            
        print "# glue rules for string-to-string"
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
            print "# insertion grammar %s ||| terminal ||| terminal for string-to-{string,tree}" %(options.glue)
            print "grammar = insertion:non-terminal=%s" %(options.glue)
        if options.deletion:
            print "# deletion grammar %s ||| terminal ||| <epsilon> for string-to-{string,tree}" %(options.glue)
            print "grammar = deletion:non-terminal=%s" %(options.glue)
        print

    if options.tree_grammar:
        print "#"
        print "# tree-grammar. For details see \"cicada --tree-grammar-list\""
        print "#"
        for indexed in options.tree_grammar:
            grammar = TreeGrammar(grammar_dir=indexed, max_span=options.max_span)
    
            for transducer in grammar:
                print transducer
        print

    if options.tree_straight or options.tree_invert:
    
        straight = "false"
        invert = "false"

        if options.tree_straight:
            straight = "true"
        if options.tree_invert:
            invert = "true"
            
        print "# glue rules for string-to-tree"
        print "# straight glue rule: %s(%s %s) ||| %s(%s %s)" %(options.goal,
                                                                non_terminal(options.goal, 1), non_terminal(options.glue, 2),
                                                                options.goal,
                                                                non_terminal(options.goal, 1), non_terminal(options.glue, 2))
        print "# inverted glue rule: %s(%s %s) ||| %s(%s %s)" %(options.goal,
                                                                non_terminal(options.goal, 1), non_terminal(options.glue, 2),
                                                                options.goal,
                                                                non_terminal(options.glue, 2), non_terminal(options.goal, 1))
        print "tree-grammar = glue:goal-source=%s,goal-target=%s,non-terminal-source=%s,non-terminal-target=%s,straight=%s,invert=%s" %(options.goal, options.goal, options.glue, options.glue, straight, invert)
        print

    if options.fallback and options.fallback_glue:
        raise ValueError, "specified both of fallback glue and fallback"

    if options.fallback:
        print "# fallback grammar with glue non-terminal for tree-to-{string,tree}"
        print "# tree-grammar = fallback:goal=%s, non-terminal=%s" %(options.goal, options.glue)
        print "# or, use this fallback grammar which will not replace source-side non-terminals with %s" %(options.glue)
        print "tree-grammar = fallback"
        print
    elif options.fallback_glue:
        print "# fallback grammar with glue non-terminal for tree-to-{string,tree}"
        print "tree-grammar = fallback:goal=%s, non-terminal=%s" %(options.goal, options.glue)
        print "# or, use this fallback grammar which will not replace source-side non-terminals with %s" %(options.glue)
        print "# tree-grammar = fallback"
        print


    ### feature-functions

    print "#"
    print "# feature functions. For details see \"cicada --feature-function-list\""
    print "#"

    if options.feature_ngram:
        print "# ngram feature. If you have multiple ngrams, you should modify name"
        print "# no-boe-eos=true implies that the forest is explicitly annotated with <s> and </s>."
        for ngram in options.feature_ngram:
            print "feature-function = ngram: name=ngram, no-bos-eos=true, file=%s" %(ngram)
        print

    if options.feature_kenlm:
        print "# kenlm feature. If you have multiple ngrams, you should modify name"
        print "# no-boe-eos=true implies that the forest is explicitly annotated with <s> and </s>."
        for ngram in options.feature_kenlm:
            print "feature-function = kenlm: name=kenlm, no-bos-eos=true, file=%s" %(ngram)
        print
    
    
    if options.feature_lexicon:
        print "# lexicon feature computes P(target-sentence | source-sentence) based on model1/viterbi/noisy-or"
        print "feature-function = lexicon: lexicon=%s" %(options.feature_lexicon)
        print

    print "feature-function = word-penalty"
    print "feature-function = rule-penalty"
    
    if options.tree or options.ghkm or options.tree_cky:
        print "feature-function = glue-tree-penalty"
    else:
        print "# feature-function = glue-tree-penalty"
    
    print "# feature-function = arity-penalty"
    print "# feature-function = non-latin-penalty"
    print "# feature-function = rule-shape"
    print
    
    ### inputs
    print "#"
    print "# inputs. We support: input-{id,bitext,sentence,lattice,forest,span,alignemnt,dependency}"
    print "#"

    if options.scfg:
        print "input-sentence = true"
    elif options.phrase:
        print "input-sentence = true"
    elif options.tree or options.ghkm:
        print "input-forest = true"
    elif options.tree_cky:
        print "input-sentence = true"
    print

    ### operations

    print "#"
    print "# operations. For details, see \"cicada --operation-list\""
    print "#"
    if options.scfg:
        print "# SCFG translation"
        print "operation = compose-cky"
        print "# Alternatively, you can use parse-cky to perform beam search during composition"
        print "# Here, we use the features found only in the grammar"
        print "# operation = parse-cky:size=%d,${weights}" %(options.beam)
    elif options.phrase:
        print "# phrase translation"
        print "operation = compose-phrase"
    elif options.tree or options.ghkm:
        print "# tree-to-{string,tree} translation"
        print "operation = compose-tree"
        print "# Alternatively, you can use parse-tree to perform beam search during composition"
        print "# Here, we use the features found only in the grammar"
        print "# operation = parse-tree:size=%d,${weights}" %(options.beam)
    elif options.tree_cky:
        print "# string-to-{string,tree} translation"
        print "operation = compose-tree-cky"
        print "# Alternatively, you can use parse-tree-cky to perform beam search during composition"
        print "# Here, we use the features found only in the grammar"
        print "# operation = parse-tree-cky:size=%d,${weights}" %(options.beam)
    print

    print "# annotate <s> and </s> to the forest, so that we can compute ngram LM scores with no-bos-eos=true"
    print "operation = push-bos-eos"
    print

    print "# cube-pruning"
    print "operation = apply:prune=true,size=%d,${weights}" %(options.beam)
    print
    
    print "# forest-pruning"
    print "# this is an example of pruning by density. Alternatives are: beam, edge or kbest. "
    print "# operation = prune:density=10,${weights}"
    print
    
    print "# remove <s> and </s>"
    print "operation = remove-bos-eos:forest=true"
    print
    
    print "# expand-ngram"
    print "# This is used to expand forest so that collecting ngrams are easier without memorizig states."
    print "# This is required when optimizing via xBLEU over forest."
    print "# operation = expand-ngram:order=4"
    print

    print "# non-MBR decoding, and output forest or kbests"
    print "# kbest=0 implies forest output, otherwise, kbest outputs"
    print "operation = output:${file},kbest=${kbest},unique=true,${weights}"
    print "# for a simple, forest output"
    print "# operation = output:${file},forest=true"
    print

    print "# MBR decoding"
    print
    print "# Optionally, prune forest"
    print "# operation = prune:density=10,${weights}"
    print
    print "# First, collect expected ngrams"
    print "# operation = expected-ngram:${weights},scale=1.2,order=4"
    print
    print "# Second, compute expected-BLEU"
    print "# Here, we rescore by bleu-expected"
    print "# The weight file, \"weights.mbr\" should contains a single line: \"bleu-expected 1\""
    print "# so that the bleu-expected feature is used during cube-pruning."
    print "# operation = apply:prune=true,size=%d,feature=\"bleu-expected:order=4\",weights=weights.mbr" %(options.beam)
    print
    print "# Third, output!"
    print "# Here, we use the same parameter employed in the previous cube-pruning"
    print "# operation = output:${file},kbest=${kbest},unique=true,weights=weights.mbr"
    print 
