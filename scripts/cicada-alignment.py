#!/usr/bin/env python
#
#  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
#
### a wrapper script (similar to phrase-extract in moses)
### we support only "extraction" meaning only step 5 and 6
### TODO: use argparse for command-lines...?

import threading
import multiprocessing

import time
import sys
import os, os.path
import string
import re
import subprocess

from optparse import OptionParser, make_option

opt_parser = OptionParser(
    option_list=[
	
    # output directory/filename prefix
    make_option("--root-dir", default="", action="store", type="string",
                metavar="DIRECTORY", help="root directory for outputs"),
    make_option("--corpus-dir", default="", action="store", type="string",
                metavar="PREFIX", help="corpus directory (default: ${root_dir}/corpus)"),
    make_option("--giza-f2e", default="", action="store", type="string",
                metavar="DIRECTORY", help="giza directory for P(f|e) (default: ${root_dir}/giza.${f}-${e})"),
    make_option("--giza-e2f", default="", action="store", type="string",
                metavar="DIRECTORY", help="giza directory for P(e|f) (default: ${root_dir}/giza.${e}-${f})"),
    make_option("--model-dir", default="", action="store", type="string",
                metavar="DIRECTORY", help="model directory (default: ${root_dir}/model)"),
    make_option("--alignment-dir", default="", action="store", type="string",
                metavar="DIRECTORY", help="alignment directory (default: ${model_dir})"),
    make_option("--lexical-dir", default="", action="store", type="string",
                metavar="DIRECTORY", help="lexical transltion table directory (default: ${model_dir)"),
    
    ### source/target flags
    make_option("--f", default="F", action="store", type="string",
                metavar="SUFFIX", help="source (or 'French')  language suffix for training corpus"),
    make_option("--e", default="E", action="store", type="string",
                metavar="SUFFIX", help="target (or 'English') language suffix for training corpus"),
    make_option("--a", default="A", action="store", type="string",
                metavar="SUFFIX", help="source-to-target alignment suffix for training corpus"),
    ### span...
    make_option("--sf", default="SF", action="store", type="string",
                metavar="SUFFIX", help="source (or 'French')  span suffix for training corpus"),
    make_option("--se", default="SE", action="store", type="string",
                metavar="SUFFIX", help="target (or 'English') span suffix for training corpus"),
    ### forest!
    make_option("--ff", default="SF", action="store", type="string",
                metavar="SUFFIX", help="source (or 'French')  forest suffix for training corpus"),
    make_option("--fe", default="SE", action="store", type="string",
                metavar="SUFFIX", help="target (or 'English') forest suffix for training corpus"),
    
    # data prefix
    make_option("--corpus", default="corpus", action="store", type="string",
                help="bilingual trainging corpus"),

    # alignment method
    make_option("--alignment", default="grow-diag-final-and", action="store", type="string",
                help="alignment methods (default: grow-diag-final-and)"),
    
    # steps
    make_option("--first-step", default=1, action="store", type="int", metavar='STEP', help="first step (default: 1)"),
    make_option("--last-step",  default=3, action="store", type="int", metavar='STEP', help="last step  (default: 3)"),
    
    ## iteratin
    make_option("--iteration-cluster", default=50, action="store", type="int", metavar='ITERATION', help="word cluter iterations (default: 50)"),
    make_option("--iteration-model1",  default=5,  action="store", type="int", metavar='ITERATION', help="Model1 iteratins (default: 5)"),
    make_option("--iteration-hmm",     default=5,  action="store", type="int", metavar='ITERATION', help="HMM iteratins    (default: 5)"),
    
    ## # of clusters
    make_option("--cluster",     default=50, action="store", type="int", metavar='CLUSTER', help="# of clusters (default: 50)"),
    
    ## training parameters
    make_option("--symmetric",   default=None, action="store_true", help="symmetric training"),
    make_option("--posterior",   default=None, action="store_true", help="posterior constrained training"),
    make_option("--variational", default=None, action="store_true", help="variational Bayes estimates"),
    
    ## options for lexicon model training
    make_option("--p0",              default=0.01, action="store", type="float", metavar='P0',    help="parameter for NULL alignment (default: 0.01)"),
    make_option("--prior-lexicon",   default=0.01, action="store", type="float", metavar="PRIOR", help="lexicon model prior (default: 0.01)"),
    make_option("--prior-alignment", default=0.01, action="store", type="float", metavar="PRIOR", help="alignment model prior (default: 0.01)"),
    make_option("--smooth-lexicon",   default=1e-20, action="store", type="float", metavar="SMOOTH", help="smoothing for lexicon model (default: 1e-20)"),
    make_option("--smooth-alignment", default=1e-20, action="store", type="float", metavar="SMOOTH", help="smoothing for alignment model (default: 1e-20)"),

    # CICADA Toolkit directory
    make_option("--cicada-dir", default="", action="store", type="string",
                metavar="DIRECTORY", help="cicada directory"),
    
    make_option("--threads", default=2, action="store", type="int",
                help="# of thrads for thread-based parallel processing"),

    ## debug messages
    make_option("--debug", default=0, action="store", type="int"),
    ])


### dump to stderr
stdout = sys.stdout
sys.stdout = sys.stderr

def run_command(command):
    fp = os.popen(command)
    while 1:
        data = fp.read(1)
        if not data: break
        stdout.write(data)

def compressed_file(file):
    if not file:
        return file
    if os.path.exists(file):
        return file
    if os.path.exists(file+'.gz'):
	return file+'.gz'
    if os.path.exists(file+'.bz2'):
	return file+'.bz2'
    (base, ext) = os.path.splitext(file)
    if ext == '.gz' or ext == '.bz2':
	if os.path.exists(base):
	    return base
    return file


class CICADA:
    def __init__(self, dir=""):

	self.dir = dir	
	if not dir: return
	
	if not os.path.exists(self.dir):
	    raise ValueError, self.dir + " does not exist"
	
	self.dir = os.path.realpath(self.dir)
        
	self.bindirs = []
	for dir in ('bin', 'progs', 'scripts'): 
	    bindir = os.path.join(self.dir, dir)
	    if os.path.exists(bindir) and os.path.isdir(bindir):
		self.bindirs.append(bindir)
        self.bindirs.append(self.dir)
	
        for binprog in ('cicada_cluster_word',
                        ## step 1
                        'cicada_lexicon_model1',
                        'cicada_lexicon_hmm',
                        'cicada_lexicon_dice',
                        ## step2
                        'cicada_alignment', 
                        ## step 3
                        ):
	    
	    for bindir in self.bindirs:
		prog = os.path.join(bindir, binprog)
		if os.path.exists(prog):
		    setattr(self, binprog, prog)
		    break
	    if not hasattr(self, binprog):
		raise ValueError, binprog + ' does not exist'
        
class Corpus:

    def __init__(self, corpus_dir="", corpus="", f="", e="", a="", sf="", se="", ff="", fe=""):

        self.source_tag = f
        self.target_tag = e
        self.alignment_tag = a

        self.source_span_tag = sf
        self.target_span_tag = se

        self.source_forest_tag = ff
        self.target_forest_tag = fe
        
        self.source = compressed_file(corpus+'.'+f)
        self.target = compressed_file(corpus+'.'+e)
        self.alignment = compressed_file(corpus+'.'+a)
        
        self.source_span = compressed_file(corpus+'.'+sf)
        self.target_span = compressed_file(corpus+'.'+se)
        
        self.source_forest = compressed_file(corpus+'.'+ff)
        self.target_forest = compressed_file(corpus+'.'+fe)

        self.corpus_dir = corpus_dir

class Cluster:

    def __init__(self, cicada=None, corpus="", name="", cluster=64, iteration=64, threads=8, debug=0):
        
        self.cicada  = cicada
        
        command = cicada.cicada_cluster_word
        
        command += " --input \"%s\""  %(corpus)
        command += " --output \"%s\"" %(name)
        command += " --cluster %d" %(cluster)
        command += " --iteration %d" %(iteration)
        command += " --threads %d" %(threads)
        
        if debug:
            command += " --debug=%d" %(debug)
        else:
            command += " --debug"
            
        self.command = command
        self.cluster = name
    
    def run(self):
        run_command(self.command)

class Prepare:
    
    def __init__(self, cicada=None, corpus=None, cluster=64, iteration=64, threads=8, debug=0):
        
        if not os.path.exists(corpus.corpus_dir):
            os.makedirs(corpus.corpus_dir)

        self.corpus = corpus
        self.source = Cluster(cicada=cicada,
                              corpus=corpus.source,
                              name=os.path.join(corpus.corpus_dir, corpus.source_tag+'.vcb.classes'),
                              cluster=cluster,
                              iteration=iteration,
                              threads=threads,
                              debug=debug)
        
        self.target = Cluster(cicada=cicada,
                              corpus=corpus.target,
                              name=os.path.join(corpus.corpus_dir, corpus.target_tag+'.vcb.classes'),
                              cluster=cluster,
                              iteration=iteration,
                              threads=threads,
                              debug=debug)
    
    def run(self):
        self.source.run()
        self.target.run()

class Giza:

    def __init__(self,
                 cicada=None,
                 corpus=None,
                 cluster=None,
                 dir_source_target="",
                 dir_target_source="",
                 prefix_source_target="",
                 prefix_target_source="",
                 iteration_model1=5,
                 iteration_hmm=5,
                 prior_lexicon=0.1,
                 prior_alignment=0.1,
                 smooth_lexicon=1e-20,
                 smooth_alignment=1e-20,
                 p0=1e-4,
                 symmetric=None,
                 posterior=None,
                 variational=None,
                 threads=8,
                 debug=0):

        if not os.path.exists(dir_source_target):
            os.makedirs(dir_source_target)
        if not os.path.exists(dir_target_source):
            os.makedirs(dir_target_source)
        
        command = ""
        
        if iteration_hmm > 0:
            command = cicada.cicada_lexicon_hmm
        elif iteration_model1 > 0:
            command = cicada.cicada_lexicon_model1
        else:
            raise ValueError, "invalid model iterations"
        
        command += " --source \"%s\"" %(corpus.source)
        command += " --target \"%s\"" %(corpus.target)

        if os.path.exists(corpus.alignment):
            command += " --alignment \"%s\"" %(corpus.alignment)
            

        if iteration_hmm > 0:
            command += " --classes-source \"%s\"" %(compressed_file(cluster.source.cluster))
            command += " --classes-target \"%s\"" %(compressed_file(cluster.target.cluster))
        
        if iteration_hmm > 0:
            self.alignment_source_target = os.path.join(dir_source_target, prefix_source_target + '.alignment.final.gz')
            self.alignment_target_source = os.path.join(dir_target_source, prefix_target_source + '.alignment.final.gz')

            command += " --output-alignment-source-target \"%s\"" %(self.alignment_source_target)
            command += " --output-alignment-target-source \"%s\"" %(self.alignment_target_source)
            
        self.lexicon_source_target = os.path.join(dir_source_target, prefix_source_target + '.lexicon.final.gz')
        self.lexicon_target_source = os.path.join(dir_target_source, prefix_target_source + '.lexicon.final.gz')
        
        command += " --output-lexicon-source-target \"%s\"" %(self.lexicon_source_target)
        command += " --output-lexicon-target-source \"%s\"" %(self.lexicon_target_source)

        self.viterbi_source_target = os.path.join(dir_source_target, prefix_source_target + '.A3.final.gz')
        self.viterbi_target_source = os.path.join(dir_target_source, prefix_target_source + '.A3.final.gz')
        
        command += " --viterbi-source-target \"%s\"" %(self.viterbi_source_target)
        command += " --viterbi-target-source \"%s\"" %(self.viterbi_target_source)

        if iteration_hmm > 0:
            command += " --iteration-hmm %d" %(iteration_hmm)
            command += " --iteration-model1 %d" %(iteration_model1)
        else:
            command += " --iteration %d" %(iteration_model1)

        if symmetric:
            command += " --symmetric"
        if posterior:
            command += " --posterior"
        if variational:
            command += " --variational-bayes"
        
        self.p0 = p0
        self.prior_lexicon = prior_lexicon
        self.prior_alignment = prior_alignment
        self.smooth_lexicon = smooth_lexicon
        self.smooth_alignment = smooth_alignment
        
        command += " --p0 %.20g" %(p0)
        command += " --prior-lexicon %.20g"   %(prior_lexicon)
        command += " --smooth-lexicon %.20g"   %(smooth_lexicon)
        if iteration_hmm > 0:
            command += " --prior-alignment %.20g" %(prior_alignment)
            command += " --smooth-alignment %.20g" %(smooth_alignment)
                    
        command += " --threads %d" %(threads)

        if debug:
            command += " --debug=%d" %(debug)
        else:
            command += " --debug"
        
        self.command = command

    def run(self):
        run_command(self.command)

class AlignmentHeuristic:

    def __init__(self,
                 cicada=None,
                 corpus=None,
                 cluster=None,
                 giza=None,
                 alignment_dir="",
                 alignment="grow-diag-final-and",
                 threads=8,
                 debug=0):

        if not os.path.exists(alignment_dir):
            os.makedirs(alignment_dir)
        
        command = cicada.cicada_alignment
        command += " --source-target \"%s\"" %(compressed_file(giza.viterbi_source_target))
        command += " --target-source \"%s\"" %(compressed_file(giza.viterbi_target_source))
        
        self.alignment = os.path.join(alignment_dir, "aligned." + alignment)

        command += " --output \"%s\"" %(self.alignment)
        
        if 'f2e' in alignment:
            command += " --f2e"
        if 'e2f' in alignment:
            command += " --e2f"
        
        if 'union' in alignment:
            command += " --union"
        elif 'intersection' in alignment:
            command += " --intersection"
        elif 'grow' in alignment:
            command += " --grow"
            if 'diag' in alignment:
                command += ' --diag'
            if 'final' in alignment:
                if 'final-and' in alignment:
                    command += " --final-and"
                else:
                    command += " --final"

        if 'itg' in alignment:
            command += " --itg"
        if 'max-match' in alignment:
            command += " --max-match"
        
        command += " --threads %d" %(threads)

        if debug:
            command += " --debug=%d" %(debug)
        else:
            command += " --debug"
        
        self.command = command

    def run(self):
        run_command(self.command)


class AlignmentPosterior:

    def __init__(self,
                 cicada=None,
                 corpus=None,
                 cluster=None,
                 giza=None,
                 alignment_dir="",
                 alignment="grow-diag-final-and",
                 threads=8,
                 debug=0):

        if not os.path.exists(alignment_dir):
            os.makedirs(alignment_dir)

        learn_hmm = None
        if hasattr(giza, 'alignment_source_target'):
            learn_hmm = 1
            
        
        command = cicada.cicada_lexicon_model1
        if learn_hmm:
            command = cicada.cicada_lexicon_hmm
           
        command += " --source \"%s\"" %(corpus.source)
        command += " --target \"%s\"" %(corpus.target)

        if learn_hmm:
            command += " --classes-source \"%s\"" %(compressed_file(cluster.source.cluster))
            command += " --classes-target \"%s\"" %(compressed_file(cluster.target.cluster))

        command += " --lexicon-source-target \"%s\"" %(compressed_file(giza.lexicon_source_target))
        command += " --lexicon-target-source \"%s\"" %(compressed_file(giza.lexicon_target_source))

        if learn_hmm:
            command += " --alignment-source-target \"%s\"" %(compressed_file(giza.alignment_source_target))
            command += " --alignment-target-source \"%s\"" %(compressed_file(giza.alignment_target_source))
        
        self.alignment = os.path.join(alignment_dir, "aligned." + alignment)
        
        command += " --viterbi-source-target \"%s\"" %(self.alignment)
        # we do not specify an alternative alignment!
        # command += " --viterbi-target-source /dev/null"

        if learn_hmm:
            command += " --iteration-model1 0"
            command += " --iteration-hmm 0"
        else:
            command += " --iteration 0"
            
        command += " --p0 %.20g" %(giza.p0)
        command += " --prior-lexicon %.20g"  %(giza.prior_lexicon)
        command += " --smooth-lexicon %.20g" %(giza.smooth_lexicon)
        if learn_hmm:
            command += " --prior-alignment %.20g"  %(giza.prior_alignment)
            command += " --smooth-alignment %.20g" %(giza.smooth_alignment)

        if 'itg' in alignment:
            command += " --itg"
        if 'max-match' in alignment:
            command += " --max-match"

        ## dump in moses mode
        command += " --moses"
        
        command += " --threads %d" %(threads)

        if debug:
            command += " --debug=%d" %(debug)
        else:
            command += " --debug"
        
        self.command = command

    def run(self):
        run_command(self.command)

(options, args) = opt_parser.parse_args()

if options.root_dir:
    if not os.path.exists(options.root_dir):
	os.makedirs(options.root_dir)

if not options.corpus_dir:
    options.corpus_dir = os.path.join(options.root_dir, "corpus")
if not options.giza_f2e:
    options.giza_f2e = os.path.join(options.root_dir, "giza.%s-%s" %(options.f, options.e))
if not options.giza_e2f:
    options.giza_e2f = os.path.join(options.root_dir, "giza.%s-%s" %(options.e, options.f))
if not options.model_dir:
    options.model_dir = os.path.join(options.root_dir, "model")
if not options.lexical_dir:
    options.lexical_dir = options.model_dir
if not options.alignment_dir:
    options.alignment_dir = options.model_dir

cicada = CICADA(options.cicada_dir)

corpus = Corpus(corpus_dir=options.corpus_dir,
                corpus=options.corpus,
                f=options.f,
                e=options.e,
                a=options.a,
                sf=options.sf,
                se=options.se,
                ff=options.ff,
                fe=options.fe)

prepare = Prepare(cicada=cicada,
                  corpus=corpus,
                  cluster=options.cluster,
                  iteration=options.iteration_cluster,
                  threads=options.threads,
                  debug=options.debug)

if options.first_step <= 1 and options.last_step >= 1:
    print "(1) preparing corpus started  @", time.ctime()
    prepare.run()
    print "(1) preparing corpus finished @", time.ctime()

giza = Giza(cicada=cicada,
            corpus=corpus,
            cluster=prepare,
            dir_source_target=options.giza_e2f,
            dir_target_source=options.giza_f2e,
            prefix_source_target=corpus.target_tag+'-'+corpus.source_tag,
            prefix_target_source=corpus.source_tag+'-'+corpus.target_tag,
            iteration_model1=options.iteration_model1,
            iteration_hmm=options.iteration_hmm,
            prior_lexicon=options.prior_lexicon,
            prior_alignment=options.prior_alignment,
            smooth_lexicon=options.smooth_lexicon,
            smooth_alignment=options.smooth_alignment,
            p0=options.p0,
            symmetric=options.symmetric,
            posterior=options.posterior,
            variational=options.variational,
            threads=options.threads,
            debug=options.debug)

## run giza++ in two directions
if options.first_step <= 2 and options.last_step >= 2:
    print "(2) running giza started  @", time.ctime()
    giza.run()
    print "(2) running giza finished @", time.ctime()

alignment=None
if "posterior" in options.alignment:
    alignment = AlignmentPosterior(cicada=cicada,
                                   corpus=corpus,
                                   cluster=prepare,
                                   giza=giza,
                                   alignment_dir=options.alignment_dir,
                                   alignment=options.alignment,
                                   threads=options.threads,
                                   debug=options.debug)
else:
    alignment = AlignmentHeuristic(cicada=cicada,
                                   corpus=corpus,
                                   cluster=prepare,
                                   giza=giza,
                                   alignment_dir=options.alignment_dir,
                                   alignment=options.alignment,
                                   threads=options.threads,
                                   debug=options.debug)

if options.first_step <= 3 and options.last_step >= 3:
    print "(3) generate word alignment started  @", time.ctime()
    alignment.run()
    print "(3) generate word alignment finished @", time.ctime()


