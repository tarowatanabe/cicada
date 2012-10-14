#!/usr/bin/env python
#
#  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
#
###
### a wrapper script for indexing models, such as grammar, rule-grammar, word-clusters and global-lexicon
### Currently, we support only grammar and tree-grammar
###

import threading
import multiprocessing

import time
import sys
import os, os.path
import string
import re
import subprocess
import shutil

import UserList
import UserString
import cStringIO

from optparse import OptionParser, make_option

opt_parser = OptionParser(
    option_list=[
    # output directory/filename prefix
    make_option("--root-dir", default="", action="store", type="string",
                metavar="DIRECTORY", help="root directory for outputs"),
    make_option("--model-dir", default="", action="store", type="string",
                metavar="DIRECTORY", help="model directory (default: ${root_dir}/model)"),
    make_option("--lexical-dir", default="", action="store", type="string",
                metavar="DIRECTORY", help="lexical transltion table directory (default: ${model_dir})"),

    make_option("--lexicon-source-target", default="", action="store", type="string",
                metavar="LEXICON", help="lexicon for P(target | source) (default: ${lexical_dir}/lex.f2n)"),
    make_option("--lexicon-target-source", default="", action="store", type="string",
                metavar="LEXICON", help="lexicon for P(source | target) (default: ${lexical_dir}/lex.n2f)"),

    ## smoothing...
    make_option("--prior", default=0.1, action="store", type="float", metavar="PRIOR", help="model prior (default: 0.1)"),
    
    ## feature/attribute names
    make_option("--feature",   default=[], action="append", type="string", help="feature definitions"),
    make_option("--attribute", default=[], action="append", type="string", help="attribute definitions"),
    
    ### options...
    make_option("--phrase", default=None, action="store_true", help="index phrase grammar"),
    make_option("--scfg",   default=None, action="store_true", help="index SCFG grammar"),
    make_option("--ghkm",   default=None, action="store_true", help="index GHKM (tree-to-string) grammar"),
    make_option("--tree",   default=None, action="store_true", help="index tree-to-tree grammar"),
    
    make_option("--cky",    default=None, action="store_true", help="CKY mode indexing for tree-grammar"),
    make_option("--reordering", default=None, action="store_true", help="reordering for phrase grammar"),
    
    ## additional feature functions
    make_option("--feature-root",               default=None, action="store_true", help="generative probability"),
    make_option("--feature-type",               default=None, action="store_true", help="observation probability"),
    make_option("--feature-singleton",          default=None, action="store_true", help="singleton features"),
    make_option("--feature-cross",              default=None, action="store_true", help="cross features"),
    make_option("--feature-unaligned",          default=None, action="store_true", help="unaligned features"),
    make_option("--feature-internal",           default=None, action="store_true", help="internal features"),
    make_option("--feature-height",             default=None, action="store_true", help="height features"),

    make_option("--feature-lexicon",            default=None, action="store_true", help="compute Model1 features"),
    make_option("--feature-model1",             default=None, action="store_true", help="compute Model1 features"),
    make_option("--feature-noisy-or",           default=None, action="store_true", help="compute noisy-or features"),
    make_option("--feature-insertion-deletion", default=None, action="store_true", help="compute insertion/deletion features"),
    
    make_option("--threshold-insertion", default=0.5, action="store", type="float", help="threshold for insertion (default: 0.5)"),
    make_option("--threshold-deletion",  default=0.5, action="store", type="float", help="threshold for deletion (default: 0.5)"),
    
    ## quantize
    make_option("--quantize", default=None, action="store_true", help="perform quantization"),
    
    ## kbest options
    make_option("--kbest", default=0, action="store", type="float",
                metavar="KBEST", help="kbest max count of rules"),
    ## max-malloc
    make_option("--max-malloc", default=4, action="store", type="float",
                metavar="MALLOC", help="maximum memory in GB (default: 4)"),

    # CICADA Toolkit directory
    make_option("--cicada-dir", default="", action="store", type="string",
                metavar="DIRECTORY", help="cicada directory"),
    # MPI Implementation.. if different from standard location...
    make_option("--mpi-dir", default="", action="store", type="string",
                metavar="DIRECTORY", help="MPI directory"),

    make_option("--threads", default=2, action="store", type="int",
                help="# of thrads for thread-based parallel processing"),
    # perform threading or MPI training    
    make_option("--mpi", default=0, action="store", type="int",
                help="# of processes for MPI-based parallel processing. Identical to --np for mpirun"),
    make_option("--mpi-host", default="", action="store", type="string",
                help="list of hosts to run job. Identical to --host for mpirun", metavar="HOSTS"),
    make_option("--mpi-host-file", default="", action="store", type="string",
                help="host list file to run job. Identical to --hostfile for mpirun", metavar="FILE"),
    make_option("--pbs", default=None, action="store_true",
                help="PBS for launching processes"),
    make_option("--pbs-queue", default="ltg", action="store", type="string",
                help="PBS queue for launching processes (default: ltg)", metavar="NAME"),

    ## debug messages
    make_option("--debug", default=0, action="store", type="int"),
    ])


def run_command(command):
    try:
        retcode = subprocess.call(command, shell=True)
        if retcode:
            sys.exit(retcode)
    except:
        raise ValueError, "subprocess.call failed: %s" %(command)

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

class Quoted:
    def __init__(self, arg):
        self.arg = arg
        
    def __str__(self):
        return '"' + str(self.arg) + '"'

class Option:
    def __init__(self, arg, value=None):
        self.arg = arg
        self.value = value

    def __str__(self,):
        option = self.arg
        
        if self.value is not None:
            if isinstance(self.value, int):
                option += " %d" %(self.value)
            elif isinstance(self.value, long):
                option += " %d" %(self.value)
            elif isinstance(self.value, float):
                option += " %.20g" %(self.value)
            else:
                option += " %s" %(str(self.value))
        return option

class Program:
    def __init__(self, *args):
        self.args = list(args[:])

    def __str__(self,):
        return ' '.join(map(str, self.args))
    
    def __iadd__(self, other):
        self.args.append(other)
        return self

class QSUB(multiprocessing.Process):
    def __init__(self, command=""):
        
        multiprocessing.Process.__init__(self)
        self.command = command
        
    def run(self):
        popen = subprocess.Popen("qsub -S /bin/sh", shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        popen.communicate(self.command)
        
class PBS:
    ### rewrite this by threads!
    def __init__(self, queue=""):
        self.queue = queue
        self.qsub = 'qsub'

        self.workers = []

    def __del__(self):
        for worker in self.workers:
            worker.join()
            
    def run(self, command="", threads=1, memory=0.0, name="name", logfile=None):
        pipe = cStringIO.StringIO()
        
        pipe.write("#!/bin/sh\n")
        pipe.write("#PBS -N %s\n" %(name))
        pipe.write("#PBS -e /dev/null\n")
        pipe.write("#PBS -o /dev/null\n")
        pipe.write("#PBS -W block=true\n")
        
        if self.queue:
            pipe.write("#PBS -q %s\n" %(self.queue))
        
        mem = ""
        if memory >= 1.0:
            mem=":mem=%dgb" %(int(memory))
        elif memory >= 0.001:
            mem=":mem=%dmb" %(int(amount * 1000))
        elif memory >= 0.000001:
            mem=":mem=%dkb" %(int(amount * 1000 * 1000))

        pipe.write("#PBS -l select=1:ncpus=%d:mpiprocs=1%s\n" %(threads, mem))
        
        # setup variables
        if os.environ.has_key('TMPDIR_SPEC'):
            pipe.write("export TMPDIR_SPEC=%s\n" %(os.environ['TMPDIR_SPEC']))
        if os.environ.has_key('LD_LIBRARY_PATH'):
            pipe.write("export LD_LIBRARY_PATH=%s\n" %(os.environ['LD_LIBRARY_PATH']))
        if os.environ.has_key('DYLD_LIBRARY_PATH'):
            pipe.write("export DYLD_LIBRARY_PATH=%s\n" %(os.environ['DYLD_LIBRARY_PATH']))
            
        pipe.write("if test \"$PBS_O_WORKDIR\" != \"\"; then\n")
        pipe.write("  cd $PBS_O_WORKDIR\n")
        pipe.write("fi\n")
        
        if logfile:
            pipe.write("%s 2> %s\n" %(command, logfile))
        else:
            pipe.write("%s\n" %(command))
            
        self.workers.append(QSUB(pipe.getvalue()))
        self.workers[-1].start();

class Threads:
    
    def __init__(self, cicada=None, threads=1):
        
        command = "%s" %(cicada.thrsh)
        command += " --threads %d" %(threads)
        
        self.popen = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE)
        self.pipe = self.popen.stdin
        
    def __del__(self):
        self.pipe.close()
        self.popen.wait()

    def run(self, command="", logfile=None):
        if logfile:
            self.pipe.write("%s 2> %s\n" %(command, logfile))
        else:
            self.pipe.write("%s\n" %(command))
        self.pipe.flush()

class MPI:
    
    def __init__(self, cicada=None, dir="", hosts="", hosts_file="", number=0):
        
	self.dir = dir
	self.hosts = hosts
        self.hosts_file = hosts_file
        self.number = number
	
        if self.dir:
            if not os.path.exists(self.dir):
                raise ValueError, self.dir + " does not exist"
            self.dir = os.path.realpath(self.dir)

        if self.hosts_file:
            if not os.path.exists(self.hosts_file):
                raise ValueError, self.hosts_file + " does no exist"
            self.hosts_file = os.path.realpath(hosts_file)

        self.bindir = self.dir
	
        for binprog in ['mpirun']:
            if self.bindir:
                prog = os.path.join(self.bindir, 'bin', binprog)
                if not os.path.exists(prog):
                    prog = os.path.join(self.bindir, binprog)
                    if not os.path.exists(prog):
                        raise ValueError, prog + " does not exist at " + self.bindir
                    
                setattr(self, binprog, prog)
            else:
                setattr(self, binprog, binprog)
        
        command = self.mpirun
        
        if self.number > 0:
            command += ' --np %d' %(self.number)
        if self.hosts:
            command += ' --host %s' %(self.hosts)
        elif self.hosts_file:
            command += ' --hostfile %s' %(self.hosts_file)

        if os.environ.has_key('TMPDIR_SPEC'):
            command += ' -x TMPDIR_SPEC'
        if os.environ.has_key('LD_LIBRARY_PATH'):
            command += ' -x LD_LIBRARY_PATH'
        if os.environ.has_key('DYLD_LIBRARY_PATH'):
            command += ' -x DYLD_LIBRARY_PATH'
        
        command += " %s" %(cicada.mpish)
        
        self.popen = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE)
        self.pipe = self.popen.stdin
        
    def __del__(self):
        self.pipe.close()
        self.popen.wait()
        
    def run(self, command="", logfile=None):
        if logfile:
            self.pipe.write("%s 2> %s\n" %(command, logfile))
        else:
            self.pipe.write("%s\n" %(command))
        self.pipe.flush()


class CICADA:
    def __init__(self, dir=""):
        bindirs = []
        
        if not dir:
            dir = os.path.abspath(os.path.dirname(__file__))
            bindirs.append(dir)
            parent = os.path.dirname(dir)
            if parent:
                dir = parent
        else:
            dir = os.path.realpath(dir)
            if not os.path.exists(dir):
                raise ValueError, dir + " does not exist"
            bindirs.append(dir)
        
	for subdir in ('bin', 'progs', 'scripts'): 
	    bindir = os.path.join(dir, subdir)
	    if os.path.exists(bindir) and os.path.isdir(bindir):
		bindirs.append(bindir)
	
        for binprog in ('cicada_index_grammar',
                        'cicada_index_tree_grammar',
                        ### indexers
                        'cicada_filter_extract',
                        'cicada_filter_extract_phrase',
                        'cicada_filter_extract_scfg',
                        'cicada_filter_extract_ghkm',
                        ### filters
                        'mpish', ### mpi-launcher
                        'thrsh', ### thread-launcher
                        ### launchers
                        ):
	    
	    for bindir in bindirs:
		prog = os.path.join(bindir, binprog)
                
                if not os.path.exists(prog): continue
                if os.path.isdir(prog): continue
                
                setattr(self, binprog, prog)
                break
            
	    if not hasattr(self, binprog):
		raise ValueError, binprog + ' does not exist'

class IndexPhrase:
    def __init__(self, cicada=None, model_dir="", cky=None, reordering=None):
        self.indexer = cicada.cicada_index_grammar
        self.filter  = cicada.cicada_filter_extract_phrase
        self.filter += " --cicada"
        
        if reordering:
            self.filter += " --reordering --bidirectional"

        self.cky = None
        self.grammar = "grammar"
        self.name = "phrase"
        
        self.counts = os.path.join(model_dir, "phrase-counts")
        self.scores = os.path.join(model_dir, "phrase-score")
        self.index  = os.path.join(model_dir, "phrase-index")
        self.base   = model_dir
        
class IndexSCFG:
    def __init__(self, cicada=None, model_dir="", cky=None, reordering=None):
        self.indexer = cicada.cicada_index_grammar
        self.filter  = cicada.cicada_filter_extract_scfg
        self.cky = None
        self.grammar = "grammar"
        self.name = "scfg"

        self.counts = os.path.join(model_dir, "scfg-counts")
        self.scores = os.path.join(model_dir, "scfg-score")
        self.index  = os.path.join(model_dir, "scfg-index")
        self.base   = model_dir

class IndexGHKM:
    def __init__(self, cicada=None, model_dir="", cky=None, reordering=None):
        self.indexer = cicada.cicada_index_tree_grammar
        self.filter  = cicada.cicada_filter_extract_ghkm
        self.cky = cky
        self.grammar = "tree-grammar"
        self.name = "ghkm"

        self.counts = os.path.join(model_dir, "ghkm-counts")
        self.scores = os.path.join(model_dir, "ghkm-score")
        self.index  = os.path.join(model_dir, "ghkm-index")
        self.base   = model_dir

class IndexTree:
    def __init__(self, cicada=None, model_dir="", cky=None, reordering=None):
        self.indexer = cicada.cicada_index_tree_grammar
        self.filter  = cicada.cicada_filter_extract_ghkm
        self.cky = cky
        self.grammar = "tree-grammar"
        self.name = "tree"

        self.counts = os.path.join(model_dir, "tree-counts")
        self.scores = os.path.join(model_dir, "tree-score")
        self.index  = os.path.join(model_dir, "tree-index")
        self.base   = model_dir

## additional features...
class Features:
    def __init__(self,
                 root=None,
                 types=None,
                 singleton=None,
                 cross=None,
                 unaligned=None,
                 internal=None,
                 height=None):
        self.root      = root
        self.types     = types
        self.singleton = singleton
        self.cross     = cross
        self.unaligned = unaligned
        self.internal  = internal
        self.height    = height

        self.options = ""
        
        if root:
            self.options += " --feature-root"
        if types:
            self.options += " --feature-type"
        if singleton:
            self.options += " --feature-singleton"
        if cross:
            self.options += " --feature-cross"
        if unaligned:
            self.options += " --feature-unaligned"
        if internal:
            self.options += " --feature-internal"
        if height:
            self.options += " --feature-height"

class Lexicon:
    def __init__(self,
                 lexicon_source_target="",
                 lexicon_target_source="",
                 lexicon=None,
                 model1=None,
                 noisy_or=None,
                 insertion_deletion=None,
                 threshold_insertion=0.1,
                 threshold_deletion=0.1):
        self.lexicon_source_target = compressed_file(lexicon_source_target)
        self.lexicon_target_source = compressed_file(lexicon_target_source)
        
        self.lexicon  = lexicon
        self.model1   = model1
        self.noisy_or = noisy_or
        self.insertion_deletion = insertion_deletion
        
        self.threshold_insertion = threshold_insertion
        self.threshold_deletion  = threshold_deletion

        self.options = ""
        
        if lexicon or model1 or noisy_or or insertion_deletion:
            self.options += " --lexicon-source-target \"%s\"" %(self.lexicon_source_target)
            self.options += " --lexicon-target-source \"%s\"" %(self.lexicon_target_source)

            if lexicon:
                self.options += " --feature-lexicon"
            if model1:
                self.options += " --feature-model1"
            if noisy_or:
                self.options += " --feature-noisy-or"
            if insertion_deletion:
                self.options += " --feature-insertion-deletion"
                self.options += " --threshold-insertion %.20g" %(self.threshold_insertion)
                self.options += " --threshold-deletion %.20g" %(self.threshold_deletion)


class Index(UserString.UserString):
    def __init__(self,
                 cicada=None,
                 indexer=None,
                 lexicon=None,
                 feats=None,
                 input="",
                 output="",
                 name="",
                 root_joint="",
                 root_source="",
                 root_target="",
                 prior=0.1,
                 kbest=0,
                 quantize=None,
                 features=[],
                 attributes=[]):

        if not input:
            raise ValueError, "no input? %s" %(input)
        if not output:
            raise ValueError, "no output? %s" %(output)
        if not root_joint:
            raise ValueError, "no root count? %s" %(root_joint)
        if not root_source:
            raise ValueError, "no root source? %s" %(root_source)
        if not root_target:
            raise ValueError, "no root target? %s" %(root_target)
        
        self.name    = "index-" + indexer.name
        self.logfile = os.path.join(indexer.base, "index-" + indexer.name + "." + name + ".log")
        
        command = ""

        if kbest > 0:
            
            self.threads = 3

            command = cicada.cicada_filter_extract
            command += " --nbest %d" %(kbest)
            command += " --input \"%s\"" %(input)
            command += " | "
            command += indexer.filter
            command += " --dirichlet-prior %g" %(prior)
            command += " --root-joint \"%s\""  %(root_joint)
            command += " --root-source \"%s\"" %(root_source)
            command += " --root-target \"%s\"" %(root_target)
            if feats:
                command += feats.options
            if lexicon:
                command += lexicon.options
            command += " | "
            command += indexer.indexer
        else:
            
            self.threads = 2

            command = indexer.filter
            command += " --dirichlet-prior %g" %(prior)
            command += " --root-joint \"%s\""  %(root_joint)
            command += " --root-source \"%s\"" %(root_source)
            command += " --root-target \"%s\"" %(root_target)
            if feats:
                command += feats.options
            if lexicon:
                command += lexicon.options
            command += " --input \"%s\"" %(input)
            command += " | "
            command += indexer.indexer
            
        if quantize:
            command += " --quantize"

        input_path='-'
        sep = ':'
        if indexer.cky:
            input_path += sep
            sep = ','
            input_path +='cky=true'
        if features:
            input_path += sep
            sep = ','
            input_path += ','.join(features)
        if attributes:
            input_path += sep
            sep = ','
            input_path += ','.join(attributes)
        
        # add debug flag
        input_path += sep
        sep = ','
        input_path +='debug=1'

        command += " --input %s" %(input_path)
        command += " --output \"%s\"" %(output)

        UserString.UserString.__init__(self, '('+command+')')
        
class Score:
    def __init__(self, input="", output="", name=""):
        self.input = input
        self.output = output
        self.name = name


class Scores(UserList.UserList):
    def __init__(self, indexer=None):

        UserList.UserList.__init__(self)

        prefix = indexer.scores
        output = indexer.index

        path_files = os.path.join(prefix, 'files');
        path_root_joint  = os.path.join(prefix, 'root-joint.gz')
        path_root_source = os.path.join(prefix, 'root-source.gz')
        path_root_target = os.path.join(prefix, 'root-target.gz')
        
        if not os.path.exists(path_files):
            raise ValueError, "no path to files: %s" %(path_files)
        if not os.path.exists(path_root_joint):
            raise ValueError, "no path to root-joint: %s" %(path_root_joint)
        if not os.path.exists(path_root_source):
            raise ValueError, "no path to root-source: %s" %(path_root_source)
        if not os.path.exists(path_root_target):
            raise ValueError, "no path to root-target: %s" %(path_root_target)
        
        self.root_joint  = path_root_joint
        self.root_source = path_root_source
        self.root_target = path_root_target
        
        if output:
            if os.path.exists(output):
                if not os.path.isdir(output):
                    os.remove(output)
                    os.makedirs(output)
            else:
                os.makedirs(output)
        
        for line in open(path_files):
            name = line.strip()
            if not name: continue
            
            path = os.path.join(prefix, name)
            if not os.path.exists(path):
                raise ValueError, "no path to scores: %s" %(path)

            root,stem = os.path.splitext(name)

            self.append(Score(path, os.path.join(output, root + '.bin'), root))
            
if __name__ == '__main__':
    (options, args) = opt_parser.parse_args()

    ### dump to stderr
    stdout = sys.stdout
    sys.stdout = sys.stderr

    ### setup defaults!
    
    if options.root_dir:
        if not os.path.exists(options.root_dir):
            os.makedirs(options.root_dir)

    if not options.model_dir:
        options.model_dir = os.path.join(options.root_dir, "model")
    if not options.lexical_dir:
        options.lexical_dir = options.model_dir
    if not options.lexicon_source_target:
        options.lexicon_source_target = os.path.join(options.lexical_dir, "lex.f2n")
    if not options.lexicon_target_source:
        options.lexicon_target_source = os.path.join(options.lexical_dir, "lex.n2f")

    cicada = CICADA(options.cicada_dir)

    indexer = None
    if options.phrase:
        indexer = IndexPhrase(cicada, model_dir=options.model_dir, cky=options.cky, reordering=options.reordering)
    elif options.scfg:
        indexer = IndexSCFG(cicada, model_dir=options.model_dir, cky=options.cky, reordering=options.reordering)
    elif options.ghkm:
        indexer = IndexGHKM(cicada, model_dir=options.model_dir, cky=options.cky, reordering=options.reordering)
    elif options.tree:
        indexer = IndexTree(cicada, model_dir=options.model_dir, cky=options.cky, reordering=options.reordering)
    else:
        raise ValueError, "no indexer?"

    scores = Scores(indexer)
    features = Features(root=options.feature_root,
                        types=options.feature_type,
                        singleton=options.feature_singleton,
                        cross=options.feature_cross,
                        unaligned=options.feature_unaligned,
                        internal=options.feature_internal,
                        height=options.feature_height)
    lexicon = Lexicon(lexicon_source_target=options.lexicon_source_target,
                      lexicon_target_source=options.lexicon_target_source,
                      lexicon=options.feature_lexicon,
                      model1=options.feature_model1,
                      noisy_or=options.feature_noisy_or,
                      insertion_deletion=options.feature_insertion_deletion,
                      threshold_insertion=options.threshold_insertion,
                      threshold_deletion=options.threshold_deletion)
    
    fp = open(os.path.join(indexer.index, "files"), 'w')

    if options.pbs:
        # we use pbs to run jobs
        pbs = PBS(queue=options.pbs_queue)
    
        for score in scores:
            index = Index(cicada=cicada,
                          indexer=indexer,
                          lexicon=lexicon,
                          feats=features,
                          input=score.input,
                          output=score.output,
                          name=score.name,
                          root_joint=scores.root_joint,
                          root_source=scores.root_source,
                          root_target=scores.root_target,
                          prior=options.prior,
                          kbest=options.kbest,
                          quantize=options.quantize,
                          features=options.feature,
                          attributes=options.attribute)
            
            fp.write(os.path.basename(score.output)+'\n')

            print str(index), '2> %s'%(index.logfile)
            
            pbs.run(command=index, threads=index.threads, memory=options.max_malloc, name=index.name, logfile=index.logfile)
    
    elif options.mpi:
        mpi = MPI(cicada=cicada,
                  dir=options.mpi_dir,
                  hosts=options.mpi_host,
                  hosts_file=options.mpi_host_file,
                  number=options.mpi)
        
        for score in scores:
            index = Index(cicada=cicada,
                          indexer=indexer,
                          lexicon=lexicon,
                          feats=features,
                          input=score.input,
                          output=score.output,
                          name=score.name,
                          root_joint=scores.root_joint,
                          root_source=scores.root_source,
                          root_target=scores.root_target,
                          prior=options.prior,
                          kbest=options.kbest,
                          quantize=options.quantize,
                          features=options.feature,
                          attributes=options.attribute)

            fp.write(os.path.basename(score.output)+'\n')

            print str(index), '2> %s'%(index.logfile)

            mpi.run(command=index, logfile=index.logfile)
    
    else:
        threads = Threads(cicada=cicada, threads=options.threads)
    
        for score in scores:
            index = Index(cicada=cicada,
                          indexer=indexer,
                          lexicon=lexicon,
                          feats=features,
                          input=score.input,
                          output=score.output,
                          name=score.name,
                          root_joint=scores.root_joint,
                          root_source=scores.root_source,
                          root_target=scores.root_target,
                          prior=options.prior,
                          kbest=options.kbest,
                          quantize=options.quantize,
                          features=options.feature,
                          attributes=options.attribute)
        
            fp.write(os.path.basename(score.output)+'\n')

            print str(index), '2> %s'%(index.logfile)

            threads.run(command=index, logfile=index.logfile)
