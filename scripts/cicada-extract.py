#!/usr/bin/env python
#
#  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
#
### a wrapper script (similar to phrase-extract in moses)
### we support "lexicon construction" and "phrase/rule extraction" meaning only step 4 through 6
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
    make_option("--model-dir", default="", action="store", type="string",
                metavar="DIRECTORY", help="model directory (default: ${root_dir}/model)"),
    make_option("--alignment-dir", default="", action="store", type="string",
                metavar="DIRECTORY", help="alignment directory (default: ${model_dir})"),
    make_option("--lexical-dir", default="", action="store", type="string",
                metavar="DIRECTORY", help="lexical transltion table directory (default: ${model_dir})"),
    make_option("--counts-dir", default="", action="store", type="string",
                metavar="DIRECTORY", help="grammar counts directory (default: ${model_dir})"),
    make_option("--score-dir", default="", action="store", type="string",
                metavar="DIRECTORY", help="grammar score directory (default: ${model_dir})"),

    make_option("--temporary-dir", default="", action="store", type="string",
                metavar="DIRECTORY", help="temporary directory"),

    ### source/target flags
    make_option("--f", default="", action="store", type="string",
                metavar="FILE-OR-SUFFIX", help="source (or 'French')  language file or suffix"),
    make_option("--e", default="", action="store", type="string",
                metavar="FILE-OR-SUFFIX", help="target (or 'English') language file or suffix"),
    make_option("--a", default="", action="store", type="string",
                metavar="FILE-OR-SUFFIX", help="alignment file or suffix"),

    ### span...
    make_option("--sf", default="", action="store", type="string",
                metavar="FILE-OR-SUFFIX", help="source (or 'French')  span file or suffix"),
    make_option("--se", default="", action="store", type="string",
                metavar="FILE-OR-SUFFIX", help="target (or 'English') span file or suffix"),
    ### forest!
    make_option("--ff", default="", action="store", type="string",
                metavar="FILE-OR-SUFFIX", help="source (or 'French')  forest file or suffix"),
    make_option("--fe", default="", action="store", type="string",
                metavar="FILE-OR-SUFFIX", help="target (or 'English') forest file or suffix"),
    
    # data prefix
    make_option("--corpus", default="", action="store", type="string",
                help="bilingual trainging corpus prefix"),

    # alignment method
    make_option("--alignment", default="grow-diag-final-and", action="store", type="string",
                help="alignment methods (default: grow-diag-final-and)"),
    
    # steps
    make_option("--first-step", default=4, action="store", type="int", metavar='STEP', help="first step (default: %default)"),
    make_option("--last-step",  default=6, action="store", type="int", metavar='STEP', help="last step  (default: %default)"),

    ## option for lexicon
    make_option("--lexicon-inverse", default=None, action="store_true", help="use inverse alignment"),
    make_option("--lexicon-prior", default=0.1, action="store", type="float", metavar="PRIOR", help="lexicon model prior (default: %default)"),
    make_option("--lexicon-variational", default=None, action="store_true", help="variational Bayes estimates"),
    make_option("--lexicon-l0",          default=None, action="store_true", help="L0 regularization"),
    make_option("--lexicon-l0-alpha", default=100, action="store", type="float", help="L0 regularization parameter (default: %default)"),
    make_option("--lexicon-l0-beta",  default=0.01, action="store", type="float", help="L0 regularization parameter (default: %default)"),

    # option for extraction
    make_option("--phrase", default=None, action="store_true", help="extract phrase"),
    make_option("--scfg",   default=None, action="store_true", help="extract SCFG"),
    make_option("--ghkm",   default=None, action="store_true", help="extract GHKM (tree-to-string)"),
    make_option("--tree",   default=None, action="store_true", help="extract tree-to-tree"),

    make_option("--non-terminal", default="[x]", action="store", type="string", help="default non-terminal for GHKM rule (default: %default)"),

    make_option("--max-sentence-length", default=0, action="store", type="int",
                metavar="LENGTH", help="maximum sentence size (default: 0 == no limit)"),
    
    make_option("--max-span-source", default=15, action="store", type="int",
                metavar="LENGTH", help="maximum source span size (default: %default)"),
    make_option("--max-span-target", default=15, action="store", type="int",
                metavar="LENGTH", help="maximum target span size (default: %default)"),
    make_option("--min-hole-source", default=1, action="store", type="int",
                metavar="LENGTH", help="minimum source hole size (default: %default)"),
    make_option("--min-hole-target", default=1, action="store", type="int",
                metavar="LENGTH", help="minimum target hole size (default: %default)"),
    make_option("--max-length", default=7, action="store", type="int",
                metavar="LENGTH", help="maximum terminal length (default: %default)"),
    make_option("--max-fertility", default=4, action="store", type="int",
                metavar="FERTILITY", help="maximum terminal fertility (default: %default)"),
    make_option("--max-nodes", default=15, action="store", type="int",
                metavar="NODES", help="maximum rule nodes (default: %default)"),
    make_option("--max-height", default=4, action="store", type="int",
                metavar="HEIGHT", help="maximum rule height (default: %default)"),
    make_option("--max-compose", default=0, action="store", type="int",
                metavar="COMPOSE", help="maximum rule composition (default: %default)"),
    make_option("--max-rank", default=2, action="store", type="int",
                metavar="RANK", help="maximum rule rank (default: %default)"),
    make_option("--max-scope", default=0, action="store", type="int",
                metavar="SCOPE", help="maximum rule scope (default: %default)"),
    make_option("--cutoff", default=0.0, action="store", type="float",
                help="cutoff counts (default: %default)"),
    make_option("--collapse-source", default=None, action="store_true",
                help="collapse source side for CKY parsing"),
    make_option("--collapse-target", default=None, action="store_true",
                help="collapse target side for CKY parsing"),
    make_option("--exhaustive", default=None, action="store_true",
                help="exhaustive extraction in SCFG, GHKM and Tree"),
    make_option("--constrained", default=None, action="store_true",
                help="constrained extraction in SCFG, GHKM and Tree"),
    make_option("--project", default=None, action="store_true",
                help="project non-terminal symbols in GHKM"),
    
    make_option("--sentential", default=None, action="store_true",
                help="extract sentential rule"),
    
    ## max-malloc
    make_option("--max-malloc", default=8, action="store", type="float",
                metavar="MALLOC", help="maximum memory in GB (default: %default)"),

    # CICADA Toolkit directory
    make_option("--cicada-dir", default="", action="store", type="string",
                metavar="DIRECTORY", help="cicada directory"),
    # MPI Implementation.. if different from standard location...
    make_option("--mpi-dir", default="", action="store", type="string",
                metavar="DIRECTORY", help="MPI directory"),

    # perform threading or MPI training    
    make_option("--mpi", default=0, action="store", type="int",
                help="# of processes for MPI-based parallel processing. Identical to --np for mpirun"),
    make_option("--mpi-host", default="", action="store", type="string",
                help="list of hosts to run job. Identical to --host for mpirun", metavar="HOSTS"),
    make_option("--mpi-host-file", default="", action="store", type="string",
                help="host list file to run job. Identical to --hostfile for mpirun", metavar="FILE"),
    
    make_option("--threads", default=2, action="store", type="int",
                help="# of thrads for thread-based parallel processing"),
    
    make_option("--pbs", default=None, action="store_true",
                help="PBS for launching processes"),
    make_option("--pbs-queue", default="ltg", action="store", type="string",
                help="PBS queue for launching processes (default: %default)", metavar="NAME"),

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


class PBS:
    def __init__(self, queue=""):
        self.queue = queue
        self.pbs = 'pbs'

    def run(self, command="", threads=1, memory=0.0, name="name", mpi=None, logfile=None):
        popen = subprocess.Popen("qsub -S /bin/sh", shell=True, stdin=subprocess.PIPE)

        pipe = popen.stdin
        
        pipe.write("#!/bin/sh\n")
        pipe.write("#PBS -S /bin/sh\n")
        pipe.write("#PBS -N %s\n" %(name))
        pipe.write("#PBS -W block=true\n")
        pipe.write("#PBS -e localhost:/dev/null\n")
        pipe.write("#PBS -o localhost:/dev/null\n")
        
        if self.queue:
            pipe.write("#PBS -q %s\n" %(self.queue))

        mem = ""
        if memory >= 1.0:
            mem=":mem=%dgb" %(int(memory))
        elif memory >= 0.001:
            mem=":mem=%dmb" %(int(memory * 1000))
        elif memory >= 0.000001:
            mem=":mem=%dkb" %(int(memory * 1000 * 1000))
        
        if mpi:
            pipe.write("#PBS -l select=%d:ncpus=%d:mpiprocs=1%s\n" %(mpi.number, threads, mem))
        else:
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

        prefix = ''
        if mpi:
            prefix = mpi.mpirun
            
            if mpi.dir:
                prefix += ' --prefix %s' %(mpi.dir)
            if os.environ.has_key('TMPDIR_SPEC'):
                prefix += ' -x TMPDIR_SPEC'
            if os.environ.has_key('LD_LIBRARY_PATH'):
                prefix += ' -x LD_LIBRARY_PATH'
            if os.environ.has_key('DYLD_LIBRARY_PATH'):
                prefix += ' -x DYLD_LIBRARY_PATH'
            prefix += ' '
        
        suffix = ''
        if logfile:
            suffix = " 2> %s" %(logfile)
        
        pipe.write(prefix + command + suffix + '\n')
        
        pipe.close()
        popen.wait()
            
class MPI:
    
    def __init__(self, dir="", hosts="", hosts_file="", number=0):
        
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
                
    def run(self, command, logfile=None):
        mpirun = self.mpirun
        if self.dir:
            mpirun += ' --prefix %s' %(self.dir)
        if self.number > 0:
            mpirun += ' --np %d' %(self.number)
        if self.hosts:
            mpirun += ' --host %s' %(self.hosts)
        elif self.hosts_file:
            mpirun += ' --hostfile %s' %(self.hosts_file)

        if os.environ.has_key('TMPDIR_SPEC'):
            mpirun += ' -x TMPDIR_SPEC'
        if os.environ.has_key('LD_LIBRARY_PATH'):
            mpirun += ' -x LD_LIBRARY_PATH'
        if os.environ.has_key('DYLD_LIBRARY_PATH'):
            mpirun += ' -x DYLD_LIBRARY_PATH'

	mpirun += ' ' + command

        if logfile:
            mpirun += " 2> %s" %(logfile)
            
	run_command(mpirun)

class QSub:
    def __init__(self, mpi=None, pbs=None):
        self.mpi = mpi
        self.pbs = pbs
        
    def run(self, command, name="name", memory=0.0, threads=1, logfile=None):
        if logfile:
            print str(command), '2> %s' %(logfile)
        else:
            print str(command)

        if self.pbs:
            self.pbs.run(str(command), name=name, memory=memory, threads=threads, logfile=logfile)
        else:
            if logfile:
                run_command(str(command) + " 2> %s" %(logfile))
            else:
                run_command(str(command))
    
    def mpirun(self, command, name="name", memory=0.0, threads=1, logfile=None):
        if not self.mpi:
            raise ValueError, "no mpi?"

        if logfile:
            print str(command), '2> %s' %(logfile)
        else:
            print str(command)

        if self.pbs:
            self.pbs.run(str(command), name=name, memory=memory, mpi=self.mpi, threads=threads, logfile=logfile)
        else:
            self.mpi.run(str(command), logfile=logfile)

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
	
        for binprog in ('cicada_alignment',
                        ## step 4
                        'cicada_lexicon', 
                        ## step 5
                        'cicada_extract_phrase', 'cicada_extract_phrase_mpi',
                        'cicada_extract_scfg',   'cicada_extract_scfg_mpi',
                        'cicada_extract_ghkm',   'cicada_extract_ghkm_mpi',
                        'cicada_extract_tree',   'cicada_extract_tree_mpi',
                        ## step6
                        'cicada_extract_score', 'cicada_extract_score_mpi',):
	    
	    for bindir in bindirs:
		prog = os.path.join(bindir, binprog)
                
                if not os.path.exists(prog): continue
                if os.path.isdir(prog): continue
                
                setattr(self, binprog, prog)
                break

	    if not hasattr(self, binprog):
		raise ValueError, binprog + ' does not exist'
        
class Corpus:

    def __init__(self, corpus="", f="", e="", a="", sf="", se="", ff="", fe=""):
        
        if not corpus:
            # Directly specify data
            self.source = compressed_file(f)
            self.target = compressed_file(e)
            self.alignment = compressed_file(a)
            
            self.source_span = compressed_file(sf)
            self.target_span = compressed_file(se)
        
            self.source_forest = compressed_file(ff)
            self.target_forest = compressed_file(fe)
        else:
            # Moses style access
            self.source = compressed_file(corpus+'.'+f)
            self.target = compressed_file(corpus+'.'+e)
            self.alignment = compressed_file(corpus+'.'+a)
            
            self.source_span = compressed_file(corpus+'.'+sf)
            self.target_span = compressed_file(corpus+'.'+se)
            
            self.source_forest = compressed_file(corpus+'.'+ff)
            self.target_forest = compressed_file(corpus+'.'+fe)

class Alignment:
    def __init__(self, corpus=None, alignment_dir="", alignment=""):
        
        if os.path.exists(corpus.alignment):
            self.alignment = corpus.alignment
        else:
            self.alignment = compressed_file(os.path.join(alignment_dir, 'aligned.'+alignment))

            if not os.path.exists(self.alignment):
                raise ValueError, "no alignment data %s" %(self.alignment)

class Lexicon:
    def __init__(self, cicada=None, corpus=None, alignment=None, lexical_dir="",
                 prior=0.1,
                 variational=None,
                 l0=None,
                 l0_alpha=10,
                 l0_beta=0.5,
                 inverse=None,
                 max_malloc=4,
                 threads=4, mpi=None, pbs=None,
                 debug=None):
        self.threads = threads
        self.mpi = mpi
        self.pbs = pbs
        
        self.max_malloc = max_malloc

        self.source_target = compressed_file(os.path.join(lexical_dir, 'lex.f2n'))
        self.target_source = compressed_file(os.path.join(lexical_dir, 'lex.n2f'))
        self.makedirs = lexical_dir
        self.data = []

        self.data.append(corpus.source)
        self.data.append(corpus.target)

        command = "%s" %(cicada.cicada_lexicon)
        
        command += " --source \"%s\"" %(corpus.source)
        command += " --target \"%s\"" %(corpus.target)
        command += " --alignment \"%s\"" %(alignment.alignment)
        
        command += " --output-source-target \"%s.gz\"" %(os.path.join(lexical_dir, 'lex.f2n'))
        command += " --output-target-source \"%s.gz\"" %(os.path.join(lexical_dir, 'lex.n2f'))

        if variational:
            command += " --variational-bayes"
        if l0:
            command += " --pgd"
        
        command += " --prior %g" %(prior)
        command += " --l0-alpha %g" %(l0_alpha)
        command += " --l0-beta %g" %(l0_beta)
        
        if inverse:
            command += " --inverse"

        command += " --threads %d" %(threads)

        if debug:
            command += " --debug=%d" %(debug)
        else:
            command += " --debug"
        
        self.command = command

    def run(self):
        for file in self.data:
            if not os.path.exists(file):
                raise ValueError, "no file: " + file

        if not os.path.exists(self.makedirs):
            os.makedirs(self.makedirs)

        logfile = os.path.join(self.makedirs, "extract-lexicon.log")

        QSub(mpi=self.mpi, pbs=self.pbs).run(self.command,
                                             threads=self.threads,
                                             name="lexicon",
                                             memory=self.max_malloc,
                                             logfile=logfile)
        
        self.source_target = compressed_file(self.source_target)
        self.target_source = compressed_file(self.target_source)

class Extract:
    def __init__(self, max_malloc=8, threads=4, mpi=None, pbs=None, makedirs=""):
        self.threads = threads
        self.max_malloc = max_malloc
        self.mpi = mpi
        self.pbs = pbs
        self.makedirs = makedirs
        self.data = []
        self.logfile = ""
        self.name = ""
        
        if not hasattr(self, 'command'):
            self.command = ""
        
    def run(self):
        for file in self.data:
            if not os.path.exists(file):
                raise ValueError, "no file: " + file

        if not os.path.exists(self.makedirs):
            os.makedirs(self.makedirs)

        if not self.name:
            self.name = "extract"

        logfile = os.path.join(self.makedirs, self.logfile)

        qsub = QSub(mpi=self.mpi, pbs=self.pbs)

        if self.mpi:
            qsub.mpirun(self.command, name=self.name, memory=self.max_malloc, logfile=logfile, threads=2)
        else:
            qsub.run(self.command, name=self.name, memory=self.max_malloc, logfile=logfile, threads=self.threads)

class ExtractPhrase(Extract):
    
    def __init__(self, cicada=None, corpus=None, alignment=None,
                 counts_dir="",
                 max_length=7, max_fertility=4,
                 max_malloc=8, threads=4, mpi=None, pbs=None,
                 debug=None):
        Extract.__init__(self, max_malloc, threads, mpi, pbs, counts_dir)
        
        self.counts = os.path.join(counts_dir, "phrase-counts")

        self.data.append(corpus.source)
        self.data.append(corpus.target)
        self.logfile = "extract-phrase.log"
        self.name = "extract-phrase"
        
        prog_name = cicada.cicada_extract_phrase
        if mpi:
            prog_name = cicada.cicada_extract_phrase_mpi
        
        command = prog_name
        
        command += " --source \"%s\"" %(corpus.source)
        command += " --target \"%s\"" %(corpus.target)
        command += " --alignment \"%s\"" %(alignment.alignment)
        
        command += " --output \"%s\"" %(self.counts)
        
        command += " --max-length %d"    %(max_length)
        command += " --max-fertility %d" %(max_fertility)
        
        command += " --max-malloc %g" %(max_malloc)
        
        if not mpi:
            command += " --threads %d" %(threads)

        if debug:
            command += " --debug=%d" %(debug)
        else:
            command += " --debug"
        
        self.command = command

class ExtractSCFG(Extract):
    
    def __init__(self, cicada=None, corpus=None, alignment=None,
                 counts_dir="",
                 max_length=7, max_fertility=4,
                 max_span_source=15, max_span_target=20,
                 min_hole_source=1, min_hole_target=1,
                 max_rank=2,
                 exhaustive=None, constrained=None, sentential=None,
                 max_malloc=8, threads=4, mpi=None, pbs=None,
                 debug=None):
        Extract.__init__(self, max_malloc, threads, mpi, pbs, counts_dir)
        
        self.counts = os.path.join(counts_dir, "scfg-counts")
        
        self.data.append(corpus.source)
        self.data.append(corpus.target)
        self.logfile = "extract-scfg.log"
        self.name = "extract-scfg"
        
        if os.path.exists(corpus.source_span) and os.path.exists(corpus.target_span):
            raise ValueError, "both of source/target span specified... which one?"
        
        prog_name = cicada.cicada_extract_scfg
        if mpi:
            prog_name = cicada.cicada_extract_scfg_mpi
        
        command = prog_name
        
        command += " --source \"%s\"" %(corpus.source)
        command += " --target \"%s\"" %(corpus.target)
        command += " --alignment \"%s\"" %(alignment.alignment)

        if os.path.exists(corpus.source_span):
            command += " --spans-source \"%s\"" %(corpus.source_span)
        if os.path.exists(corpus.target_span):
            command += " --spans-target \"%s\"" %(corpus.target_span)
        
        command += " --output \"%s\"" %(self.counts)
        
        command += " --max-length %d"    %(max_length)
        command += " --max-fertility %d" %(max_fertility)
        command += " --max-span-source %d"      %(max_span_source)
        command += " --max-span-target %d"      %(max_span_target)
        command += " --min-hole-source %d"      %(min_hole_source)
        command += " --min-hole-target %d"      %(min_hole_target)
        command += " --max-rank %d" %(max_rank)
        
        if exhaustive:
            command += " --exhaustive"
        if constrained:
            command += " --constrained"
        if sentential:
            command += " --sentential"
        
        command += " --max-malloc %g" %(max_malloc)

        if not mpi:
            command += " --threads %d" %(self.threads)        
        if debug:
            command += " --debug=%d" %(debug)
        else:
            command += " --debug"
        
        self.command = command

class ExtractGHKM(Extract):
    
    def __init__(self, cicada=None, corpus=None, alignment=None,
                 counts_dir="",
                 non_terminal="", max_sentence_length=0, max_nodes=15, max_height=4, max_compose=0, max_scope=0,
                 cutoff=0.0,
                 collapse_source=None,
                 collapse_target=None,
                 exhaustive=None,
                 constrained=None,
                 project=None,
                 max_malloc=8, threads=4, mpi=None, pbs=None,
                 debug=None):
        Extract.__init__(self, max_malloc, threads, mpi, pbs, counts_dir)
        
        self.counts = os.path.join(counts_dir, "ghkm-counts")

        if os.path.exists(corpus.source_forest) and os.path.exists(corpus.target_forest):
            raise ValueError, "both source and target forest.. we can extract string-to-tree or tree-to-string"
        
        tree_to_string = 1
        if os.path.exists(corpus.target_forest):
            tree_to_string = None
            
        if tree_to_string:
            self.data.append(corpus.source_forest)
            self.data.append(corpus.target)
        else:
            self.data.append(corpus.target_forest)
            self.data.append(corpus.source)
            
        self.logfile = "extract-ghkm.log"
        self.name = "extract-ghkm"
        
        prog_name = cicada.cicada_extract_ghkm
        if mpi:
            prog_name = cicada.cicada_extract_ghkm_mpi
        
        command = prog_name
        
        if tree_to_string:
            command += " --source \"%s\"" %(corpus.source_forest)
            command += " --target \"%s\"" %(corpus.target)
        else:
            ## strig-to-tree extraction...!
            command += " --source \"%s\"" %(corpus.target_forest)
            command += " --target \"%s\"" %(corpus.source)
            command += " --inverse"
            command += " --swap"
            
        command += " --alignment \"%s\"" %(alignment.alignment)
        
        command += " --output \"%s\"" %(self.counts)
        
        if non_terminal:
            if non_terminal[0] != '[' or non_terminal[-1] != ']':
                raise ValueError, "invalid non-terminal: %s" %(non_terminal)

            command += " --non-terminal \"%s\"" %(non_terminal)

        if max_sentence_length > 0:
            command += " --max-sentence-length %d"   %(max_sentence_length)
        
        command += " --max-nodes %d"   %(max_nodes)
        command += " --max-height %d"  %(max_height)
        command += " --max-compose %d" %(max_compose)
        command += " --max-scope %d"   %(max_scope)

        command += " --cutoff %.20g" %(cutoff)
        
        if collapse_source:
            command += " --collapse-source"
        if collapse_target:
            command += " --collapse-target"
        if exhaustive:
            command += " --exhaustive"
        if constrained:
            command += " --constrained"
        if project:
            command += " --project"
        
        command += " --max-malloc %g" %(max_malloc)

        if not mpi:
            command += " --threads %d" %(self.threads)
        if debug:
            command += " --debug=%d" %(debug)
        else:
            command += " --debug"

        self.command = command


class ExtractTree(Extract):
    
    def __init__(self, cicada=None, corpus=None, alignment=None,
                 counts_dir="",
                 max_sentence_length=0, max_nodes=15, max_height=4, max_compose=0, max_scope=0,
                 cutoff=0.0,
                 collapse_source=None,
                 collapse_target=None,
                 exhaustive=None,
                 constrained=None,
                 max_malloc=8, threads=4, mpi=None, pbs=None,
                 debug=None):
        Extract.__init__(self, max_malloc, threads, mpi, pbs, counts_dir)
        
        self.counts = os.path.join(counts_dir, "tree-counts")

        self.data.append(corpus.source_forest)
        self.data.append(corpus.target_forest)
        self.logfile = "extract-tree.log"
        self.name = "extract-tree"
        
        prog_name = cicada.cicada_extract_tree
        if mpi:
            prog_name = cicada.cicada_extract_tree_mpi
        
        command = prog_name
        
        command += " --source \"%s\"" %(corpus.source_forest)
        command += " --target \"%s\"" %(corpus.target_forest)
        command += " --alignment \"%s\"" %(alignment.alignment)
        
        command += " --output \"%s\"" %(self.counts)
        
        command += " --max-nodes %d"   %(max_nodes)
        command += " --max-height %d"  %(max_height)
        command += " --max-compose %d" %(max_compose)
        command += " --max-scope %d"   %(max_scope)
        
        command += " --cutoff %.20g" %(cutoff)
        
        if collapse_source:
            command += " --collapse-source"
        if collapse_target:
            command += " --collapse-target"
        if exhaustive:
            command += " --exhaustive"
        if constrained:
            command += " --constrained"
        
        command += " --max-malloc %g" %(max_malloc)
        
        if not mpi:
            command += " --threads %d" %(self.threads)
            
        if debug:
            command += " --debug=%d" %(debug)
        else:
            command += " --debug"
        
        self.command = command

class ExtractScore(Extract):
    
    def __init__(self, cicada=None,
                 counts_dir="",
                 score_dir="",
                 temporary_dir="",
                 phrase=None, scfg=None, ghkm=None, tree=None,
                 max_malloc=8, threads=4, mpi=None, pbs=None,
                 debug=None):
        Extract.__init__(self, max_malloc, threads, mpi, pbs, score_dir)
        
        option = ""
        if phrase:
            self.counts = os.path.join(counts_dir, "phrase-counts")
            self.scores = os.path.join(score_dir,  "phrase-score")
            self.logfile = "extract-score.phrase.log"
            option = " --score-phrase"
        elif scfg:
            self.counts = os.path.join(counts_dir, "scfg-counts")
            self.scores = os.path.join(score_dir,  "scfg-score")
            self.logfile = "extract-score.scfg.log"
            option = " --score-scfg"
        elif ghkm:
            self.counts = os.path.join(counts_dir, "ghkm-counts")
            self.scores = os.path.join(score_dir,  "ghkm-score")
            self.logfile = "extract-score.ghkm.log"
            option = " --score-ghkm"
        elif tree:
            self.counts = os.path.join(counts_dir, "tree-counts")
            self.scores = os.path.join(score_dir,  "tree-score")
            self.logfile = "extract-score.tree.log"
            option = " --score-ghkm"
        else:
            raise ValueError, "no count type?"

        if not os.path.exists(self.counts):
            raise ValueError, "no counts? %s" %(self.counts)

        self.name = "extract-score"
        
        prog_name = cicada.cicada_extract_score
        if mpi:
            prog_name = cicada.cicada_extract_score_mpi
        
        command = prog_name
        
        command += " --input \"%s\"" %(self.counts)
        command += " --output \"%s\"" %(self.scores)
        command += option
        command += " --max-malloc %g" %(max_malloc)
        
        if temporary_dir:
            command += " --temporary \"%s\"" %(temporary_dir)
        
        if mpi:
            command += " --prog \"%s\"" %(prog_name)
        else:
            command += " --threads %d" %(self.threads)
            
        if debug:
            command += " --debug=%d" %(debug)
        else:
            command += " --debug"

        self.command = command

if __name__ == '__main__':
    (options, args) = opt_parser.parse_args()

    ### dump to stderr
    stdout = sys.stdout
    sys.stdout = sys.stderr

    if options.root_dir:
        if not os.path.exists(options.root_dir):
            os.makedirs(options.root_dir)

    if not options.model_dir:
        options.model_dir = os.path.join(options.root_dir, "model")
    if not options.lexical_dir:
        options.lexical_dir = options.model_dir
    if not options.alignment_dir:
        options.alignment_dir = options.model_dir
    if not options.counts_dir:
        options.counts_dir = options.model_dir
    if not options.score_dir:
        options.score_dir = options.model_dir

    if not options.temporary_dir:
        if os.environ.has_key('TMPDIR_SPEC') and os.environ['TMPDIR_SPEC']:
            options.temporary_dir = os.environ['TMPDIR_SPEC']
    else:
        os.environ['TMPDIR_SPEC'] = options.temporary_dir
        
    cicada = CICADA(options.cicada_dir)

    mpi = None
    if options.mpi_host or options.mpi_host_file or options.mpi > 0:
        mpi = MPI(dir=options.mpi_dir,
                  hosts=options.mpi_host,
                  hosts_file=options.mpi_host_file,
                  number=options.mpi)
        
    pbs = None
    if options.pbs:
        pbs = PBS(queue=options.pbs_queue)

    corpus = Corpus(corpus=options.corpus,
                    f=options.f,
                    e=options.e,
                    a=options.a,
                    sf=options.sf,
                    se=options.se,
                    ff=options.ff,
                    fe=options.fe)

    alignment = Alignment(corpus=corpus, alignment_dir=options.alignment_dir, alignment=options.alignment)

    lexicon = Lexicon(cicada=cicada, corpus=corpus, alignment=alignment,
                      lexical_dir=options.lexical_dir,
                      prior=options.lexicon_prior,
                      variational=options.lexicon_variational,
                      l0=options.lexicon_l0,
                      l0_alpha=options.lexicon_l0_alpha,
                      l0_beta=options.lexicon_l0_beta,
                      inverse=options.lexicon_inverse,
                      threads=options.threads, mpi=mpi, pbs=pbs,
                      debug=options.debug)

    if options.first_step <= 4 and options.last_step >= 4:
        print "(4) generate lexical translation table started  @", time.ctime()
        lexicon.run()
        print "(4) generate lexical translation table finished @", time.ctime()

    if options.first_step <= 5 and options.last_step >= 5:
        extract = None
        if options.phrase:
            extract = ExtractPhrase(cicada=cicada, corpus=corpus, alignment=alignment,
                                    counts_dir=options.counts_dir,
                                    max_length=options.max_length,
                                    max_fertility=options.max_fertility,
                                    max_malloc=options.max_malloc, threads=options.threads, mpi=mpi, pbs=pbs,
                                    debug=options.debug)
        elif options.scfg:
            extract = ExtractSCFG(cicada=cicada, corpus=corpus, alignment=alignment,
                                  counts_dir=options.counts_dir,
                                  max_length=options.max_length,
                                  max_fertility=options.max_fertility,
                                  max_span_source=options.max_span_source,
                                  max_span_target=options.max_span_target,
                                  min_hole_source=options.min_hole_source,
                                  min_hole_target=options.min_hole_target,
                                  max_rank=options.max_rank,
                                  exhaustive=options.exhaustive,
                                  constrained=options.constrained,
                                  sentential=options.sentential,
                                  max_malloc=options.max_malloc, threads=options.threads, mpi=mpi, pbs=pbs,
                                  debug=options.debug)
        elif options.ghkm:
            extract = ExtractGHKM(cicada=cicada, corpus=corpus, alignment=alignment,
                                  counts_dir=options.counts_dir,
                                  non_terminal=options.non_terminal,
                                  max_sentence_length=options.max_sentence_length,
                                  max_nodes=options.max_nodes,
                                  max_height=options.max_height,
                                  max_compose=options.max_compose,
                                  max_scope=options.max_scope,
                                  cutoff=options.cutoff,
                                  collapse_source=options.collapse_source,
                                  collapse_target=options.collapse_target,
                                  exhaustive=options.exhaustive,
                                  constrained=options.constrained,
                                  project=options.project,
                                  max_malloc=options.max_malloc, threads=options.threads, mpi=mpi, pbs=pbs,
                                  debug=options.debug)
        elif options.tree:
            extract = ExtractTree(cicada=cicada, corpus=corpus, alignment=alignment,
                                  counts_dir=options.counts_dir,
                                  max_sentence_length=options.max_sentence_length,
                                  max_nodes=options.max_nodes,
                                  max_height=options.max_height,
                                  max_compose=options.max_compose,
                                  max_scope=options.max_scope,
                                  cutoff=options.cutoff,
                                  collapse_source=options.collapse_source,
                                  collapse_target=options.collapse_target,
                                  exhaustive=options.exhaustive,
                                  constrained=options.constrained,
                                  max_malloc=options.max_malloc, threads=options.threads, mpi=mpi, pbs=pbs,
                                  debug=options.debug)
        else:
            raise ValueError, "no count type?"

        print "(5) extract phrase table started @", time.ctime()
        extract.run()
        print "(5) extract phrase table finished @", time.ctime()

    if options.first_step <= 6 and options.last_step >= 6:
        score = ExtractScore(cicada=cicada,
                             counts_dir=options.counts_dir,
                             score_dir=options.score_dir,
                             temporary_dir=options.temporary_dir,
                             phrase=options.phrase, scfg=options.scfg, ghkm=options.ghkm, tree=options.tree,
                             max_malloc=options.max_malloc, threads=options.threads, mpi=mpi, pbs=pbs,
                             debug=options.debug)
    
        print "(6) score phrase table started @", time.ctime()
        score.run()
        print "(6) score phrase table finished @", time.ctime()
