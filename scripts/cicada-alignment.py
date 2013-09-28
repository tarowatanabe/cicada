#!/usr/bin/env python
#
#  Copyright(C) 2010-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
#
### a wrapper script (similar to phrase-extract in moses)
### we support only "extraction" meaning only step 5 and 6
### TODO: use argparse for command-lines...?

import threading
import multiprocessing

import time
import sys
import os, os.path
import stat
import string
import re
import subprocess

### for find_executable!
import distutils.spawn

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
    make_option("--alignment", default="posterior-itg", action="store", type="string",
                help="alignment methods (default: %default)"),
    
    # steps
    make_option("--first-step", default=1, action="store", type="int", metavar='STEP', help="first step (default: %default)"),
    make_option("--last-step",  default=3, action="store", type="int", metavar='STEP', help="last step  (default: %default)"),
    
    ## iteratin
    make_option("--iteration-cluster", default=50, action="store", type="int", metavar='ITERATION', help="word cluter iterations (default: %default)"),
    make_option("--iteration-model1",  default=5,  action="store", type="int", metavar='ITERATION', help="Model1 iteratins (default: %default)"),
    make_option("--iteration-hmm",     default=5,  action="store", type="int", metavar='ITERATION', help="HMM iteratins    (default: %default)"),
    make_option("--iteration-model4",  default=5,  action="store", type="int", metavar='ITERATION', help="Model4 iteratins    (default: %default)"),
    
    ## # of clusters
    make_option("--cluster",     default=50, action="store", type="int", metavar='CLUSTER', help="# of clusters (default: %default)"),
    
    ## training parameters
    make_option("--symmetric",   default=None, action="store_true", help="symmetric training"),
    make_option("--posterior",   default=None, action="store_true", help="posterior constrained training"),
    make_option("--dynamic",     default=None, action="store_true", help="dynamically recompute base alignment"),
    make_option("--variational", default=None, action="store_true", help="variational Bayes estimates"),
    make_option("--l0",          default=None, action="store_true", help="L0 regularization"),
    
    ## options for lexicon model training
    make_option("--p0",              default=0.01, action="store", type="float", metavar='P0',    help="parameter for NULL alignment (default: %default)"),
    make_option("--insertion-p1",    default=0.01, action="store", type="float", metavar='P1',    help="parameter for NULL insertion (default: %default)"),
    make_option("--prior-lexicon",   default=0.01, action="store", type="float", metavar="PRIOR", help="lexicon model prior (default: %default)"),
    make_option("--prior-alignment", default=0.01, action="store", type="float", metavar="PRIOR", help="alignment model prior (default: %default)"),
    make_option("--prior-distortion", default=0.01, action="store", type="float", metavar="PRIOR", help="distortion model prior (default: %default)"),
    make_option("--prior-fertility", default=0.01, action="store", type="float", metavar="PRIOR", help="fertility model prior (default: %default)"),
    
    make_option("--smooth-lexicon",   default=1e-100, action="store", type="float", metavar="SMOOTH", help="lower-bound parameter for lexicon model (default: %default)"),
    make_option("--smooth-alignment", default=1e-100, action="store", type="float", metavar="SMOOTH", help="lower-bound parameter for alignment model (default: %default)"),
    make_option("--smooth-distortion", default=1e-100, action="store", type="float", metavar="SMOOTH", help="lower-bound parameter for distortion model (default: %default)"),
    make_option("--smooth-fertility", default=1e-100, action="store", type="float", metavar="SMOOTH", help="lower-bound parameter for fertility model (default: %default)"),

    make_option("--l0-alpha", default=100, action="store", type="float", help="L0 regularization parameter (default: %default)"),
    make_option("--l0-beta",  default=0.01, action="store", type="float", help="L0 regularization parameter (default: %default)"),
    
    # CICADA Toolkit directory
    make_option("--cicada-dir", default="", action="store", type="string",
                metavar="DIRECTORY", help="cicada directory"),
    
    make_option("--threads", default=2, action="store", type="int",
                help="# of thrads for thread-based parallel processing"),

    ## max-malloc
    make_option("--max-malloc", default=8, action="store", type="float",
                metavar="MALLOC", help="maximum memory in GB (default: %default)"),
    ## PBS
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

        # how to find binary location...?
        if not distutils.spawn.find_executable('qsub'):
            raise ValueError, "no qsub in your executable path?"

    def run(self, command="", threads=1, memory=0.0, name="name", mpi=None, logfile=None):

        popen = subprocess.Popen(['qsub', '-S', '/bin/sh'], stdin=subprocess.PIPE)
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
	
        for binprog in ('cicada_cluster_word',
                        ## step 1
                        'cicada_alignment_model4',
                        'cicada_alignment_model1',
                        'cicada_alignment_hmm',
                        ## step2
                        'cicada_alignment', 
                        ## step 3
                        ):
	    
	    for bindir in bindirs:
                prog = os.path.join(bindir, binprog)
                
                if not os.path.exists(prog): continue
                if os.path.isdir(prog): continue

                setattr(self, binprog, prog)
                break
            
	    if not hasattr(self, binprog):
		raise ValueError, binprog + ' does not exist'
        
class Corpus:

    def __init__(self, corpus_dir="", corpus="", f="", e="", a="", sf="", se="", ff="", fe=""):
        print corpus

        if not corpus:
            self.source_tag = 'src'
            self.target_tag = 'trg'
            self.alignment_tag = 'algn'
            
            self.source_span_tag = 'span-src'
            self.target_span_tag = 'span-trg'
            
            self.source_forest_tag = 'forest-src'
            self.target_forest_tag = 'forest-trg'
            
            self.source = compressed_file(f)
            self.target = compressed_file(e)
            self.alignment = compressed_file(a)
            
            self.source_span = compressed_file(sf)
            self.target_span = compressed_file(se)
            
            self.source_forest = compressed_file(ff)
            self.target_forest = compressed_file(fe)
        else:
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

    def __init__(self, cicada=None, corpus="", name="", cluster=64, iteration=64, threads=8, mpi=None, pbs=None, debug=0):
        
        self.cicada  = cicada
        
        self.mpi = mpi
        self.pbs = pbs
        self.threads = threads
        
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
        QSub(mpi=self.mpi, pbs=self.pbs).run(self.command,
                                             threads=self.threads,
                                             name="cluster",
                                             logfile=self.cluster+'.log')

class Prepare:
    
    def __init__(self, cicada=None, corpus=None, cluster=64, iteration=64, threads=8, mpi=None, pbs=None, debug=0):
        
        if not os.path.exists(corpus.corpus_dir):
            os.makedirs(corpus.corpus_dir)

        self.corpus = corpus
        self.source = Cluster(cicada=cicada,
                              corpus=corpus.source,
                              name=os.path.join(corpus.corpus_dir, corpus.source_tag+'.vcb.classes'),
                              cluster=cluster,
                              iteration=iteration,
                              threads=threads,
                              mpi=mpi,
                              pbs=pbs,
                              debug=debug)
        
        self.target = Cluster(cicada=cicada,
                              corpus=corpus.target,
                              name=os.path.join(corpus.corpus_dir, corpus.target_tag+'.vcb.classes'),
                              cluster=cluster,
                              iteration=iteration,
                              threads=threads,
                              mpi=mpi,
                              pbs=pbs,
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
                 iteration_model4=5,
                 prior_lexicon=0.1,
                 prior_alignment=0.1,
                 prior_distortion=0.1,
                 prior_fertility=0.1,
                 smooth_lexicon=1e-20,
                 smooth_alignment=1e-20,
                 smooth_distortion=1e-20,
                 smooth_fertility=1e-20,
                 p0=1e-4,
                 insertion_p1=1e-4,
                 symmetric=None,
                 posterior=None,
                 dynamic=None,
                 variational=None,
                 l0=None,
                 l0_alpha=10,
                 l0_beta=0.5,
                 threads=8,
                 mpi=None,
                 pbs=None,
                 debug=0):

        self.mpi = mpi
        self.pbs = pbs
        self.threads = threads
        
        if not os.path.exists(dir_source_target):
            os.makedirs(dir_source_target)
        if not os.path.exists(dir_target_source):
            os.makedirs(dir_target_source)
        
        command = ""
        if iteration_model4 > 0:
            command = cicada.cicada_alignment_model4
        elif iteration_hmm > 0:
            command = cicada.cicada_alignment_hmm
        elif iteration_model1 > 0:
            command = cicada.cicada_alignment_model1
        else:
            raise ValueError, "invalid model iterations"
        
        command += " --source \"%s\"" %(corpus.source)
        command += " --target \"%s\"" %(corpus.target)

        if os.path.exists(corpus.alignment):
            command += " --alignment \"%s\"" %(corpus.alignment)

        if iteration_hmm > 0 or iteration_model4 > 0:
            command += " --classes-source \"%s\"" %(compressed_file(cluster.source.cluster))
            command += " --classes-target \"%s\"" %(compressed_file(cluster.target.cluster))
        
        if iteration_model4 > 0:
            self.distortion_source_target = os.path.join(dir_source_target, prefix_source_target + '.distortion.final.gz')
            self.distortion_target_source = os.path.join(dir_target_source, prefix_target_source + '.distortion.final.gz')

            command += " --output-distortion-source-target \"%s\"" %(self.distortion_source_target)
            command += " --output-distortion-target-source \"%s\"" %(self.distortion_target_source)            

        if iteration_model4 > 0:
            self.fertility_source_target = os.path.join(dir_source_target, prefix_source_target + '.fertility.final.gz')
            self.fertility_target_source = os.path.join(dir_target_source, prefix_target_source + '.fertility.final.gz')

            command += " --output-fertility-source-target \"%s\"" %(self.fertility_source_target)
            command += " --output-fertility-target-source \"%s\"" %(self.fertility_target_source)
        
        if iteration_hmm > 0 or iteration_model4 > 0:
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

        self.posterior_source_target = os.path.join(dir_source_target, prefix_source_target + '.posterior.final.gz')
        self.posterior_target_source = os.path.join(dir_target_source, prefix_target_source + '.posterior.final.gz')
        
        command += " --posterior-source-target \"%s\"" %(self.viterbi_source_target)
        command += " --posterior-target-source \"%s\"" %(self.viterbi_target_source)
        

        if iteration_model4 > 0:
            command += " --iteration-model4 %d" %(iteration_model4)
            command += " --iteration-hmm %d" %(iteration_hmm)
            command += " --iteration-model1 %d" %(iteration_model1)
        elif iteration_hmm > 0:
            command += " --iteration-hmm %d" %(iteration_hmm)
            command += " --iteration-model1 %d" %(iteration_model1)
        else:
            command += " --iteration %d" %(iteration_model1)

        if symmetric:
            command += " --symmetric"
        if posterior:
            command += " --posterior"
        if dynamic and iteration_model4 > 0:
            command += " --dynamic"

        if variational:
            command += " --variational-bayes"
        if l0:
            command += " --pgd"
        
        command += " --l0-alpha %g" %(l0_alpha)
        command += " --l0-beta %g" %(l0_beta)
        
        self.p0 = p0
        self.insertion_p1 = insertion_p1
        self.prior_lexicon    = prior_lexicon
        self.prior_alignment  = prior_alignment
        self.prior_distortion  = prior_distortion
        self.prior_fertility  = prior_fertility
        self.smooth_lexicon   = smooth_lexicon
        self.smooth_alignment = smooth_alignment
        self.smooth_distortion = smooth_distortion
        self.smooth_fertility = smooth_fertility
        
        command += " --p0 %.20g" %(p0)
        command += " --prior-lexicon %.20g"   %(prior_lexicon)
        command += " --smooth-lexicon %.20g"   %(smooth_lexicon)
        if iteration_hmm > 0 or iteration_model4 > 0:
            command += " --prior-alignment %.20g" %(prior_alignment)
            command += " --smooth-alignment %.20g" %(smooth_alignment)

        if iteration_model4 > 0:
            command += " --insertion-source-target %.20g" %(insertion_p1)
            command += " --insertion-target-source %.20g" %(insertion_p1)

            command += " --prior-distortion %.20g" %(prior_distortion)
            command += " --smooth-distortion %.20g" %(smooth_distortion)
            command += " --prior-fertility %.20g" %(prior_fertility)
            command += " --smooth-fertility %.20g" %(smooth_fertility)
            
                    
        command += " --threads %d" %(threads)

        if debug:
            command += " --debug=%d" %(debug)
        else:
            command += " --debug"
        
        self.command = command

    def run(self, logfile):
        QSub(mpi=self.mpi, pbs=self.pbs).run(self.command,
                                             threads=self.threads,
                                             name="giza",
                                             logfile=logfile)

class Alignment:

    def __init__(self,
                 cicada=None,
                 corpus=None,
                 cluster=None,
                 giza=None,
                 alignment_dir="",
                 alignment="grow-diag-final-and",
                 threads=8,
                 mpi=None,
                 pbs=None,
                 debug=0):

        self.mpi = mpi
        self.pbs = pbs
        self.threads = threads

        if not os.path.exists(alignment_dir):
            os.makedirs(alignment_dir)

        posterior_mode = None
        posterior_threshold = 0.0
        if 'posterior' in alignment:
            posterior_mode = 1
            result = re.compile(r"posterior-(.+)").search(alignment)
            if result:
                try:
                    posterior_threshold = float(result.group(1))
                except:
                    posterior_threshold = 0.0
        
        command = cicada.cicada_alignment

        if posterior_mode:
            command += " --source-target \"%s\"" %(compressed_file(giza.posterior_source_target))
            command += " --target-source \"%s\"" %(compressed_file(giza.posterior_target_source))
            command += " --posterior"
            
            if posterior_threshold > 0.0:
                command += " --posterior-threshold %.20g" %(posterior_threshold)

        else:
            command += " --source-target \"%s\"" %(compressed_file(giza.viterbi_source_target))
            command += " --target-source \"%s\"" %(compressed_file(giza.viterbi_target_source))

        if os.path.exists(corpus.source_span):
            command += " --span-source \"%s\"" %(corpus.source_span)
        if os.path.exists(corpus.target_span):
            command += " --span-target \"%s\"" %(corpus.target_span)
        
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
        if 'closure' in alignment:
            command += " --closure"
        
        command += " --threads %d" %(threads)

        if debug:
            command += " --debug=%d" %(debug)
        else:
            command += " --debug"
        
        self.command = command

    def run(self):
        QSub(mpi=self.mpi, pbs=self.pbs).run(self.command,
                                             threads=self.threads,
                                             name="align-heu",
                                             logfile=self.alignment+'.log')



class Aligner:

    def __init__(self,
                 cicada=None,
                 cluster=None,
                 giza=None,
                 alignment_dir="",
                 threads=8,
                 debug=0):

        if not os.path.exists(alignment_dir):
            os.makedirs(alignment_dir)
            
        learn_hmm = None
        learn_model4 = None
        if hasattr(giza, 'alignment_source_target'):
            learn_hmm = 1
        if hasattr(giza, 'distortion_source_target'):
            learn_model4 = 1

        command = cicada.cicada_alignment_model1
        if learn_hmm:
            command = cicada.cicada_alignment_hmm
        if learn_model4:
            command = cicada.cicada_alignment_model4
        command += " \\\n"

        if learn_hmm or learn_model4:
            command += " --classes-source \"%s\"" %(os.path.realpath(compressed_file(cluster.source.cluster)))
            command += " \\\n"
            command += " --classes-target \"%s\"" %(os.path.realpath(compressed_file(cluster.target.cluster)))
            command += " \\\n"

        command += " --lexicon-source-target \"%s\"" %(os.path.realpath(compressed_file(giza.lexicon_source_target)))
        command += " \\\n"
        command += " --lexicon-target-source \"%s\"" %(os.path.realpath(compressed_file(giza.lexicon_target_source)))
        command += " \\\n"

        if learn_hmm or learn_model4:
            command += " --alignment-source-target \"%s\"" %(os.path.realpath(compressed_file(giza.alignment_source_target)))
            command += " \\\n"
            command += " --alignment-target-source \"%s\"" %(os.path.realpath(compressed_file(giza.alignment_target_source)))
            command += " \\\n"

        if learn_model4:
            command += " --insertion-source-target %.20g" %(giza.insertion_p1)
            command += " \\\n"
            command += " --insertion-target-source %.20g" %(giza.insertion_p1)
            command += " \\\n"
            command += " --distortion-source-target \"%s\"" %(os.path.realpath(compressed_file(giza.distortion_source_target)))
            command += " \\\n"
            command += " --distortion-target-source \"%s\"" %(os.path.realpath(compressed_file(giza.distortion_target_source)))
            command += " \\\n"
            command += " --fertility-source-target \"%s\"" %(os.path.realpath(compressed_file(giza.fertility_source_target)))
            command += " \\\n"
            command += " --fertility-target-source \"%s\"" %(os.path.realpath(compressed_file(giza.fertility_target_source)))
            command += " \\\n"

        if learn_model4:
            command += " --iteration-model4 0"
            command += " \\\n"
            command += " --iteration-model1 0"
            command += " \\\n"
            command += " --iteration-hmm 0"
            command += " \\\n"
        elif learn_hmm:
            command += " --iteration-model1 0"
            command += " \\\n"
            command += " --iteration-hmm 0"
            command += " \\\n"
        else:
            command += " --iteration 0"
            command += " \\\n"
            
        command += " --p0 %.20g" %(giza.p0)
        command += " \\\n"
        command += " --prior-lexicon %.20g"  %(giza.prior_lexicon)
        command += " \\\n"
        command += " --smooth-lexicon %.20g" %(giza.smooth_lexicon)
        command += " \\\n"
        if learn_hmm or learn_model4:
            command += " --prior-alignment %.20g"  %(giza.prior_alignment)
            command += " \\\n"
            command += " --smooth-alignment %.20g" %(giza.smooth_alignment)
            command += " \\\n"

        if learn_model4:
            command += " --prior-distortion %.20g"  %(giza.prior_distortion)
            command += " \\\n"
            command += " --smooth-distortion %.20g" %(giza.smooth_distortion)
            command += " \\\n"
            command += " --prior-fertility %.20g"  %(giza.prior_fertility)
            command += " \\\n"
            command += " --smooth-fertility %.20g" %(giza.smooth_fertility)
            command += " \\\n"
        
        if debug:
            command += " --debug=%d" %(debug)
            command += " \\\n"
        else:
            command += " --debug"
            command += " \\\n"
        
        self.command = command

    def run(self, fp):
        fp.write("#!/bin/sh\n")
        fp.write("\n")
        fp.write("exec ")
        fp.write(self.command)
        fp.write(" \"$@\"")

if __name__ == '__main__':
    (options, args) = opt_parser.parse_args()

    ### dump to stderr
    stdout = sys.stdout
    sys.stdout = sys.stderr

    if options.root_dir:
        if not os.path.exists(options.root_dir):
            os.makedirs(options.root_dir)

    if not options.corpus_dir:
        options.corpus_dir = os.path.join(options.root_dir, "corpus")
    if not options.model_dir:
        options.model_dir = os.path.join(options.root_dir, "model")
    if not options.alignment_dir:
        options.alignment_dir = options.model_dir

    cicada = CICADA(options.cicada_dir)

    pbs = None
    if options.pbs:
        pbs = PBS(queue=options.pbs_queue)

    corpus = Corpus(corpus_dir=options.corpus_dir,
                    corpus=options.corpus,
                    f=options.f,
                    e=options.e,
                    a=options.a,
                    sf=options.sf,
                    se=options.se,
                    ff=options.ff,
                    fe=options.fe)

    if not options.giza_f2e:
        options.giza_f2e = os.path.join(options.root_dir, "giza.%s-%s" %(corpus.source_tag, corpus.target_tag))
    if not options.giza_e2f:
        options.giza_e2f = os.path.join(options.root_dir, "giza.%s-%s" %(corpus.target_tag, corpus.source_tag))

    prepare = Prepare(cicada=cicada,
                      corpus=corpus,
                      cluster=options.cluster,
                      iteration=options.iteration_cluster,
                      threads=options.threads,
                      pbs=pbs,
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
                iteration_model4=options.iteration_model4,
                prior_lexicon=options.prior_lexicon,
                prior_alignment=options.prior_alignment,
                prior_distortion=options.prior_distortion,
                prior_fertility=options.prior_fertility,
                smooth_lexicon=options.smooth_lexicon,
                smooth_alignment=options.smooth_alignment,
                smooth_distortion=options.smooth_distortion,
                smooth_fertility=options.smooth_fertility,
                p0=options.p0,
                insertion_p1=options.insertion_p1,
                symmetric=options.symmetric,
                dynamic=options.dynamic,
                posterior=options.posterior,
                variational=options.variational,
                l0=options.l0,
                l0_alpha=options.l0_alpha,
                l0_beta=options.l0_beta,
                threads=options.threads,
                pbs=pbs,
                debug=options.debug)

    ## run giza++ in two directions
    if options.first_step <= 2 and options.last_step >= 2:
        # dump aligner...
        aligner = Aligner(cicada=cicada,
                          cluster=prepare,
                          giza=giza,
                          alignment_dir=options.alignment_dir,
                          threads=options.threads,
                          debug=options.debug)
    
        aligner_path = os.path.join(options.alignment_dir, "aligner.sh")
        aligner.run(open(aligner_path, 'w'))
        os.chmod(aligner_path, os.stat(aligner_path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        
        print "(2) running giza started  @", time.ctime()
        giza.run(os.path.join(options.alignment_dir, 'giza.log'))
        print "(2) running giza finished @", time.ctime()
    
    alignment = Alignment(cicada=cicada,
                          corpus=corpus,
                          cluster=prepare,
                          giza=giza,
                          alignment_dir=options.alignment_dir,
                          alignment=options.alignment,
                          threads=options.threads,
                          pbs=pbs,
                          debug=options.debug)

    if options.first_step <= 3 and options.last_step >= 3:
        print "(3) generate word alignment started  @", time.ctime()
        alignment.run()
        print "(3) generate word alignment finished @", time.ctime()

