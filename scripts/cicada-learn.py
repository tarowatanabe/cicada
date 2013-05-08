#!/usr/bin/env python
#
#  Copyright(C) 2012-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
#

import threading
import multiprocessing

import time
import sys
import os, os.path
import string
import re
import subprocess
import cStringIO

from optparse import OptionParser, make_option

opt_parser = OptionParser(
    option_list=[
	
    # output directory/filename prefix
    make_option("--root-dir", default="", action="store", type="string",
                metavar="DIRECTORY", help="root directory for outputs"),
    make_option("--prefix", default="learn", action="store", type="string",
                metavar="PREFIX", help="prefix for outputs (default: %default)"),

    
    make_option("--devset", default="", action="store", type="string",
                metavar="FILE", help="training data"),
    make_option("--refset", default="", action="store", type="string",
                metavar="FILE", help="reference translations"),

    make_option("--config", default="", action="store", type="string",
                metavar="CONFIG", help="cicada config file"),

    make_option("--iteration", default=10, action="store", type="int",
                metavar="ITERATION", help="# of iterations (default: %default)"),
    make_option("--iteration-first", default=1, action="store", type="int",
                metavar="ITERATION", help="The first iteration (default: %default)"),
    make_option("--weights", default="", action="store", type="string",
                metavar="FILE", help="initial weights"),
    
    make_option("--C", default=1e-5, action="store", type="float",
                metavar="C", help="hyperparameter (default: %default)"),
    make_option("--regularize-l1", action="store_true",
                metavar="REGULARIZER", help="L1 regularization"),
    make_option("--regularize-l2", action="store_true",
                metavar="REGULARIZER", help="L2 regularization"),
    make_option("--scorer", default="bleu:order=4,exact=true", action="store", type="string",
                metavar="SCORER", help="scorer for oracle computation (default: %default)"),
    make_option("--scorer-cube", default=400, action="store", type="int",
                metavar="SIZE", help="cube size for oracle computation (default: %default)"),
    make_option("--learn", default="xbleu", action="store", type="string",
                metavar="LEARN", help="learning algorithms from [softmax, svm, linear, pegasos, mira, cw, arow, nherd, cp, mcp, xbleu] (default: %default)"),
    make_option("--learn-options", default="", action="store", type="string",
                metavar="OPTION", help="additional learning options"),
    
    make_option("--kbest", default=0, action="store", type="int",
                metavar="KBEST", help="kbest size (default: %default)"),
    make_option("--forest", default=None, action="store_true",
                help="forest based learning"),
    make_option("--merge", default=None, action="store_true",
                help="perform kbest merging"),
    make_option("--interpolate", default=0.0, action="store", type="float",
                help="perform weights interpolation"),
        
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
    make_option("--mpi-options", default="", action="store", type="string",
                metavar="OPTION", help="additional MPI options"),
    make_option("--threads", default=1, action="store", type="int",
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
        self.qsub = 'qsub'
        
        # how to find binary location...?

    def run(self, command="", name="name", memory=0.0, mpi=None, threads=1, logfile=None):

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

            if os.environ.has_key('TMPDIR_SPEC'):
                prefix += ' -x TMPDIR_SPEC'
            if os.environ.has_key('LD_LIBRARY_PATH'):
                prefix += ' -x LD_LIBRARY_PATH'
            if os.environ.has_key('DYLD_LIBRARY_PATH'):
                prefix += ' -x DYLD_LIBRARY_PATH'
            
            prefix += ' ' + mpi.options
            prefix += ' '
        
        suffix = ''
        if logfile:
            suffix = " 2> %s" %(logfile)

        pipe.write(prefix + command + suffix + '\n')
        
        pipe.close()
        popen.wait()
        
class MPI:
    def __init__(self, dir="", hosts="", hosts_file="", number=0, options=""):
        
	self.dir = dir
	self.hosts = hosts
        self.hosts_file = hosts_file
        self.number = number
        self.options = options
	
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
                if not os.path.exists(prog) or os.path.isdir(prog):
                    prog = os.path.join(self.bindir, binprog)
                    if not os.path.exists(prog) or os.path.isdir(prog):
                        raise ValueError, prog + " does not exist at " + self.bindir
                    
                setattr(self, binprog, prog)
            else:
                setattr(self, binprog, binprog)
                
    def run(self, command, logfile=None):
        mpirun = self.mpirun
        if self.number > 0:
            mpirun += ' --np %d' %(self.number)
        if self.hosts:
            mpirun += ' --host %s' %(self.hosts)
        elif self.hosts_file:
            mpirun += ' --hostfile "%s"' %(self.hosts_file)

        if os.environ.has_key('TMPDIR_SPEC'):
            mpirun += ' -x TMPDIR_SPEC'
        if os.environ.has_key('LD_LIBRARY_PATH'):
            mpirun += ' -x LD_LIBRARY_PATH'
        if os.environ.has_key('DYLD_LIBRARY_PATH'):
            mpirun += ' -x DYLD_LIBRARY_PATH'

        mpirun += ' ' + self.options
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
            self.pbs.run(str(command), name=name, memory=memory, mpi=self.mpi, logfile=logfile)
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
	
        for binprog in ('cicada',
                        'cicada_mpi',
                        'cicada_learn',
                        'cicada_learn_mpi',
                        'cicada_learn_kbest',
                        'cicada_learn_kbest_mpi',
                        'cicada_oracle',
                        'cicada_oracle_mpi',
                        'cicada_oracle_kbest',
                        'cicada_oracle_kbest_mpi',
                        'cicada_eval',
                        'cicada_filter_config',
                        'cicada_filter_weights',):
	    
	    for bindir in bindirs:
		prog = os.path.join(bindir, binprog)
                
                if not os.path.exists(prog): continue
                if os.path.isdir(prog): continue
                
                setattr(self, binprog, prog)
                break

	    if not hasattr(self, binprog):
		raise ValueError, binprog + ' does not exist'

def learn_algorithms(command):
    pattern = re.compile(r"\s+--learn-(\S+)\s*")

    algs = []
    
    for line in cStringIO.StringIO(subprocess.check_output([command, "--help"])):
        result = pattern.search(line)
        if result:
            algs.append(result.group(1))
    return algs
    

if __name__ == '__main__':
    (options, args) = opt_parser.parse_args()
    
    ### dump to stderr
    stdout = sys.stdout
    sys.stdout = sys.stderr

    ### config
    if not os.path.exists(options.config):
        raise ValueError, "no config file: %s" %(options.config)
    
    ### root-dir
    if options.root_dir:
        if not os.path.exists(options.root_dir):
            os.makedirs(options.root_dir)

    ### regularizer
    regularizer = "--regularize-l2"
    if options.regularize_l1 and options.regularize_l2:
        raise ValueError, "we do not support both L1 and L2 regularizers"
    if options.regularize_l1:
        regularizer = "--regularize-l1"
    if options.regularize_l2:
        regularizer = "--regularize-l2"

    ### cicada
    cicada = CICADA(dir=options.cicada_dir)
    
    ### MPI
    mpi = None
    if options.mpi_host or options.mpi_host_file or options.mpi > 0:
        mpi = MPI(dir=options.mpi_dir,
                  hosts=options.mpi_host,
                  hosts_file=options.mpi_host_file,
                  number=options.mpi,
                  options=options.mpi_options)
    
    ### PBS
    pbs = None
    if options.pbs:
        pbs = PBS(queue=options.pbs_queue)
    
    ### QSUB
    qsub = QSub(mpi=mpi, pbs=pbs)
    
    ### iterations
    if options.iteration_first <= 0:
        options.iteration_first = 1
    if options.iteration_first > options.iteration:
        raise ValueError, "invalid iterations"

    ### defaults to forest...
    if not options.forest and options.kbest <= 0:
        options.forest = 1
        options.kbest = 0
    if options.forest and options.kbest > 0:
        raise ValueError, "forest-mode or kbest-mode?"

    ## kbest or forest oracle..
    cicada_oracle     = cicada.cicada_oracle
    cicada_oracle_mpi = cicada.cicada_oracle_mpi
    if not options.forest:
        cicada_oracle     = cicada.cicada_oracle_kbest
        cicada_oracle_mpi = cicada.cicada_oracle_kbest_mpi

    learn_mpi        = None
    cicada_learn     = None
    cicada_learn_mpi = None
    
    if options.forest:
        learn_forest     = learn_algorithms(cicada.cicada_learn)
        learn_forest_mpi = learn_algorithms(cicada.cicada_learn_mpi)
        
        if options.learn not in learn_forest or options.learn not in learn_forest_mpi:
            raise ValueError, "%s is not supported by forest learner" %(options.learn)

        if not mpi and options.learn not in learn_forest:
            raise ValueError, "%s is not supported by non-mpi-forest learner" %(options.learn)

        cicada_learn     = cicada.cicada_learn
        cicada_learn_mpi = cicada.cicada_learn_mpi
        
        if mpi and options.learn in learn_forest_mpi:
            learn_mpi = 1

    else:
        learn_kbest     = learn_algorithms(cicada.cicada_learn_kbest)
        learn_kbest_mpi = learn_algorithms(cicada.cicada_learn_kbest_mpi)

        if options.learn not in learn_kbest or options.learn not in learn_kbest_mpi:
            raise ValueError, "learner %s is not supported by kbest learner" %(options.learn)
        
        if not mpi and options.learn not in learn_kbest:
            raise ValueError, "%s is not supported by non-mpi-kbest learner" %(options.learn)

        cicada_learn     = cicada.cicada_learn_kbest
        cicada_learn_mpi = cicada.cicada_learn_kbest_mpi

        if mpi and options.learn in learn_kbest_mpi:
            learn_mpi = 1

    interpolate = None
    if options.interpolate > 0.0 and options.interpolate < 1.0:
        interpolate = 1

    weights_config = 'weights-one=true'
    if options.weights:
        if not os.path.exists(options.weights):
            raise ValueError, "no initial weights %s" %(options.weights)
        
        weights_config = "weights=%s" %(optins.weights)
    else:
        weights_file = os.path.join(options.root_dir, options.prefix + ".0.weights")
        
        open(weights_file, 'w').close()

        weights_config = "weights=%s" %(weights_file)
    
    weiset = []
    tstset = []
    orcset = []
    
    for iter in range(1, options.iteration_first):
        prefix = options.prefix + ".%d" %(iter)

        weights = os.path.join(options.root_dir, prefix + ".weights")
        decoded = os.path.join(options.root_dir, prefix + ".forest")
        oracle  = os.path.join(options.root_dir, prefix + ".forest.oracle")
        if not options.forest:
            decoded = os.path.join(options.root_dir, prefix + ".kbest")
            oracle  = os.path.join(options.root_dir, prefix + ".kbest.oracle")
        
        weiset.append(weights)
        tstset.append(decoded)
        orcset.append(oracle)

    for iter in range(options.iteration_first, options.iteration+1):
        print "iteration: %d" %(iter)
        
        ## setup output files
        prefix = options.prefix + ".%d" %(iter)
        
        weights       = os.path.join(options.root_dir, prefix + ".weights")
        weights_learn = os.path.join(options.root_dir, prefix + ".weights.learn")
        if len(weiset) > 0:
            weights_config = 'weights=%s' %(weiset[-1])
        
        decoded = os.path.join(options.root_dir, prefix + ".forest")
        oracle  = os.path.join(options.root_dir, prefix + ".forest.oracle")
        if not options.forest:
            decoded = os.path.join(options.root_dir, prefix + ".kbest")
            oracle  = os.path.join(options.root_dir, prefix + ".kbest.oracle")
        onebest = os.path.join(options.root_dir, prefix + ".1best")
        
        weiset.append(weights)
        tstset.append(decoded)
        orcset.append(oracle) 
        
        config = os.path.join(options.root_dir, prefix + ".config")
        mteval = os.path.join(options.root_dir, prefix + ".eval")
        
        print "generate config file %s @ %s" %(config, time.ctime())
        
        qsub.run(Program(cicada.cicada_filter_config,
                         Option('--weights', weights_config),
                         Option('--kbest', options.kbest),
                         Option('--file', "directory=%s" %(decoded)),
                         Option('--input', Quoted(options.config)),
                         Option('--output', Quoted(config))),
                 name="config")
        
        print "decode %s @ %s" %(decoded, time.ctime())
        
        if mpi:
            qsub.mpirun(Program(cicada.cicada_mpi,
                                Option('--input', Quoted(options.devset)),
                                Option('--config', Quoted(config)),
                                Option('--debug')),
                        name="decode",
                        memory=options.max_malloc,
                        threads=options.threads,
                        logfile=Quoted(decoded+'.log'))
        else:
            qsub.run(Program(cicada.cicada,
                             Option('--input', Quoted(options.devset)),
                             Option('--config', Quoted(config)),
                             Option('--threads', options.threads),
                             Option('--debug')),
                     name="decode",
                     memory=options.max_malloc,
                     threads=options.threads,
                     logfile=Quoted(decoded+'.log'))
        
        if options.forest:
            print "1best %s @ %s" %(onebest, time.ctime())
                
            if mpi:
                qsub.mpirun(Program(cicada.cicada_mpi,
                                    Option('--input', Quoted(decoded)),
                                    Option('--input-forest'),
                                    Option('--input-directory'),
                                    Option('--operation', 'output:kbest=1,%s,file=%s' %(weights_config, onebest)),
                                    Option('--debug')),
                            name="onebest",
                            memory=options.max_malloc,
                            threads=options.threads,
                            logfile=Quoted(onebest+'.log'))
            else:
                qsub.run(Program(cicada.cicada,
                                 Option('--input', Quoted(decoded)),
                                 Option('--input-forest'),
                                 Option('--input-directory'),
                                 Option('--operation', 'output:kbest=1,%s,file=%s' %(weights_config, onebest)),
                                 Option('--threads', options.threads),
                                 Option('--debug')),
                         name="onebest",
                         memory=options.max_malloc,
                         threads=options.threads,
                         logfile=Quoted(onebest+'.log'))
            
            print "evaluate %s @ %s" %(mteval, time.ctime())
            
            qsub.run(Program(cicada.cicada_eval,
                             Option('--refset', Quoted(options.refset)),
                             Option('--tstset', Quoted(onebest)),
                             Option('--output', Quoted(mteval)),
                             Option('--scorer', options.scorer)),
                     name="evaluate")
        else:
            print "evaluate %s @ %s" %(mteval, time.ctime())
            
            qsub.run(Program(cicada.cicada_eval,
                             Option('--refset', Quoted(options.refset)),
                             Option('--tstset', Quoted(decoded)),
                             Option('--output', Quoted(mteval)),
                             Option('--scorer', options.scorer)),
                     name="evaluate")
        
        print "oracle %s @ %s" %(oracle, time.ctime())
        
        ### set up the oracle to compute...
        oracle_tstset = Option('--tstset', Quoted(decoded))
        if options.merge:
            oracle_tstset = Option('--tstset', ' '.join(map(lambda x: str(Quoted(x)), tstset)))

        oracle_scorer_cube = ''
        oracle_forest = ''
        if options.forest:
            oracle_scorer_cube = Option('--scorer-cube', options.scorer_cube)
            oracle_forest = Option('--forest')
                    
        if mpi:
            qsub.mpirun(Program(cicada_oracle_mpi,
                                Option('--refset', Quoted(options.refset)),
                                oracle_tstset,
                                Option('--output', Quoted(oracle)),
                                Option('--scorer', options.scorer),
                                oracle_scorer_cube,
                                oracle_forest,
                                Option('--directory'),
                                Option('--debug')),
                        name="oracle",
                        memory=options.max_malloc,
                        threads=options.threads,
                        logfile=Quoted(oracle+'.log'))
        else:
            qsub.run(Program(cicada_oracle,
                             Option('--refset', Quoted(options.refset)),
                             oracle_tstset,
                             Option('--output', Quoted(oracle)),
                             Option('--scorer', options.scorer),
                             oracle_cube_size,
                             oracle_forest,
                             Option('--directory'),
                             Option('--threads', options.threads),
                             Option('--debug'),),
                     name="oracle",
                     memory=options.max_malloc,
                     threads=options.threads,
                     logfile=Quoted(oracle+'.log'))
        
        print "learn %s @ %s" %(weights, time.ctime())

        ### training data, oracle data
        learn_input  = Option('--input', ' '.join(map(lambda x: str(Quoted(x)), tstset)))
        learn_oracle = Option('--oracle', ' '.join(map(lambda x: str(Quoted(x)), orcset)))
        learn_unite  = ''
        if options.merge:
            learn_oracle = Option('--oracle', Quoted(oracle))
            learn_unite  = Option('--unite')
        
        learn_weights = ''
        if len(weiset) > 1:
            learn_weights = Option('--weights', Quoted(weiset[-2]))
            
        # if we interpolate, first, generate weights_learn then, dump...
        learn_output = Option('--output', Quoted(weights))
        if interpolate and len(weiset) > 1:
            learn_output = Option('--output', Quoted(weights_learn))

        learn_algorithm = Option('--learn-' + options.learn)

        if learn_mpi and mpi:
            qsub.mpirun(Program(cicada_learn_mpi,
                                learn_input,
                                learn_oracle,
                                Option('--refset', options.refset),
                                Option('--scorer', options.scorer),
                                learn_unite,
                                learn_output,
                                learn_weights,
                                learn_algorithm,
                                options.learn_options,
                                Option('--C', options.C),
                                Option(regularizer),
                                Option('--debug', 2),),
                        name="learn",
                        memory=options.max_malloc,
                        threads=options.threads,
                        logfile=Quoted(weights+'.log'))
        else:
            qsub.run(Program(cicada_learn,
                             learn_input,
                             learn_oracle,
                             Option('--refset', options.refset),
                             Option('--scorer', options.scorer),
                             learn_unite,
                             learn_output,
                             learn_weights,
                             learn_algorithm,
                             options.learn_options,
                             Option('--C', options.C),
                             Option(regularizer),
                             Option('--threads', options.threads),
                             Option('--debug', 2),),
                     name="learn",
                     memory=options.max_malloc,
                     threads=options.threads,
                     logfile=Quoted(weights+'.log'))
            
        if interpolate:
            print "interpolate %s @ %s" %(weights, time.ctime())
            
            qsub.run(Program(cicada.cicada_filter_weights,
                             Option('--output', Quoted(weights)),
                             Option(weiset[-2] + ":scale=%g" %(1.0 - options.interpolate)),
                             Option(weights_learn + ":scale=%g" %(options.interpolate))),
                     name="interpolate")
