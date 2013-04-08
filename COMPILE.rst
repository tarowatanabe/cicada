Requirements

expgram: ngram language model training/indexing

Others:
Boost library     (http://www.boost.org/)
MPI (Open MPI)    (http://www.open-mpi.org/)
    REMARK: Under Linux, it is recommended not to use memory managers in open-mpi, which may conflict with your
    	    favorite memory managers, such as jemalloc/tcmalloc. To disable this, for instance, you
	    can edit "openmpi-mca-params.conf" and add(or disable) mca by, "memory = ^ptmalloc2"
	    For details, consult open-mpi documentation.

ICU               (http://site.icu-project.org/)

Optional:
	msgpack: http://msgpack.org
		 NOTE: msgpack-0.5.7 has a bug in which deletion may be called twice!
		       Grab the latest copy from https://github.com/msgpack/msgpack-c

	srilm:   Language model training (http://www-speech.sri.com/projects/srilm/)
	         For indexing, you still need "expgram"

	GIZA++:  alignment model training (http://code.google.com/p/giza-pp/)
	         or, you can uses moses (http://www.statmt.org/moses/) for alignment training with GIZA++
	         or, try berkeley aligner (http://code.google.com/p/berkeleyaligner/)
		 or, try postcat (http://www.seas.upenn.edu/~strctlrn/CAT/CAT.html).
		 
		 Remark that cicada can "align words" using symmetized training of berkeley aligner and/or posterior
		 constrained training of postcat with parameter smoothing via naive-Bayes or L0-norm.
	
	For better memory management:

	gperftools (http://code.google.com/p/gperftools/)
	jemalloc  (http://www.canonware.com/jemalloc/)

	   For Linux, you should install one of them for better memory performance
	   and to measure how many bytes malloced, since mallinfo is "broken" for memory more than 4GB.

