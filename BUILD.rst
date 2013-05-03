Requirements

expgram: ngram language model training/indexing

Others:
Boost library     (http://www.boost.org/)
MPI (Open MPI)    (http://www.open-mpi.org/)
   We strongly recommend open-mpi since it is regularly tested.

   If you encounter a problem like, "mca_btl_tcp_frag_recv: readv
   failed: Connection timed out" then, try this patch:

*** ompi/mca/btl/tcp/btl_tcp_frag.c.org	2012-04-03 23:30:11.000000000 +0900
--- ompi/mca/btl/tcp/btl_tcp_frag.c	2013-05-02 10:43:21.571867286 +0900
***************
*** 201,206 ****
--- 201,207 ----
  	switch(opal_socket_errno) {
  	case EINTR:
  	    continue;
+ 	case ETIMEDOUT:
  	case EWOULDBLOCK:
  	    return false;
  	case EFAULT:

   This will force open-mpi to re-reading buffer again, even after
   timeout.


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

