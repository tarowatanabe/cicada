===================
How to Build cicada
===================

Get the cutting-edge source code from `github.com <http://github.com/tarowatanabe/cicada>`_:

.. code:: bash

  git clone https://github.com/tarowatanabe/cicada.git

Or, grab the stable tar archive from `cicada <http://www2.nict.go.jp/univ-com/multi_trans/cicada>`_.

Compile
-------

.. code:: bash

   ./autogen.sh (required when you get the code by git clone)
   ./configure
   make
   make install (optional)

You can set several options. For details see the requirement section.
::

  --enable-snappy         enable snappy
  --enable-msgpack        enable msgpack
  --enable-jemalloc       enable jemalloc
  --enable-tcmalloc       enable tcmalloc
  --enable-profiler       enable profiling via google's libprofiler
  --enable-static-boost   Prefer the static boost libraries over the shared
                          ones [no]

  --with-kenlm-order=ORDER
                          kenlm's max order [default=6]
  --with-zlib=DIR         zlib in DIR
  --with-bzlib=DIR        bzlib in DIR
  --with-lzma=DIR         lzma in DIR
  --with-blas=DIR         blas in DIR
  --with-snappy=DIR       snappy in DIR
  --with-msgpack=DIR      msgpack in DIR
  --with-jemalloc=DIR     jemalloc in DIR
  --with-tcmalloc=DIR     tcmalloc in DIR
  --with-profiler=DIR     profiler in DIR
  --with-boost=DIR        prefix of Boost 1.42 [guess]

Requirements
------------

- expgram: http://www2.nict.go.jp/univ-com/multi_trans/expgram

  This is our ngram toolkit for ngram language model training/indexing.
  Although we support `kenlm` which can estimate ngram language model
  and perform indexing, expgram is smaller and performs better rest-cost
  estimation for CKY-style decoding, but at the cost of slower speed.

- Boost library: http://www.boost.org/

  The minimum requirement is boost version 1.42. Prior to this
  version, there were a couple of serious bugs which prevent us from
  running correctly.

- MPI (Open MPI): http://www.open-mpi.org/

  We strongly recommend open-mpi since it is regularly tested.
  The MPI libraries are automatically detected by the `configure`
  script by finding either `mpic++`, `mpicxx` or `mpiCC`. Thus, those
  mpi specific compilers should be on the executable path.

- ICU: http://site.icu-project.org/
   
  The `configure` script relies on `icu-config` installed by the ICU
  library. Thus, `icu-config` must be in the executable path.

- Optional:

   + snappy: http://code.google.com/p/snappy/

   + msgpack: http://msgpack.org

     NOTE: msgpack-0.5.7 has a bug in which deletion may be called twice!
     Grab the latest copy from https://github.com/msgpack/msgpack-c

   + srilm:   Language model training (http://www-speech.sri.com/projects/srilm/)

   + GIZA++:  alignment model training (http://code.google.com/p/giza-pp/)

     You can uses moses (http://www.statmt.org/moses/) for alignment training with GIZA++
     or, berkeley aligner (http://code.google.com/p/berkeleyaligner/),
     postcat (http://www.seas.upenn.edu/~strctlrn/CAT/CAT.html).

     Remark that cicada can "align words" using symmetized training of berkeley aligner and/or posterior
     constrained training of postcat with parameter smoothing via naive-Bayes or L0-norm.

   + For better memory management:

     * gperftools: http://code.google.com/p/gperftools/
     * jemalloc: http://www.canonware.com/jemalloc/

     For Linux, you should install one of them for better memory performance
     and to measure how many bytes malloced, since mallinfo is
     "broken" for memory more than 4GB.
     They are configured by `--with-{jemalloc,tcmalloc}` and should be
     enabled using `--enable-{jemalloc,tcmalloc}`

   + docutils: http://docutils.sourceforge.net

     Manpages are written in reStructuredText format, and if you want
     manpages, you need to install docutils.
