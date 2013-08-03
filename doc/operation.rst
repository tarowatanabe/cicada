Operations
==========

``cicada`` and its MPI version ``cicada_mpi`` process inputs in
parallel, and each input is processed by a sequence of operations
specified by the ``--operation`` options.

Input
-----

Input to cicada is controlled by the ``--input-*`` flags, and each
line may consists of various components split by a delimiter ``|||``:

::

  [id |||] [lattice or sentence] [||| forest] [||| span] [||| alignment] [||| dependency] [||| sentence]*


The id is optional, and if an id is not supplied here, cicada assumes
a numerical order starting from 0.
By default, cicada assumes sentence input, but if ``--input-lattice``
option is specified, lattice is assumed, instead.
It is followed by optional entries, `forest` (``--input-forest``),
`span` (``--input-span``), `alignment` (``--input-alignment``),
`dependency` (``--input-dependency``) and a set of `sentence` (``--input-bitext``).

The input can be stdin (when ``--input`` is `-`) or a file with each
line consisting of an input described above.
Or, a directory can be specified. When ``--input-directory`` is
specified, we always try a directory input, and if the input file
specified by the ``--input`` option is not a directory, cicada will quit
with error.
The directory should contain a file ``id.gz`` with each file contains
one line of an input, and id specifies the id of the input.

Output
------

``output`` is a part of operations, and multiple ``output`` instances
are supported, but we do not support output to different directories
or different files, but they are output to the same directory or file
pointed to by the ``directory`` or ``file`` options.
``kbest`` specifies the k-best derivations using the ``weights`` or
``weights-one``. Different types of yields are supported by the
``yield`` option. ``unique`` option controls whether to output unique
k-bests or allow duplicates. When ``kbest`` is not set or set to zero,
we will output forest by default, but you can out put both or either
of ``forest`` and/or ``lattice``. Other outputs, like alignment is
also supported by ``alignment`` option. The ``graphviz`` output ouput
forest or lattice in a dot format for visualization by graphviz.

:: 

  output: kbest or hypergraph output
        kbest=<kbest size> zero for hypergraph output (default)
        insertion-prefix=<prefix attatched to inserted word>
        unique=[true|false] unique translation
        weights=weight file for feature
        weights-one=[true|false] one initialize weight
        yield=[sentence|string|terminal-pos|derivation|tree|treebank|graphviz|alignment|span]
	yield  for kbest
        graphviz=[true|false] dump in graphviz format
        debinarize=[true|false] debinarize k-best trees
        statistics=[true|false] dump various statistics (size etc.)
        lattice=[true|false] dump lattice
        forest=[true|false] dump forest
	span=[true|false] dump spans
        alignment=[true|false] dump alignment
        dependency=[true|false] dump dependency
        bitext=[true|false] dump bitext
        no-id=[true|false] do not output id
        directory=directory for output
        file=file for output


Operations
----------

The full list of available options and their options can be listed by
``cicada --operation-list``. Here is a quick summary of the
operations:

apply:
    feature application
  
binarize:
   perform binarization (monolingual tree)

debinarize:
   de-binarize forest

compose-:
   composition

parse-:
   parse
   

generate-earley: re-generation from tree
	depth: depth of rule pattern (= vertial Markovization + 1. <= 0 for infinity)
	width: width of rule pattern (= horitonzal Markovization. < 0 for infinity)

intersect: compute intersection
	lattice=[true|false] intersect with lattice
	target=[true|false] intersect with one of target


expand-ngram: expand hypergrpah by ngram
	order=<ngram order>

expected-ngram: expected ngram computation
	order=<ngram order>
	bos-eos=[true|false] include <s> and </s> in ngrams
	weights=weight file for feature
	weights-one=[true|false] one initialized weight
	scale=scaling for score
normalize: feature value normalizer
	prefix=feature name prefix

push-:
   push-weights

remove-:
    remove symbols

sort-:
   sort forest
   
span-forest:
    annotate terminal span

clear:
    clear data structure

verify:
    verify hypergrpah

viterbi:
    compute viterbi tree



