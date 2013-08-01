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

``output`` is a part of operations, 


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



