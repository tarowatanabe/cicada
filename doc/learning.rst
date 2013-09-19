Learning
========

One of the important procedures in machine translation is to fit
parameters toward a given tuning data, which is very close to the
actual test data. Basically we supports three kinds of parameter
learning strategies:

k-best merging batch learning:
   In MERT, PRO or batch-MIRA learning, k-bests or pruned forests are
   generated in each iteration, and k-bests (or forests) are merged
   across iterations. The merged k-bests (or forests) are used as a
   training data for parameter optimization. We support MERT by
   ``cicada-mert.py``, and other objectives, like PRO or xBLEU by
   ``cicada-learn.py``.

Online learning: 
   In ``cicada-learn-online.py``, training data is split into
   mini-batch, and the mini-batch local k-best translations (or
   forests) are generated. Parameters are optimized using the local
   translations as a training data, and perform updates to the
   previously learned parameters.

Batch learning: 
   If non-local features, such as ngram language models, are not
   integrated in the model, we do not have to perform decoding in each
   iteration. ``cicada-maxent.py`` constructs forests, computes oracle
   forests (or uses gold-standard forests), and optimizes toward the
   oracles.

Reference Translations
----------------------

During learning, we need a set of reference translations which looks
like following:

::

   0 ||| first reference
   0 ||| second reference
   1 ||| first reference for the second input
   1 ||| second reference for the second input

See `doc/eval.rst` for details.

k-best merging batch learning
-----------------------------

In the k-best merging batch learning approach, k-bests are generated
in each iteration. This is managed by providing a template
configuration file. The training script, ``cicada-mert.py`` and
``cicada-learn.py`` replace following strings:

- ``${weights}``: Parameters for decoding.
- ``${kbest}``: k-best size for decoding. 0 implies forest output, not k-best.
- ``${file}``: Output file or directory.

MERT
````

MERT was a standard but not actively maintained.

.. code:: bash

  cicada-mert.py \
	  --srcset <source data> \
	  --refset <reference transaltion data> \
	  --config <configuration file>

By default, forest-MERT is performed. If the forest is too large, the
you can prune forests by appropriate operations or tweak beam
size. For details, see `doc/operation.rst`. If you prefer k-best MERT,
then, use ``--kbest <k-best-size>`` option. One of the benefits of
MERT is that it can use any evaluation metrics. You can use
``--scorer`` option to change the criterion. See `doc/eval.rst` for
details.

Others (Recommended)
````````````````````

Alternatively, following non-MERT is better in practice:

.. code:: bash

  cicada-learn.py \
	  --srcset <source data> \
	  --refset <reference transaltion data> \
	  --config <configuration file> \
	  --merge 

By default, we use xBLEU as our objective, which can be modified by
``--learn`` option. ``--merge`` implies k-best merging learning, and
without this option, the k-best translations generated in each
iteration is treated as a separate data set. Another important option
is ``--regularize-{l1,l2}`` which specify the strength of
optimization. Larger value implies less-overfit to training data, but
a smaller value implies fitting to the training data. By default, we
use no regularization, but usually, it is recommended to set
``--regularize-l2 1e-5`` which prefers fitting to the training data.

Online learning
---------------

When tuning data is very large, you can try online-learning which is
very fast in practice:

.. code:: bash

  cicada-learn-online.py \
	  --srcset <source data> \
	  --refset <reference transaltion data> \
	  --config <configuration file>

One of the major difference for setting online learning is the
configuration parameters. Since we do not interchange decoding and
optimization step, the configuration file should not contain any
variables, such as ``${weights}``, and remove ``output``
operation. The parameters for each operation, such as ``apply`` or
``prune`` are automatically set by decoder when no weights are
provided. By default, we use forests-based optimization if no
``--kbest`` options is provided. The use of ``--asynchronous`` flag
implies asynchronous parameter merging across workers which performed
the best in our internal studies.

Batch learning
--------------

There will be a situation when no non-local features, such as ngram
language models, are not integrated in the model, and uses the
features defined in a synchronous grammar. In this case, you can use a
simple batch learning:

.. code:: bash

  cicada-maxent.py \
	  --srcset <source data> \
	  --refset <reference transaltion data> \
	  --config <configuration file> \
	  --compose compose-cky

Which computes forest based from the source data using the
``compose-cky`` operation using the grammar specified in the given
configuration file. By default, we use `softmax` as an objective.

Precompute Forests
------------------

One of the resource demanding operations, both in terms of time and
memory, is "composition" operation. Since this composition step is not
affected by the parameters which are optimized during tuning, the
composed forests can be precomputed given a tuning data using the
following configuration, for instance, for SCFG:

::

   operation = compose-cky
   operation = output:directory=[output directory],forest=true

Then, the configuration file for tuning can avoid ``operation =
compose-cky`` and use ``input-forest =  true``  to load the
pre-composed forests. The ``--srcset`` option for tuning script can
use the ``[output directory]``.

Parallel Learning
-----------------

We support learning in parallel using either pthreads or MPI,
controlled by:

--threads        # of threads
--mpi            # of MPI jobs
--mpi-host       comma delimited list of hosts for MPI jobs
--mpi-host-file  host file for use with MPI jobs
--pbs            Run under PBS
--pbs-queue      PBS queue name

Remark that some objectives are implemented only by MPI or by
threads (sorry for this inconvenience!). If you see an error
message like "... is not supported by ...", you can try different
parallel learning strategies.
