Word alignment
==============

The cicada toolkit can align words given bilingual data. Basically,
you can perform word alignment by running the script
`cicada-alignment.py` with ``--f`` for the source side and ``--e`` for
the target side.

.. code:: bash

  cicada-alignment.py --f <source data> --e <target data>

This will result in four directories and files:

::

  corpus/
      src.vcb.classes     (word classes for the source side)
      trg.vcb.classes     (word classes for the target side)
  giza.src-trg/         (parameters for P(source | target))
      src-trg.A3.final.gz         (Viterbi alignment)
      src-trg.posterior.final.gz  (posterior for P(source | target))
      src-trg.alignment.final.gz  (alignment model for HMM)
      src-trg.distortion.final.gz (distortion model for Model4)
      src-trg.fertility.final.gz  (fertility model for Model4)
      src-trg.lexicon.final.gz    (lexicon model for Model1/HMM/Model4)
  giza.trg-src/         (parameters for P(target | source))
  model/
      aligned.posterior-itg (word alignment)
      aligner.sh            (word aligner script. See **Forced Alignment** section)

Word alignment training is performed in three steps:

1. Source and target words are clustered, and clustering
   results are placed in the `corpus` directory (or specified by
   ``--corpus-dir`` option).
2. Word alignment models are learned and Viterbi alignment is
   computed, then, all the parameters are put on two directories,
   `giza.src-trg`  and `giza.trg-src` directories (or specified by
   ``--giza-f2e`` and ``--giza-e2f`` options).
3. Two models in opposite directions are combined, and put in `model`
   directory (or specified by ``--model-dir``).

The file ``model/aligned.posterior-itg`` is the final alignment for
the given bilingual data.

Details
-------

The supported alignment models are:

- IBM Model 1 [1]_ (``--iteration-model1 5``)
- HMM [2]_         (``--iteration-hmm 5``)
- IBM Model 4 [1]_ (``--iteration-model4 5``)

which are controlled by specifying the number of iterations. If you
want to train only HMM, then, set ``--iteration-model4`` to ``0``.

Two directions are integrated during the learning process either by:

- Naive symmetric learning [2]_ (``--symmetric``)
- Posterior regularized learning [3]_ (``--symmetric`` and ``--posterior``, which are recommended)

The parameters are smoothed by:

- Native variational Bayes estimate (``--variational``, which is recommended)
- L0 regularization [5]_ (``--l0``, probably better than ``--variational`` but slow.)

After the parameter estimation, we can produce the final word
alignment by specifying ``--alignment`` option:

- Simple combination heuristics from two Viterbi alignments of two
  directions (``grow-diag-final`` etc.)
- Produce word alignment which are higher than a certain threshold
  given the posterior probabilities in two directions
  (``posterior-0.1`` etc.)
- Produce ITG constrained alignment from the combined posteriors
  (``posterior-itg``, which is the default parameter)
- Produce one-to-one alignment using the Hungarian algorithm from
  posteriors (``posterior-max-match``).

Forced Alignment
----------------

The aligner script, ``model/aligner.sh`` can be used to perform word
alignment for held-out test data:

.. code:: bash

  model/aligner.sh \
	  --source <source test data> \
	  --target <target test data> \
	  --viterbi-source-target <alignment for P(target | source)> \
	  --viterbi-target-source <alignment for P(source | target)> \
	  --itg

which computes ITG constrained word alignment, for instance.
The word alignment can be heuristically combined by first computing
Viterbi alignment in two directions:

.. code:: bash

  model/aligner.sh \
	  --source <source test data> \
	  --target <target test data> \
	  --viterbi-source-target <alignment for P(target | source)> \
	  --viterbi-target-source <alignment for P(source | target)>

Then, use ``cicada_alignment`` to merge:

.. code:: bash

  cicada_alignment \
	  --source-target <alignment for P(target | source)> \
	  --target-source <alignment for P(source | target)> \
	  --grow \
	  --diag \
	  --final-and

which applies ``grow-diag-final-and`` heuristic. You can also try
``--itg`` heuristic or ``--max-match`` to compute one-to-one alignment
using Hungarian algorithm. Alternatively, combined word alignment can
be estimated by first generating posteriors:

.. code:: bash

  model/aligner.sh \
	  --source <source test data> \
	  --target <target test data> \
	  --posterior-source-target <posteriors for P(target | source)> \
	  --posterior-target-source <posteriors for P(source | target)>

The posterior probabilities are encoded by a variant of base64 for
faster output:

::

 ((BAAAAAAAAAAA, BAAAAAAAAAAA, BAAAAAAAAAAA, BAAAAAAAAAAA, BAAAAAAAAAAA, BAAAAAAAAAAA, BAAAAAAAAAAA),
  (BOt6deP8ZNjY, BxCaxPXuuVj8, BzUiQOgdyuD8, B8P8ytb+X7D8, BkJgHdwmqgz8, BFPk16J3ibSw, B+vJ6/wDegDU),
  (Bn73fUUsEFTU, Bee1KAQ2y6yI, BnQAlHbua4iI, BTiE5AoyPtDU, BAAAAAAAA8D8, B3nsVxMxdmTU, BRycFC+VWays),
  (BPKvgtOMSwD4, BcUjVdQHvdTU, B6JMjnTJlrz4, BkJgHdwmqgz8, Bq/NX3lGx7z8, BSw5D7lK0yDU, BbQNS59z6Fz4),
  (B/1Rz1KsjejU, Bid6TCsHARSw, Bgt9eQq6vQCw, BIgSwLz5sJDU, BAAAAAAAA8D8, BHl5G3plkvTU, BKizUASpb+TQ),
  (BPy4hl7Gf7T8, BUSLB/IT0ezY, BHWMehlVFnT8, BOqZCBC+0nj8, B8lfLAERPJjY, B64eBzGTZuDY, ByC55kUgQkD8),
  (BWBsRiKhJcTU, B+jlodVAA1Sw, B1yI+lLOX2DU, BqrfGuOPsFzU, BGMwuXgo5DTY, BAAAAAAAA8D8, BblJQxB4PtjU),
  (BZ8syrJdB0D4, BZVlC9YIBSzY, BUAc0pdEqlDU, BqmaIdrR0Fj4, BZQcbJifTXjU, BuObcWrOq6TU, BQ4SA3vf/7z8))

The column indicates the target (or foreign) sentence position and
each row corresponds to the source (or English) sentence position
including NULL word. Each value is prefixed by 'B' followed by 11 byte
base64 encoded string. To uncover double value, for instance in
python, you can try:

.. code:: python

  struct.unpack('d', base64.b64decode(('BOt6deP8ZNjY' + 'A')[1:])[:-1])

in which we will pad extra byte ('A') and strip off the prefix ('B'),
then strip off the final byte of decoded byte string.

Then, use a threshold to combine them:

.. code:: bash

  cicada_alignment \
	  --source-target <posterior for P(target | source)> \
	  --target-source <posterior for P(source | target)> \
	  --posterior \
	  --posterior-threshold 0.2

which annotate a word pair as aligned when the square root of the
product of its posterior probabilities in two directions is higher
than 0.2 [2]_. As in Viterbi alignment combination, you can try
``--itg`` to estimate ITG constrained alignment or ``--max-match``
to apply Hungarian algorithm to compute one-to-one alignment.

Visualization
-------------

Word alignment can be visualized by `cicada_filter_alignment`:

.. code:: bash

  cicada_fiter_alignment 
	--source <source file>
	--target <target file>
	--alignment  <alignment file>
	--alignment2 <secondary alignment file> (optional)
	--inverse (inverse alignment, optional)
	--visualize (required for visualization!)

where: 

- Blue points indicate intersection.
- Green points indicate word alignment only in the primary alignment.
- Yello points indicate word alignment only in the secondary alignment.


References
----------

.. [1]	 Peter F. Brown, Vincent J. Della Pietra, Stephen A. Della
	 Pietra, and Robert L. Mercer. The mathematics of statistical
	 machine translation: parameter estimation. Comput. Linguist.,
	 19:263-311, June 1993.

.. [2]	 Percy Liang, Ben Taskar, and Dan Klein. Alignment by
	 agreement. In Proceedings of the Human Language Technology
	 Conference of the NAACL, Main Conference, pages 104-111, New
	 York City, USA, June 2006. Association for Computational
	 Linguistics.

.. [3]	 Kuzman Ganchev, João V. Graça, and Ben Taskar. Better
	 alignments = better translations? In Proceedings of ACL-08:
	 HLT, pages 986-993, Columbus, Ohio, June 2008. Association
	 for Computational Linguistics.

.. [4]	 Arne Mauser, Saša Hasan, and Hermann Ney. Extending
	 statistical machine translation with discriminative and
	 trigger-based lexicon models. In Proceedings of the 2009
	 Conference on Empirical Methods in Natural Language
	 Processing, pages 210-218, Singapore,
	 August 2009. Association for Computational Linguistics.

.. [5]	 Ashish Vaswani, Liang Huang, and David Chiang. Smaller
	 alignment models for better translations: Unsupervised word
	 alignment with the l0-norm. In Proceedings of the 50th Annual
	 Meeting of the Association for Computational Linguistics
	 (Volume 1: Long Papers), pages 311-319, Jeju Island, Korea,
	 July 2012. Association for Computational Linguistics.

