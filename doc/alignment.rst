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
      src-trg.alignment.final.gz  (alignment model)
      src-trg.distortion.final.gz (distortion model)
      src-trg.fertility.final.gz  (fertility model)
      src-trg.lexicon.final.gz    (lexicon model)
  giza.trg-src/         (parameters for P(target | source))
  model/
      aligned.posterior-itg (word alignment)
      aligner.sh            (word aligner script)

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
the given bilingual data. The aligner script, ``model/aligner.sh`` can
be used to perform word alignment for held-out test data:

.. code:: bash

  model/aligner.sh \
	  --source <source test data> \
	  --target <target test data> \
	  --viterbi-source-target <alignment for source to target> \
	  --viterbi-target-source <alignment for target to source> \
	  --itg

which computes ITG constrained word alignment.

Details
-------

The supported alignment models are:

- IBM Model 1 [1]_ (``--iteration-model1 5``)
- HMM [2]_         (``--iteration-hmm 5``)
- IBM Model 4 [3]_ (``--iteration-model4 5``)

Two directions are integrated during the learning process either by:

- Native symmetric learning [2]_ (``--symmetric``)
- Posterior constrained learning [3]_ (``--symmetric`` and ``--posterior``, which are recommended)

The parameters are smoothed by:

- Native variational Bayes estimate (``--variational``, which is recommended)
- L0 regularization [5]_ (``--l0``, probably better than ``--variational`` but slow.)

After the parameter estimation, we can produce the final word
alignment by specifying ``--alignment`` option:

- Simple heuristics (``grow-diag-final`` etc.)
- Produce ITG constrained alignment from posteriors (``posterior-itg``, which is the default parameter)
- Produce one-to-one alignment using the Hungarian algorithm from
  posteriors (``posterior-max-match``).

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

