Evaluation
==========

Following metrics are implemented:

- BLEU: IBM BLEU [1]_
- BLEUS: smoothed BLEU+1 metric [6]_
- CDER: CDer (wer + shift in one side) [7]_
- WER: word error rate [2]_
- InvWER: Inversion word error rate [8]_
- PER: position independent word error rate [3]_
- TER: translation error rate [4]_
- RIBES: RIBES [5]_
- SB: skip bigram  [6]_
- WLCS: (weighted) longest common subsequence  [6]_
- SK: string kernel

For the details of available options, see ``cicada --scorer-list``.

Evaluator
---------

You can use `cicada_eval` to score your translations:

.. code:: bash

  cicada_eval \
	  --tstset <decoder k-best output> \
	  --refset <reference data>

The ``--tstset`` can be either plain text or directory which contains
k-best translations. Bootstrap resampling [9]_ (``--bootstrap``) and
sign test [10]_ (``--sign``) are also implemented to measure the
significance.

The scorer can be set by ``--scorer`` option. For instance, BLEU can
be measured by:

.. code:: bash

  --scorer bleu:order=4

If you want to use RIBES, then:

.. code:: bash

  --scorer ribes

Two or more metrics can be linearly combined:

.. code:: bash

  --scorer 'combined:metric="bleu:order=4",weight=0.5,metric=ribes,weight=0.5'

The combination assume **reward** score, not loss score. Thus, if you
want to integrate error metrics, such as TER, the weights should be
negatives. Optionally, you can specify simple `tokenizer` option in
each evaluation metric. For the list of tokenizers, see ``cicada --tokenizer-list``.


Test and Reference Data
-----------------------

reference/test set format is as follows:

::

  segment-id |||  sentence (||| .... some information, such as features etc. ...)

which may look like followings:

::

   0 ||| first reference
   0 ||| second reference
   1 ||| first reference for the second input
   1 ||| second reference for the second input



References
----------

.. [1]	 Kishore Papineni, Salim Roukos, Todd Ward, and Wei-Jing
	 Zhu. Bleu: a method for automatic evaluation of machine
	 translation. In Proceedings of 40th Annual Meeting of the
	 Association for Computational Linguistics, pages 311-318,
	 Philadelphia, Pennsylvania, USA, July 2002. Association for
	 Computational Linguistics.

.. [2]	 Sonja Nießen, Franz Josef Och, Gregor Leusch, and Hermann
	 Ney. An evaluation tool for machine translation: Fast
	 evaluation for mt research. In Proceedings of the Second
	 International Conference on Language Resources and Evaluation
	 (LREC 2000), Athens, Greece, 2000.

.. [3]	 C. Tillmann, S. Vogel, H. Ney, A. Zubiaga,
	 and H. Sawaf. Accelerated dp based search for statistical
	 translation. In In European Conf. on Speech Communication and
	 Technology, pages 2667-2670, 1997.

.. [4]	 Matthew Snover, Bonnie Dorr, Richard Schwartz, Linnea
	 Micciulla, and John Makhoul. A study of translation edit rate
	 with targeted human annotation. In In Proceedings of
	 Association for Machine Translation in the Americas, pages
	 223-231, 2006.

.. [5]	 Hideki Isozaki, Tsutomu Hirao, Kevin Duh, Katsuhito Sudoh,
	 and Hajime Tsukada. Automatic evaluation of translation
	 quality for distant language pairs. In Proceedings of the
	 2010 Conference on Empirical Methods in Natural Language
	 Processing, pages 944-952, Cambridge, MA,
	 October 2010. Association for Computational Linguistics.

.. [6]	 Chin-Yew Lin and Franz Josef Och. Automatic evaluation of
	 machine translation quality using longest common subsequence
	 and skip-bigram statistics. In Proceedings of the 42nd
	 Meeting of the Association for Computational Linguistics
	 (ACL'04), Main Volume, pages 605-612, Barcelona, Spain,
	 July 2004.

.. [7]	 Gregor Leusch, Nicola Ueffing, and Hermann Ney. Cder:
	 Efficient mt evaluation using block movements. In In
	 Proceedings of EACL, pages 241-248, 2006.

.. [8]	 Gregor Leusch, Nicola Ueffing, Hermann Ney. A novel
	 string-to-string distance measure with applications to
	 machine translation evaluation. In Proceedings of MT
	 Summit IX, pages 240-247, 2003.

.. [9]	 Philipp Koehn. Statistical significance tests for machine
	 translation evaluation. In Dekang Lin and Dekai Wu, editors,
	 Proceedings of EMNLP 2004, pages 388-395, Barcelona, Spain,
	 July 2004. Association for Computational Linguistics.

.. [10]	 Michael Collins, Philipp Koehn, and Ivona Kučerová. Clause
	 restructuring for statistical machine translation. In ACL
	 '05: Proceedings of the 43rd Annual Meeting on Association
	 for Computational Linguistics, pages 531-540, Morristown, NJ,
	 USA, 2005. Association for Computational Linguistics.
