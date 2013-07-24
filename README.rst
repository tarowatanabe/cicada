========
 cicada
========

A statistical machine translation toolkit based on a semiring parsing
framework [1]_. Based on the generic framework, we can

   - learn model(s): tree-to-string, string-to-tree, string-to-string (with or without latent tree),
     word alignment, grammar for parsing
   - translate sentence, lattice and/or parsed tree
   - align between {hypergraph,lattice}-to-lattice (currently, we assume penn-treebank style hypergraph
     and sentence-like lattice)
   - align two sentences using symmetized IBM Model1/HMM/IBM Model4.
   - (dependency) parse lattices (or sentences)
   - analyze forest/tree/lattice

The cicada toolkit is developed at Multilingual Translation
Laboratory, Universal Communication Institute, National Institute of
Information and Communications Technology (NICT).

Remark: cicada is 蝉(CJK UNIFIED IDEOGRAPH-8749), (or セミ) in Japanese, pronounced "SEMI"

Quick Start
-----------

Compile
```````

Get the source code from `cicada <http://www2.nict.go.jp/univ-com/multi_trans/cicada>`_,
or from `github.com <http://github.com/tarowatanabe/cicada>`_, and
simply follow the GNU standard pipiline. For details, see `BUILD.rst`.

.. code:: bash

   ./autogen.sh (required when you get the code by git clone)
   ./configure
   make
   make install (optional)

Run
```

You can find a sample grammar file at *samples* directory together with
*ngram* language model. Here is an example run (Note that \\ indicates
shell's newline).

.. code:: bash

   ./progs/cicada \
      --input samples/scfg/input.txt \
      --grammar samples/scfg/grammar.bin \
      --grammar "glue:straight=true,inverted=false,non-terminal=[x],goal=[s]" \
      --grammar "insertion:non-terminal=[x]" \
      --feature-function "ngram:file=samples/scfg/ngram.bin" \
      --feature-function word-penalty \
      --feature-function rule-penalty \
      --operation compose-cky \
      --operation apply:prune=true,size=100,weights=samples/scfg/weights \
      --operation output:file=-,kbest=10,weights=samples/scfg/weights \
      --debug

This sample means:

  - Input is `samples/scfg/input.txt`
  - Three grammars:

    - The SCFG file is `samples/scfg/grammar.bin` which is a
      binary version of `samples/scfg/grammar.bz2`.
    - Additional grammar is a "glue grammar" which consists of two rules
      of "[s] -> <[x], [x]>" and "[s] -> <[s,1] [x,1], [s,1] [x,1]>"
    - Another additional grammar is an "insertion grammar" which simply
      copies the input string to output string, "[x] -> <word-x, word-x>"

  - Three feature functions:

    - 5-gram language model from `samples/scfg/ngram.bin` which is a
      binary version of `samples/scfg/ngram.bz2`.
    - word penalty feature which penalize by the number of words in
      the target side.
    - rule penalty feature which penaltize by the number of words in a
      derivation.
    - In addition, there exist features already defined for each
      hierarchical phrase pair. For example, see `samples/scfg/grammar.bz2`.

  - Actual operation:

    1. Input is composed by CKY algorithm (compose-cky) which result
       in a hypergraph.
    2. Cube-pruning (apply) to apply feature functions using 100 as a
       histogram pruning threshold using the weights at
       `samples/scfg/weights`.
    3. 10-best derivations are computed and output at
       `-` (stdout) using `samples/scfg/weights` as a
       weight vector to compute the score for each derivation.

In depth
````````

In order to train a model, see `doc/training.rst` which describes how
to create your own {tree,string}-to-{tree,string} models, and tune
parameters.


Descriptions
------------

Basically, we have four distinct structures:

   - lattice: a representation of graph implemented as a
     two-dimentional array (see `doc/lattice.rst`).
   - grammar: a collection of WFST implemented as a trie structure
     (see `doc/grammar.rst`).
   - tree-grammar: a collectin of WFSTT (tree-transducer) implemented
     as a (nested) trie structure (see `doc/tree-grammar.rst`).
   - hypergraph: a compact representation of set of trees (or forest)
     (see `doc/hypergraph.rst`).

Translation/parsing can be carried out by:

   - A lattice (or sentence) is composed with a grammar, generating a
     hypergraph [2]_ [24]_.
   - A lattice (or sentence) is composed with a tree-grammar,
     generating a hypergraph [27]_.
   - A lattice (or sentence) is composed with a phrasal grammar,
     generating a phrasal hypergraph [4]_.
   - A hypergraph/forest (or parse-tree) is composed with a phrasal
     grammar, generating another hypergraph [3]_.
   - A hypergraph/forest (or parse-tree) is composed with a tree
     grammar, generating another hypergraph [4]_.

Alignment can be carried out by:

   - A lattice is composed with dictionary, generating alignment
     hypergraph, or
   - A hypergraph is composed with dictinary, generating alignment
     hypergraph [20]_.
   - In order to support word alignment training, we can learn
     Model1/HMM/Model4 by symmetized learning [22]_ or
     symmetric posterior constrained learning [23]_ with smoothing via
     variational Bayes or via L0 prior.

     Final combined alignment can be generated either by heuristic
     (AKA grow-diag-final-and etc.) or by ITG or max-matching from
     posterior probabilities.
     Also, lexicon model can be discriminatively trained [28]_.
     For details of the training process, please refer to
     `doc/training.rst` and `doc/alignment.rst`.

Dependency parsing can be carried out by:

   - A lattice is dependency parsed by arc-standard, arc-eager, hybrid, degree2,
     which generates derivation hypergraph.
   - Forests are rescored by dependency features (TODO).
     We support dependency projection [32]_ with Model1/HMM posterior
     probabilies so that we can train arbitrary dependency parses
     after projections.

After the hypergraph generation, you can:

   - Additional features are evaluated to generate another hypergraph [4]_.
     cicada implementes cube-pruning [4]_, cube-growing [4]_,
     incremental [18]_ and exact (and stateless-inside-algorithm)
     methods.

     * cube-growing employs coarse-heuristics [11]_, such as lower-order
       ngrams etc.
     * cube-pruning implements algorithm 2 of faster cube pruning [31]_.

   - Perform variational decoding for hypergraph [10]_ or MBR decoding for hypergraph [12]_
     based on the expected ngram-counts over forest [13]_.
   - K-best sentences are generated from hypergraph [5]_.
   - Generate oracle translations (BLEU only).

Or, you can combine outputs from multiple systems by [29]_:

   - Perform parsing over nbests (use your favorite parser, such as
     Berkeley parser/Stanford parser etc.)
   - Generate context-free confusion forest by combining trees (not confusion network!)
     It is performed by collecting rules from parse trees, and
     generate by Earley algorithm
   - Generate k-best translations after feature application etc.

Or, a conventional system combination strategy of [14]_:

   - Create lattice from n-best list by incremental merging
   - Construct hypergraph by linear grammar (grammar-glue-straight + grammar-insertion)
   - Generate k-best translations after feature application etc.

Monolingual grammar learning is implemented:

   - A simple PCFG by simply extracting rules.
   - Learn latent annotated PCFG by split/merge process with an EM
     algorihtm [25]_.
   - Also, learn coarse grammars from the latent annotated PCFG for
     coarse-to-fine parsing [26]_.

Phrase/synchronou-rule/tree-to-string/string-to-tree extraction/scoring are implemented:

   - A conventional phrase extract algorithm in Moses.
   - A conventional hierarchical phrase extraction algorithm in Hiero
     with or without syntax augmentation [15]_.
   - Tree-to-string/strint-to-tree extractin from forest [16]_ [27]_.
   - Tree-to-tree rule extraction from forest [17]_ (experimental).
   - max-scope constraints to limit the grammar size [34]_.
   - After count extraction, you can perform map/reduce to compute
     model scores [19]_.
   - Then, prune your model based on Fisher's exact test [38]_.

Various learning components are implemented:

   - Large feature set from input lattice/hypergraph on large training
     data via MaxEnt (optimized by LBFGS) [3]_
   - Large/small featuer set from kbests on large/small traning data
     via MaxEnt (LBFGS)/liblinear [30]_
   - Large feature set on small devset with MIRA [6]_ [7]_, but with
     hypergraph
   - Small feature set on small devset learned by hypergraph-MERT [8]_
   - Small/large feature set on small devset learned by
     hypergraph-MaxEnt (optimized by LBFGS or SGD) + softmax-margin [9]_
   - Small/large feature set learned by iteratively construncting
     training samples with rank-learning.

     * optimization by LBFGS/liblinear etc. (similar to [33]_, but differ in kbest handling).
     * larger batching with optimized updates [37]_.
     * We have a script-based implementation + single-binary implementation for efficiency

   - xBLEU objective learned either by L-BFGS or SGD, which directly
     maximize expected-BLEU (not BLEU expectaiton) [35]_.
     Now, this is a recommended optimization method (either kbest or hypergraph learning)
   - We support feature selection by kbest-feature merging [36]_
   - Asynchronous online learning employed in [6]_.

Feature functions:

   -  The ngram language model feaature supports expgram [39]_ and
      kenlm [40]_.

Word clustering tool is also included to support word alignment
learning + translation [20]_.

References
----------

.. [1]	 Zhifei Li and Jason Eisner. First- and second-order
	 expectation semirings with applications to minimum-risk
	 training on translation forests. In Proceedings of the 2009
	 Conference on Empirical Methods in Natural Language
	 Processing, pages 40-51, Singapore, August 2009. Association
	 for Computational Linguistics.

.. [2]	 Christopher Dyer, Smaranda Muresan, and Philip
	 Resnik. Generalizing word lattice translation. In Proceedings
	 of ACL-08: HLT, pages 1012-1020, Columbus, Ohio,
	 June 2008. Association for Computational Linguistics.

.. [3]	 Chris Dyer and Philip Resnik. Context-free reordering,
	 finite-state translation. In Human Language Technologies: The
	 2010 Annual Conference of the North American Chapter of the
	 Association for Computational Linguistics, pages 858-866, Los
	 Angeles, California, June 2010. Association for Computational
	 Linguistics.

.. [4]	 Liang Huang and David Chiang. Forest rescoring: Faster
	 decoding with integrated language models. In Proceedings of
	 the 45th Annual Meeting of the Association of Computational
	 Linguistics, pages 144-151, Prague, Czech Republic,
	 June 2007. Association for Computational Linguistics.

.. [5]	 Liang Huang and David Chiang. Better k-best parsing. In
	 Proceedings of the Ninth International Workshop on Parsing
	 Technology, pages 53-64, Vancouver, British Columbia,
	 October 2005. Association for Computational Linguistics.

.. [6]	 David Chiang, Kevin Knight, and Wei Wang. 11,001 new features
	 for statistical machine translation. In Proceedings of Human
	 Language Technologies: The 2009 Annual Conference of the
	 North American Chapter of the Association for Computational
	 Linguistics, pages 218-226, Boulder, Colorado,
	 June 2009. Association for Computational Linguistics.

.. [7]	 Taro Watanabe, Jun Suzuki, Hajime Tsukada, and Hideki
	 Isozaki. Online large-margin training for statistical machine
	 translation. In Proceedings of the 2007 Joint Conference on
	 Empirical Methods in Natural Language Processing and
	 Computational Natural Language Learning (EMNLP-CoNLL), pages
	 764-773, Prague, Czech Republic, June 2007. Association for
	 Computational Linguistics.

.. [8]	 Shankar Kumar, Wolfgang Macherey, Chris Dyer, and Franz
	 Och. Efficient minimum error rate training and minimum
	 bayes-risk decoding for translation hypergraphs and
	 lattices. In Proceedings of the Joint Conference of the 47th
	 Annual Meeting of the ACL and the 4th International Joint
	 Conference on Natural Language Processing of the AFNLP, pages
	 163-171, Suntec, Singapore, August 2009. Association for
	 Computational Linguistics.

.. [9]	 Kevin Gimpel and Noah A. Smith. Softmax-margin crfs: Training
	 log-linear models with cost functions. In Human Language
	 Technologies: The 2010 Annual Conference of the North
	 American Chapter of the Association for Computational
	 Linguistics, pages 733-736, Los Angeles, California,
	 June 2010. Association for Computational Linguistics.

.. [10]	 Zhifei Li, Jason Eisner, and Sanjeev Khudanpur. Variational
	 decoding for statistical machine translation. In Proceedings
	 of the Joint Conference of the 47th Annual Meeting of the ACL
	 and the 4th International Joint Conference on Natural
	 Language Processing of the AFNLP, pages 593-601, Suntec,
	 Singapore, August 2009. Association for Computational
	 Linguistics.

.. [11]	 David Vilar and Hermann Ney. On lm heuristics for the cube
	 growing algorithm. In Annual Conference of the European
	 Association for Machine Translation, pages 242-249,
	 Barcelona, Spain, May 2009.

.. [12]	 John DeNero, David Chiang, and Kevin Knight. Fast consensus
	 decoding over translation forests. In Proceedings of the
	 Joint Conference of the 47th Annual Meeting of the ACL and
	 the 4th International Joint Conference on Natural Language
	 Processing of the AFNLP, pages 567-575, Suntec, Singapore,
	 August 2009. Association for Computational Linguistics.

.. [13]	 John DeNero, Shankar Kumar, Ciprian Chelba, and Franz
	 Och. Model combination for machine translation. In Human
	 Language Technologies: The 2010 Annual Conference of the
	 North American Chapter of the Association for Computational
	 Linguistics, pages 975-983, Los Angeles, California,
	 June 2010. Association for Computational Linguistics.

.. [14]	 Antti-Veikko Rosti, Bing Zhang, Spyros Matsoukas, and Richard
	 Schwartz. Incremental hypothesis alignment with flexible
	 matching for building confusion networks: BBN system
	 description for WMT09 system combination task. In Proceedings
	 of the Fourth Workshop on Statistical Machine Translation,
	 pages 61-65, Athens, Greece, March 2009. Association for
	 Computational Linguistics.

.. [15]	 Andreas Zollmann and Stephan Vogel. New parameterizations and
	 features for pscfg-based machine translation. In Proceedings
	 of the 4th Workshop on Syntax and Structure in Statistical
	 Translation, pages 110-117, Beijing, China,
	 August 2010. Coling 2010 Organizing Committee.

.. [16]	 Haitao Mi and Liang Huang. Forest-based translation rule
	 extraction. In Proceedings of the 2008 Conference on
	 Empirical Methods in Natural Language Processing, pages
	 206-214, Honolulu, Hawaii, October 2008. Association for
	 Computational Linguistics.

.. [17]	 Yang Liu, Yajuan Lü, and Qun Liu. Improving tree-to-tree
	 translation with packed forests. In Proceedings of the Joint
	 Conference of the 47th Annual Meeting of the ACL and the 4th
	 International Joint Conference on Natural Language Processing
	 of the AFNLP, pages 558-566, Suntec, Singapore,
	 August 2009. Association for Computational Linguistics.

.. [18]	 Liang Huang and Haitao Mi. Efficient incremental decoding for
	 tree-to-string translation. In Proceedings of the 2010
	 Conference on Empirical Methods in Natural Language
	 Processing, pages 273-283, Cambridge, MA,
	 October 2010. Association for Computational Linguistics.

.. [19]	 Chris Dyer, Aaron Cordova, Alex Mont, and Jimmy Lin. Fast,
	 easy, and cheap: Construction of statistical machine
	 translation models with MapReduce. In Proceedings of the
	 Third Workshop on Statistical Machine Translation, pages
	 199-207, Columbus, Ohio, June 2008. Association for
	 Computational Linguistics.

.. [20]	 Jason Riesa and Daniel Marcu. Hierarchical search for word
	 alignment. In Proceedings of the 48th Annual Meeting of the
	 Association for Computational Linguistics, pages 157-166,
	 Uppsala, Sweden, July 2010. Association for Computational
	 Linguistics.

.. [21]	 Jakob Uszkoreit and Thorsten Brants. Distributed word
	 clustering for large scale class-based language modeling in
	 machine translation. In Proceedings of ACL-08: HLT, pages
	 755-762, Columbus, Ohio, June 2008. Association for
	 Computational Linguistics.

.. [22]	 Percy Liang, Ben Taskar, and Dan Klein. Alignment by
	 agreement. In Proceedings of the Human Language Technology
	 Conference of the NAACL, Main Conference, pages 104-111, New
	 York City, USA, June 2006. Association for Computational
	 Linguistics.

.. [23]	 Kuzman Ganchev, João V. Graça, and Ben Taskar. Better
	 alignments = better translations? In Proceedings of ACL-08:
	 HLT, pages 986-993, Columbus, Ohio, June 2008. Association
	 for Computational Linguistics.

.. [24]	 Dan Klein and Christopher D. Manning. Parsing and
	 hypergraphs. In IN IWPT, pages 123-134, 2001.

.. [25]	 Slav Petrov, Leon Barrett, Romain Thibaux, and Dan
	 Klein. Learning accurate, compact, and interpretable tree
	 annotation. In Proceedings of the 21st International
	 Conference on Computational Linguistics and 44th Annual
	 Meeting of the Association for Computational Linguistics,
	 pages 433-440, Sydney, Australia, July 2006. Association for
	 Computational Linguistics.

.. [26]	 Slav Petrov and Dan Klein. Improved inference for
	 unlexicalized parsing. In Human Language Technologies 2007:
	 The Conference of the North American Chapter of the
	 Association for Computational Linguistics; Proceedings of the
	 Main Conference, pages 404-411, Rochester, New York,
	 April 2007. Association for Computational Linguistics.

.. [27]	 Michel Galley, Mark Hopkins, Kevin Knight, and Daniel
	 Marcu. What's in a translation rule? In Daniel Marcu Susan
	 Dumais and Salim Roukos, editors, HLT-NAACL 2004: Main
	 Proceedings, pages 273-280, Boston, Massachusetts, USA, May
	 2 - May 7 2004. Association for Computational Linguistics.

.. [28]	 Arne Mauser, Saša Hasan, and Hermann Ney. Extending
	 statistical machine translation with discriminative and
	 trigger-based lexicon models. In Proceedings of the 2009
	 Conference on Empirical Methods in Natural Language
	 Processing, pages 210-218, Singapore,
	 August 2009. Association for Computational Linguistics.

.. [29]	 Taro Watanabe and Eiichiro Sumita. Machine translation system
	 combination by confusion forest. In Proceedings of the 49th
	 Annual Meeting of the Association for Computational
	 Linguistics: Human Language Technologies, pages 1249-1257,
	 Portland, Oregon, USA, June 2011. Association for
	 Computational Linguistics.

.. [30]	 Rong-En Fan, Kai-Wei Chang, Cho-Jui Hsieh, Xiang-Rui Wang,
	 and Chih-Jen Lin. LIBLINEAR: A library for large linear
	 classification. Journal of Machine Learning Research,
	 9:1871-1874, 2008.

.. [31]	 Andrea Gesmundo and James Henderson. Faster Cube Pruning. In
	 Marcello Federico, Ian Lane, Michael Paul, and François Yvon,
	 editors, Proceedings of the seventh International Workshop on
	 Spoken Language Translation (IWSLT), pages 267-274, 2010.

.. [32]	 Wenbin Jiang and Qun Liu. Dependency parsing and projection
	 based on word-pair classification. In Proceedings of the 48th
	 Annual Meeting of the Association for Computational
	 Linguistics, pages 12-20, Uppsala, Sweden,
	 July 2010. Association for Computational Linguistics.

.. [33]	 Mark Hopkins and Jonathan May. Tuning as ranking. In
	 Proceedings of the 2011 Conference on Empirical Methods in
	 Natural Language Processing, pages 1352-1362, Edinburgh,
	 Scotland, UK., July 2011. Association for Computational
	 Linguistics.

.. [34]	 Mark Hopkins and Greg Langmead. SCFG decoding without
	 binarization. In Proceedings of the 2010 Conference on
	 Empirical Methods in Natural Language Processing, pages
	 646-655, Cambridge, MA, October 2010. Association for
	 Computational Linguistics.

.. [35]	 Antti-Veikko Rosti, Bing Zhang, Spyros Matsoukas, and Richard
	 Schwartz. Expected bleu training for graphs: Bbn system
	 description for wmt11 system combination task. In Proceedings
	 of the Sixth Workshop on Statistical Machine Translation,
	 pages 159-165, Edinburgh, Scotland, July 2011. Association
	 for Computational Linguistics.

.. [36]	 Patrick Simianer, Stefan Riezler, and Chris Dyer. Joint
	 feature selection in distributed stochastic learning for
	 large-scale discriminative training in smt. In Proceedings of
	 the 50th Annual Meeting of the Association for Computational
	 Linguistics (Volume 1: Long Papers), pages 11-21, Jeju
	 Island, Korea, July 2012. Association for Computational
	 Linguistics.

.. [37]	 Taro Watanabe. Optimized online rank learning for machine
	 translation. In Proceedings of the 2012 Conference of the
	 North American Chapter of the Association for Computational
	 Linguistics: Human Language Technologies, pages 253-262,
	 Montréal, Canada, June 2012. Association for Computational
	 Linguistics.

.. [38]	 Howard Johnson, Joel Martin, George Foster, and Roland
	 Kuhn. Improving translation quality by discarding most of the
	 phrasetable. In Proceedings of the 2007 Joint Conference on
	 Empirical Methods in Natural Language Processing and
	 Computational Natural Language Learning (EMNLP-CoNLL), pages
	 967-975, Prague, Czech Republic, June 2007. Association for
	 Computational Linguistics.

.. [39]	 Taro Watanabe, Hajime Tsukada, and Hideki Isozaki. A succinct
	 n-gram language model. In Proceedings of the ACL-IJCNLP 2009
	 Conference Short Papers, pages 341-344, Suntec, Singapore,
	 August 2009. Association for Computational Linguistics.

.. [40]	 Kenneth Heafield. Kenlm: Faster and smaller language model
	 queries. In Proceedings of the Sixth Workshop on Statistical
	 Machine Translation, pages 187-197, Edinburgh, Scotland,
	 July 2011. Association for Computational Linguistics.

