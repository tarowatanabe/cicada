Word alignment in cicada

You can try "cicada-alignment.py" which performs word alignment training in two directions.

Usage: cicada-alignment.py [options]

Options:
  --root-dir=DIRECTORY  root directory for outputs
  --corpus-dir=PREFIX   corpus directory (default: ${root_dir}/corpus)
  --giza-f2e=DIRECTORY  giza directory for P(f|e) (default:
                        ${root_dir}/giza.${f}-${e})
  --giza-e2f=DIRECTORY  giza directory for P(e|f) (default:
                        ${root_dir}/giza.${e}-${f})
  --model-dir=DIRECTORY
                        model directory (default: ${root_dir}/model)
  --alignment-dir=DIRECTORY
                        alignment directory (default: ${model_dir})
  --f=FILE-OR-SUFFIX    source (or 'French')  language file or suffix
  --e=FILE-OR-SUFFIX    target (or 'English') language file or suffix
  --a=FILE-OR-SUFFIX    alignment file or suffix
  --sf=FILE-OR-SUFFIX   source (or 'French')  span file or suffix
  --se=FILE-OR-SUFFIX   target (or 'English') span file or suffix
  --ff=FILE-OR-SUFFIX   source (or 'French')  forest file or suffix
  --fe=FILE-OR-SUFFIX   target (or 'English') forest file or suffix
  --corpus=CORPUS       bilingual trainging corpus prefix
  --alignment=ALIGNMENT
                        alignment methods (default: grow-diag-final-and)
  --first-step=STEP     first step (default: 1)
  --last-step=STEP      last step  (default: 3)
  --iteration-cluster=ITERATION
                        word cluter iterations (default: 50)
  --iteration-model1=ITERATION
                        Model1 iteratins (default: 5)
  --iteration-hmm=ITERATION
                        HMM iteratins    (default: 5)
  --iteration-model4=ITERATION
                        Model4 iteratins    (default: 5)
  --cluster=CLUSTER     # of clusters (default: 50)
  --symmetric           symmetric training
  --posterior           posterior constrained training
  --variational         variational Bayes estimates
  --l0                  L0 regularization
  --p0=P0               parameter for NULL alignment (default: 0.01)
  --insertion-p0=P0     parameter for insertion (default: 0.01)
  --prior-lexicon=PRIOR
                        lexicon model prior (default: 0.01)
  --prior-alignment=PRIOR
                        alignment model prior (default: 0.01)
  --prior-distortion=PRIOR
                        distortion model prior (default: 0.01)
  --prior-fertility=PRIOR
                        fertility model prior (default: 0.01)
  --smooth-lexicon=SMOOTH
                        lower-bound parameter for lexicon model (default:
                        1e-100)
  --smooth-alignment=SMOOTH
                        lower-bound parameter for alignment model (default:
                        1e-100)
  --smooth-distortion=SMOOTH
                        lower-bound parameter for distortion model (default:
                        1e-100)
  --smooth-fertility=SMOOTH
                        lower-bound parameter for fertility model (default:
                        1e-100)
  --l0-alpha=L0_ALPHA   L0 regularization parameter (default: 100)
  --l0-beta=L0_BETA     L0 regularization parameter (default: 0.01)
  --cicada-dir=DIRECTORY
                        cicada directory
  --threads=THREADS     # of thrads for thread-based parallel processing
  --max-malloc=MALLOC   maximum memory in GB (default: 8)
  --pbs                 PBS for launching processes
  --pbs-queue=NAME      PBS queue for launching processes (default: ltg)
  --debug=DEBUG         
  -h, --help            show this help message and exit

Currently, we implemented:

cicada_lexicon_dice
	Dice

cicada_lexicon_model1
	IBM Model1[1] + symmetric learning[2] + posterior constrained learning [3]
	+ naive variational Bayes
	+ L0 regularization [5]
	+ ITG constrained posterior alignment
	+ MaxMatch posterior alignment (using Hungarian algorithm)

cicada_lexicon_hmm
	HMM Model[2] + symmetric learning[2] + posterior constrained learning [3]
	+ naive variational Bayes
	+ L0 regularization [5]
	+ ITG constrained posteriors
	+ MaxMatch posteriors (using Hungarian algorithm)

cicada_lexicon_model4
	IBM Model4[1] + symmetric learning[2] + posterior constrained learning [3]
	+ naive variational Bayes
	+ L0 regularization [5]
	+ ITG constrained posteriors
	+ MaxMatch posteriors (using Hungarian algorithm)

cicada_lexicon_global{,_mpi}
	Trigger based global lexicon[4]

cicada_lexicon
	Learn lexicon model from word aligned data, primarily used for lexical-probabilities for phrase/rule scoring
	+ naive variational Bayes
	+ L0 regularization [5]

cicada_alignment
	Heuristic word alignment combiner (as in Moses grow-diagl-final-and etc.)
	+ ITG constrained combiner
	+ MaxMatch combiner (using Hungarian algorithm)

	Posterior word alignment matrix combiner
  	+ ITG constraint
	+ MaxMatch (using Hungarian algorithm)	

cicada_pyp_itg_learn
    Direct ITG word alignment training (Highly exprerimental)

We also support alignment visulization by cicada_filter_alignment:

cicada_fiter_alignment 
	--source <source language file>
	--target <target langauge file>
	--alignment  <alignment file>
	--alignment2 <secondary alignment file> (optional)
	--inverse (inverse alignment, optional)
	--visualize (required for visualization!)

blue point indicates "intersection"
green point indicates word alignment only in the primary alignment
yello point indicates word alignment only in the secondary alignment


[1]
@article{Brown:1993:MSM:972470.972474,
 author = {Brown, Peter F. and Pietra, Vincent J. Della and Pietra, Stephen A. Della and Mercer, Robert L.},
 title = {The mathematics of statistical machine translation: parameter estimation},
 journal = {Comput. Linguist.},
 issue_date = {June 1993},
 volume = {19},
 issue = {2},
 month = {June},
 year = {1993},
 issn = {0891-2017},
 pages = {263--311},
 numpages = {49},
 url = {http://portal.acm.org/citation.cfm?id=972470.972474},
 acmid = {972474},
 publisher = {MIT Press},
 address = {Cambridge, MA, USA},
} 

[2]
@InProceedings{liang-taskar-klein:2006:HLT-NAACL06-Main,
  author    = {Liang, Percy  and  Taskar, Ben  and  Klein, Dan},
  title     = {Alignment by Agreement},
  booktitle = {Proceedings of the Human Language Technology Conference of the NAACL, Main Conference},
  month     = {June},
  year      = {2006},
  address   = {New York City, USA},
  publisher = {Association for Computational Linguistics},
  pages     = {104--111},
  url       = {http://www.aclweb.org/anthology/N/N06/N06-1014}
}

[3]
@InProceedings{ganchev-gracca-taskar:2008:ACLMain,
  author    = {Ganchev, Kuzman  and  Gra\c{c}a, Jo\~{a}o V.  and  Taskar, Ben},
  title     = {Better Alignments = Better Translations?},
  booktitle = {Proceedings of ACL-08: HLT},
  month     = {June},
  year      = {2008},
  address   = {Columbus, Ohio},
  publisher = {Association for Computational Linguistics},
  pages     = {986--993},
  url       = {http://www.aclweb.org/anthology/P/P08/P08-1112}
}

[4]
@InProceedings{mauser-hasan-ney:2009:EMNLP,
  author    = {Mauser, Arne  and  Hasan, Sa{\v{s}}a  and  Ney, Hermann},
  title     = {Extending Statistical Machine Translation with Discriminative and Trigger-Based Lexicon Models},
  booktitle = {Proceedings of the 2009 Conference on Empirical Methods in Natural Language Processing},
  month     = {August},
  year      = {2009},
  address   = {Singapore},
  publisher = {Association for Computational Linguistics},
  pages     = {210--218},
  url       = {http://www.aclweb.org/anthology/D/D09/D09-1022}
}

[5]
@InProceedings{vaswani-huang-chiang:2012:ACL2012,
  author    = {Vaswani, Ashish  and  Huang, Liang  and  Chiang, David},
  title     = {Smaller Alignment Models for Better Translations: Unsupervised Word Alignment with the l0-norm},
  booktitle = {Proceedings of the 50th Annual Meeting of the Association for Computational Linguistics (Volume 1: Long Papers)},
  month     = {July},
  year      = {2012},
  address   = {Jeju Island, Korea},
  publisher = {Association for Computational Linguistics},
  pages     = {311--319},
  url       = {http://www.aclweb.org/anthology/P12-1033}
}
