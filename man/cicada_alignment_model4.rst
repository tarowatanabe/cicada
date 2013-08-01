=======================
cicada_alignment_model4
=======================

:Author: Taro Watanabe <taro.watanabe@nict.go.jp>
:Date: 2013-8-1
:Manual section: 1

SYNOPSIS
--------

**cicada_alignment_model4** [*options*]

OPTIONS
-------

  **--source** `arg`                           source file

  **--target** `arg`                           target file

  **--alignment** `arg`                        alignment file

  **--dependency-source** `arg`                source dependency file

  **--dependency-target** `arg`                target dependency file

  **--span-source** `arg`                      source span file

  **--span-target** `arg`                      target span file

  **--classes-source** `arg`                   source classes file

  **--classes-target** `arg`                   target classes file

  **--lexicon-source-target** `arg`            lexicon model for P(target | source)

  **--lexicon-target-source** `arg`            lexicon model for P(source | target)

  **--alignment-source-target** `arg`          alignment model for P(target | source)

  **--alignment-target-source** `arg`          alignment model for P(source | target)

  **--distortion-source-target** `arg`         distortion model for P(target | source)

  **--distortion-target-source** `arg`         distortion model for P(source | target)

  **--fertility-source-target** `arg`          fertility model for P(target | source)

  **--fertility-target-source** `arg`          fertility model for P(source | target)

  **--insertion-source-target** `arg (=0.01)`  insertion model for P(target | source)

  **--insertion-target-source** `arg (=0.01)`  insertion model for P(source | target)

  **--length-source-target** `arg (=1)`        length model for P(target | source)

  **--length-target-source** `arg (=1)`        length model for P(source | target)

  **--output-lexicon-source-target** `arg`     lexicon model output for P(target | 
                                        source)

  **--output-lexicon-target-source** `arg`     lexicon model output for P(source | 
                                        target)

  **--output-alignment-source-target** `arg`   alignment model output for P(target | 
                                        source)

  **--output-alignment-target-source** `arg`   alignment model output for P(source | 
                                        target)

  **--output-distortion-source-target** `arg`  distortion model output for P(target | 
                                        source)

  **--output-distortion-target-source** `arg`  distortion model output for P(source | 
                                        target)

  **--output-fertility-source-target** `arg`   fertility model output for P(target | 
                                        source)

  **--output-fertility-target-source** `arg`   fertility model output for P(source | 
                                        target)

  **--viterbi-source-target** `arg`            viterbi for P(target | source)

  **--viterbi-target-source** `arg`            viterbi for P(source | target)

  **--projected-source** `arg`                 source dependnecy projected from target

  **--projected-target** `arg`                 target dependency projected from source

  **--posterior-source-target** `arg`          posterior for P(target | source)

  **--posterior-target-source** `arg`          posterior for P(source | target)

  **--posterior-combined** `arg`               posterior for P(source | target) 
                                        P(target | source)

  **--iteration-model1** `arg (=5)`            max Model1 iteration

  **--iteration-hmm** `arg (=5)`               max HMM iteration

  **--iteration-model4** `arg (=5)`            max Model4 iteration

  **--symmetric** symmetric training

  **--posterior** posterior constrained training

  **--variational-bayes** variational Bayes estimates

  **--pgd** projected gradient descent

  **--dynamic** dynamically re-compute alignment

  **--itg** ITG alignment

  **--max-match** maximum matching alignment

  **--moses** Moses alignment foramt

  **--permutation** permutation

  **--hybrid** hybrid projective dependency parsing

  **--degree2** degree2 non-projective dependency 
                                        parsing

  **--mst** MST non-projective dependency parsing

  **--single-root** single root dependency

  **--p0** `arg (=0.01)`                       parameter for NULL alignment

  **--prior-lexicon** `arg (=0.01)`            Dirichlet prior for variational Bayes

  **--smooth-lexicon** `arg (=1e-100)`         smoothing parameter for uniform 
                                        distribution

  **--prior-alignment** `arg (=0.01)`          Dirichlet prior for variational Bayes

  **--smooth-alignment** `arg (=1e-100)`       smoothing parameter for uniform 
                                        distribution

  **--prior-distortion** `arg (=0.01)`         Dirichlet prior for variational Bayes

  **--smooth-distortion** `arg (=1e-100)`      smoothing parameter for uniform 
                                        distribution

  **--prior-fertility** `arg (=0.01)`          Dirichlet prior for variational Bayes

  **--smooth-fertility** `arg (=1e-100)`       smoothing parameter for uniform 
                                        distortion

  **--l0-alpha** `arg (=100)`                  L0 regularization

  **--l0-beta** `arg (=0.01)`                  L0 regularization

  **--threshold** `arg (=0)`                   write with beam-threshold (<= 0.0 
                                        implies no beam)

  **--threads** `arg`                          # of threads

  **--debug** `[=arg(=1)]`                     debug level

  **--help** help message


