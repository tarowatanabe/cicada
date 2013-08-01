==================
cicada_learn_kbest
==================

:Author: Taro Watanabe <taro.watanabe@nict.go.jp>
:Date: 2013-8-1
:Manual section: 1

OPTIONS
-------

  **--kbest** `arg`                      kbest path(s)

  **--input** `arg`                      input path(s) (an alias for --kbest)

  **--oracle** `arg`                     oracle kbest path

  **--refset** `arg`                     reference set file(s)

  **--weights** `arg`                    initial parameter

  **--weights-history** `arg`            parameter history

  **--output** `arg`                     output parameter

  **--output-objective** `arg`           output final objective

  **--iteration** `arg (=100)`           max # of iterations

  **--learn-softmax** Softmax objective

  **--learn-xbleu** xBLEU objective

  **--learn-linear** liblinear algorithm

  **--learn-svm** structural SVM

  **--optimize-lbfgs** LBFGS optimizer

  **--optimize-cg** CG optimizer

  **--solver** `arg`                     liblinear solver type (default: 1)
                                   0: L2-regularized logistic regression 
                                      (primal)
                                   1: L2-regularized L2-loss support vector 
                                      classification (dual)
                                   2: L2-regularized L2-loss support vector 
                                      classification (primal)
                                   3: L2-regularized L1-loss support vector 
                                      classification (dual)
                                   5: L1-regularized L2-loss support vector 
                                      classification
                                   6: L1-regularized logistic regression
                                   7: L2-regularized logistic regression (dual)
                                  11: L2-regularized L2-loss epsilon support 
                                      vector regression (primal)
                                  12: L2-regularized L2-loss epsilon support 
                                      vector regression (dual)
                                  13: L2-regularized L1-loss epsilon support 
                                      vector regression (dual)
                                  

  **--regularize-l1** L1-regularization

  **--regularize-l2** L2-regularization

  **--C** `arg (=1)`                     regularization constant

  **--scale** `arg (=1)`                 scaling for weight

  **--eta0** `arg`                       \eta_0 for decay

  **--eps** `arg`                        tolerance for liblinear

  **--order** `arg (=4)`                 ngram order for xBLEU

  **--annealing** annealing

  **--quenching** quenching

  **--temperature** `arg (=0)`           temperature

  **--temperature-start** `arg (=1000)`  start temperature for annealing

  **--temperature-end** `arg (=0.001)`   end temperature for annealing

  **--temperature-rate** `arg (=0.5)`    annealing rate

  **--quench-start** `arg (=0.01)`       start quench for annealing

  **--quench-end** `arg (=100)`          end quench for annealing

  **--quench-rate** `arg (=10)`          quenching rate

  **--loss-margin** direct loss margin

  **--softmax-margin** softmax margin

  **--line-search** perform line search in each iteration

  **--mert-search** perform one-dimensional mert

  **--sample-vector** perform samling

  **--direct-loss** compute loss by directly treating hypothesis 
                                  score

  **--conservative-loss** conservative loss

  **--scale-fixed** fixed scaling

  **--scorer** `arg (=bleu:order=4)`     error metric

  **--scorer-list** list of error metric

  **--unite** unite kbest sharing the same id

  **--threads** `arg`                    # of threads

  **--debug** `[=arg(=1)]`               debug level

  **--help** help message


