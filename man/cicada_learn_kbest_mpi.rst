======================
cicada_learn_kbest_mpi
======================

:Author: Taro Watanabe <taro.watanabe@nict.go.jp>
:Date: 2013-8-1
:Manual section: 1

OPTIONS
-------

  **--kbest** `arg`                      kbest path

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

  **--learn-mira** online MIRA algorithm

  **--learn-nherd** online NHERD algorithm

  **--learn-arow** online AROW algorithm

  **--learn-cw** online CW algorithm

  **--learn-pegasos** online Pegasos algorithm

  **--learn-cp** cutting plane algorithm

  **--learn-mcp** MERT by cutting plane algorithm

  **--optimize-lbfgs** LBFGS optimizer

  **--optimize-cg** CG optimizer

  **--optimize-sgd** SGD optimizer

  **--regularize-l1** L1-regularization

  **--regularize-l2** L2-regularization

  **--C** `arg (=1)`                     regularization constant

  **--scale** `arg (=1)`                 scaling for weight

  **--eta0** `arg`                       \eta_0 for decay

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

  **--mert-search-local** perform local one-dimensional mert

  **--sample-vector** perform samling

  **--direct-loss** compute loss by directly treating hypothesis 
                                  score

  **--conservative-loss** conservative loss

  **--scale-fixed** fixed scaling

  **--scorer** `arg (=bleu:order=4)`     error metric

  **--scorer-list** list of error metric

  **--unite** unite kbest sharing the same id

  **--debug** `[=arg(=1)]`               debug level

  **--help** help message


