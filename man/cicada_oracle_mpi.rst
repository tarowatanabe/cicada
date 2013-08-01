=================
cicada_oracle_mpi
=================

:Author: Taro Watanabe <taro.watanabe@nict.go.jp>
:Date: 2013-8-1
:Manual section: 1

OPTIONS
-------

  **--tstset** `arg`                           test set file(s) (in hypergraph format)

  **--refset** `arg`                           reference set file(s)

  **--output** `arg (="-")`                    output file

  **--forest** output by forest

  **--directory** output in directory

  **--scorer** `arg (=bleu:order=4,exact=true)` 
                                        error metric

  **--scorer-cube** `arg (=200)`               cube pruning size

  **--scorer-beam** `arg (=1.0000000000000001e-05)` 
                                        beam pruning size

  **--max-iteration** `arg`                    # of hill-climbing iteration

  **--min-iteration** `arg`                    # of hill-climbing iteration

  **--apply-exact** exact application

command line options:

  **--debug** `[=arg(=1)]`     debug level

  **--help** help message


