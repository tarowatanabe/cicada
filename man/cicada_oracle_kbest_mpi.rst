=======================
cicada_oracle_kbest_mpi
=======================

:Author: Taro Watanabe <taro.watanabe@nict.go.jp>
:Date: 2013-8-1
:Manual section: 1

SYNOPSIS
--------

**cicada_oracle_kbest_mpi** [*options*]

OPTIONS
-------

  **--tstset** `arg`                           test set file(s) (in kbest format)

  **--refset** `arg`                           reference set file(s)

  **--output** `arg (="-")`                    output file

  **--segment** initialize by segment score

  **--directory** output in directory

  **--scorer** `arg (=bleu:order=4,exact=true)` 
                                        error metric

  **--max-iteration** `arg`                    # of hill-climbing iteration

  **--min-iteration** `arg`                    # of hill-climbing iteration

command line options:

  **--debug** `[=arg(=1)]`     debug level

  **--help** help message


