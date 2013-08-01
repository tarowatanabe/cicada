================
cicada_alignment
================

:Author: Taro Watanabe <taro.watanabe@nict.go.jp>
:Date: 2013-8-1
:Manual section: 1

OPTIONS
-------

  **--source-target** `arg`                    P(target | source) viterbi alignment

  **--target-source** `arg`                    P(source | target) viterbi alignment

  **--span-source** `arg`                      source span data

  **--span-target** `arg`                      target span data

  **--input** `arg`                            input alignment

  **--output** `arg (="-")`                    output alignment

  **--posterior** alignment computation using posteriors

  **--posterior-threshold** `arg (=0.10000000000000001)` 
                                        threshold for posterior

  **--f2e** source target

  **--e2f** target source

  **--itg** itg

  **--max-match** max-match

  **--intersection** intersection

  **--union** union

  **--grow** grow

  **--diag** diag

  **--final** final

  **--final-and** final-and

  **--closure** closure

  **--invert** invert alignment

  **--moses** moses alignment (not GIZA++ alignment)

  **--prob-null** `arg (=0.01)`                NULL probability

  **--prob-union** `arg (=0.5)`                union probability

  **--prob-intersection** `arg (=1)`           intersection probability

  **--threads** `arg`                          # of threads

  **--debug** `[=arg(=1)]`                     debug level

  **--help** help message


