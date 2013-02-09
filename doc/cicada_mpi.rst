========
 cicada
========

---------------------------------------------
a hypergraph toolkit for machine translation
---------------------------------------------

:Author: Taro Watanabe <taro.watanabe@nict.go.jp>
:Date:   2013-2-8
:Manual section: 1

SYNOPSIS
--------

**cicada** [*options*]

DESCRIPTION
-----------



CONFIGURATION OPTIONS
---------------------

Input Options
~~~~~~~~~~~~~

  --input arg                    input file
  --input-id                     id-prefixed input
  --input-bitext                 target sentence prefixed input
  --input-sentence               sentence input
  --input-lattice                lattice input
  --input-forest                 forest input
  --input-span                   span input
  --input-alignment              alignment input
  --input-dependency             dependency input
  --input-directory              input in directory

Grammar Options
~~~~~~~~~~~~~~~

  --goal arg                     goal symbol
  --grammar arg                  grammar specification(s)
  --grammar-list                 list of available grammar specifications
  --tree-grammar arg             tree grammar specification(s)
  --tree-grammar-list            list of available tree grammar specifications

Feature Function Options
~~~~~~~~~~~~~~~~~~~~~~~~

  --feature-function             feature function(s)
  --feature-function-list        list of available feature function(s)
  --output-feature-function arg  output feature function(s)

Operation Options
~~~~~~~~~~~~~~~~~

  --operation arg                operations
  --operation-list               list of available operation(s)

Other List Options
~~~~~~~~~~~~~~~~~~

  --eval-list                    list of available evaluators
  --format-list                  list of available formatters
  --signature-list               list of available signatures
  --stemmer-list                 list of available stemmers
  --tokenizer-list               list of available tokenizers
  --matcher-list                 list of available matchers

COMMANDLINE OPTIONS
-------------------

  --config arg          configuration file
  --debug               debug level
  --help                help message

EXAMPLES
--------



SEE ALSO
--------
