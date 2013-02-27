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

  **--input** `arg (="-")`             input file

  **--input-id** id-prefixed input

  **--input-bitext** target sentence prefixed input

  **--input-sentence** sentence input

  **--input-lattice** lattice input

  **--input-forest** forest input

  **--input-span** span input

  **--input-alignment** alignment input

  **--input-dependency** dependency input

  **--input-directory** input in directory

Grammar Options
~~~~~~~~~~~~~~~

  **--goal** `arg (=[s])`              goal symbol

  **--grammar** `arg`                  grammar specification(s)

  **--grammar-list** list of available grammar specifications

  **--tree-grammar** `arg`             tree grammar specification(s)

  **--tree-grammar-list** list of available grammar specifications

Feature Function Options
~~~~~~~~~~~~~~~~~~~~~~~~

  **--feature-function** `arg`         feature function(s)

  **--feature-function-list** list of available feature function(s)

  **--output-feature-function** `arg`  output feature function(s)

Operation Options
~~~~~~~~~~~~~~~~~

  **--operation** `arg`                operations

  **--operation-list** list of available operation(s)


Other List Options
~~~~~~~~~~~~~~~~~~

  **--eval-list** list of available evaluators

  **--format-list** list of available formatters

  **--signature-list** list of available signatures

  **--stemmer-list** list of available stemmers

  **--tokenizer-list** list of available tokenizers

  **--matcher-list** list of available matchers


COMMANDLINE OPTIONS
-------------------

  **--config** `arg`           configuration file

  **--threads** `arg`          # of threads

  **--debug** `[=arg(=1)]`     debug level

  **--help** help message

EXAMPLES
--------


INPUTS
------

Input to cicada is controlled by the **--input-\*** flags, and each
line may consists of various components split by a delimiter "|||":

::

[id |||] [latttice or sentence] [||| forest] [||| span] [||| alignment] [||| dependency] [||| sentence]*

The id is optional, and if an id is not supplied here, cicada assumes
a numerical order starting from 0.
By default, cicada assumes sentence input, but if **--input-lattice**
option is specified, lattice is assumed, instead.
It is followed by optional entries, `forest`, `span`, `alignment`,
`dependency` and a set of `sentence`.


SEE ALSO
--------
