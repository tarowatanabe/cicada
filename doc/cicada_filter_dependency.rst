=========================
 cicada_filter_dependency
=========================

---------------------------------
a filter for denendency structure
---------------------------------

:Author: Taro Watanabe <taro.watanabe@nict.go.jp>
:Date:   2013-2-28
:Manual section: 1

SYNOPSIS
--------

**cicada_filter_dependency** [*options*]

DESCRIPTION
-----------



OPTIONS
-------

  **--input** `arg (="-")`         input file

  **--output** `arg (="-")`        output file

  **--map** `arg`                  map file

  **--goal** `arg (=[s])`          goal symbol

  **--non-terminal** `arg (=[x])`  non-terminal symbol

  **--mst** tranform MST dependency

  **--conll** tranform CoNLL dependency

  **--malt** tranform Malt dependency

  **--cabocha** tranform Cabocha dependency

  **--khayashi** tranform KHayashi dependency

  **--khayashi-forest** tranform KHayashi Forest dependency

  **--cicada** tranform CICADA dependency

  **--projective** project into projective dependency

  **--relation** assing relation to POS

  **--unescape** unescape terminal symbols, such as -LRB-, \* etc.

  **--normalize** normalize category, such as [,] [.] etc.

  **--func** use function word as a head in bunsetsu

  **--forest** output as a forest

  **--head** output hypergraph with explicit head

  **--span** output as a set of spans

  **--binarize** output binarized spans

  **--category** output spans with category labels

  **--leaf** output leafs (temrinals)

  **--leaf-pos** output leaf-pos (temrinal+POS)

  **--help** help message


SEE ALSO
--------


