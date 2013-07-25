===========================
 cicada_filter_penntreebank
===========================

-------------------------------------------------
a filter for tranforming penntreebank style trees
-------------------------------------------------

:Author: Taro Watanabe <taro.watanabe@nict.go.jp>
:Date:   2013-7-25
:Manual section: 1

SYNOPSIS
--------

**cicada_filter_penntreebank** [*options*]

DESCRIPTION
-----------



OPTIONS
-------

  **--input** `arg (="-")`     input file

  **--output** `arg (="-")`    output

  **--map** `arg (="")`        map terminal symbols

  **--fix-terminal** fix fragmented terminals

  **--replace-root** `arg`     replace root symbol

  **--unescape** unescape terminal symbols, such as -LRB-, \* etc.

  **--normalize** normalize category, such as [,] [.] etc.

  **--remove-none** remove -NONE-

  **--remove-cycle** remove cycle unary rules

  **--remove-bracket** remove bracket terminal (which may be confued with 
                        non-terminal!)

  **--collapse** collapse unary rules

  **--add-bos-eos** add [BOS]/[EOS] and <s>/</s>

  **--stemmer** `arg`          stemming for terminals

  **--leaf** output leaf nodes

  **--leaf-pos** output leaf/pos nodes

  **--rule** output rules

  **--treebank** output treebank

  **--span** output spans

  **--binarize** span: perform binarization

  **--category** span: added category to span

  **--unary-top** span: use top-most category for unary rules

  **--unary-bottom** span: use bottom-most category for unary rules

  **--unary-root** span: use single category for root

  **--exclude-terminal** span: no terminal in span

  **--skip** skip invalid penntree

  **--validate** validate penntreebank

  **--debug** `[=arg(=1)]`     debug level

  **--help** help message


SEE ALSO
--------


