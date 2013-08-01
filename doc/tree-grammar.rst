Tree-Grammar
============

We can use tree-transduce grammar as in:
::

  [X]([X,1] [X](a b c)) ||| [Y](A B [Y,1]) ||| feature1=0.5 feature2=-40 ||| attr="good rule"

``\``, ``(`` and ``)`` should be escaped by ``\\``, ``\(`` and ``\)``.

In cicada file-based grammar can be loaded in additin to pre-defined
fallback-rules by the ``--tree-grammar`` options:

file-name: indexed tree grammar or plain text tree grammar
	max-span=[int] maximum span (<=0 for no-constraint)
	cky|cyk=[true|false] indexing for CKY|CYK parsing/composition
	key-value=[true|false] store key-value format of features/attributes
	populate=[true|false] "populate" by pre-fetching
	feature-prefix=[prefix for feature name] add prefix to the default feature name: tree-rule-table
	attribute-prefix=[prefix for attribute name] add prefix to the default attribute name: tree-rule-table
	feature0=[feature-name]
	feature1=[feature-name]
	...
	attribute0=[attribute-name]
	attribute1=[attribute-name]
	...

fallback: fallback source-to-target, tree-to-{string,tree} transfer rule
	goal=[default goal label] target side goal
	non-terminal=[defaut non-terminal] target side non-terminal


For defailt, see ``cicada --tree-grammar-list``
