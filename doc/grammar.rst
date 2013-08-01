Grammar
=======

An example of synchronous CFG rules look like followings:
::

  [lhs] ||| [x,1] abc bb ||| ddd [x,1] ||| feature=1.0 feature2=5.6 ||| attr1=1.6
  [lhs] ||| [x,1] abc [x,2] bb ||| ddd [x,1] fff [x,2] ||| feature=-1 feature2=8 ||| attr2="bad"

where ``|||`` is used as a delimiter for the shared left hand side,
source and target rules. In addition to features, you can specify
additional meta data, called attributes which can take either
integer/float or JSON string. When indexing synchronous grammars as a
binary data, it is recommended to use a simple feature attributes
without "keys":
::

  [lhs] ||| [x,1] abc bb ||| ddd [x,1] ||| 1.0 5.6 ||| 1.6
  [lhs] ||| [x,1] abc [x,2] bb ||| ddd [x,1] fff [x,2] ||| -1 8 ||| 0.0


In cicada, various file-based gramamr can be loaded in addition to
pre-defined glue-rules using the ``--grammar`` options:

file-name: indexed grammar or plain text grammar
	max-span=[int] maximum span (<=0 for no-constraint)
	key-value=[true|false] store key-value format of features/attributes
	populate=[true|false] "populate" by pre-fetching
	feature-prefix=[prefix for feature name] add prefix to the default feature name: rule-table
	attribute-prefix=[prefix for attribute name] add prefix to the default attribute name: rule-table
	feature0=[feature-name]
	feature1=[feature-name]
	...
	attribute0=[attribute-name]
	attribute1=[attribute-name]
	...

glue: glue rules
	goal=[goal non-terminal]
	non-terminal=[default non-terminal]
	fallback=[file] the list of fallback non-terminal
	straight=[true|false] straight glue-rule
	inverted=[true|false] inverted glue-rule

insertion: terminal insertion rule
	non-terminal=[defaut non-terminal]
	fallback=[file] the list of fallback non-terminal

deletion: terminal deletion rule
	non-terminal=[defaut non-terminal]
	fallback=[file] the list of fallback non-terminal

pair: terminal pair rule (for alignment composition)
	non-terminal=[defaut non-terminal]

pos: terminal pos rule (for POS annotated input) 

unknown: pos assignment by signature
	signature=[signature for OOV]
	file=[file-name] lexical grammar
	character=[file-name] character model

format: ICU's number/date format rules
	non-terminal=[default non-terminal]
	format=[formtter spec]
	remove-space=[true|false] remove space (like Chinese/Japanese)

For detail see ``cicada --grammar-list``
