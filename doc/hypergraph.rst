hypergraph
==========

JSON hypergraph foramt description and tools to convert from/to the hypergraph.

Format
------

The hypergraph, or forest in short, is represented by the `JSON data format <http://www.json.org>`_.
Strings must be escaped (see JSON specification). You may insert
spaces or newlines at arbitrary positions, but
one-line-per-single-hypergraph is prefered for easier preprocessing.
We assume topologically sorted hypergraph (or, node-id is ordered by post-traversal order), but
input/output can handle node ids in any orders.

.. code::

  hypergraph   ::= {"rules" : [ rule *(, rule)],
	            "nodes":  [ node *(, node)],
		    "goal": JSON-INT}

  rule         ::= "non-terminal ||| +(non-terminal | terminal)"
  non-terminal ::= [.+] | [.+,[0-9]+] (json escaped string)
  terminal     ::= .+                 (json escaped string)

  node         ::= [edge *(, edge)]
  edge         ::= {"tail" : [ JSON-INT *(, JSON-INT) ],
                    "feature" : features,
                    "attribute" : attributes,
                    "rule" : JSON-INT}

  features     ::= {JSON-STRING : JSON-FLOAT *(, JSON-STRING : JSON-FLOAT)}
  attributes   ::= {JSON-STRING : value *(, JSON-STRING : value)}
  value        ::= JSON-STRING | JSON-FLOAT | JSON-INT


The hypergraph consists of three fields.
The "rules" field is a list of rules associated with each edges.
Each rule-id starts from one, and refered by its rule-id by each
hyperedge.
Zero is reserved for no-rule, or errornous hypergraph.
The "nodes" field is a list of nodes and each node is refered by its
node id.
The node id starts from zero.
The "goal" field points to the final goal node id of the hypergrpah.
Unorder the topologically sorted order, the goal field usually points
to the last node in the node list.
No goal implies an invalid hypergraph.

Each rule is a JSON string, thus some characters employed in rhs
and/or lhs should be escaped, i.e. " to \", \ to \\ etc.
The lhs of a rule consists of a single non-terminal, whereby its rhs
consists of a sequence of non-terminal or terminal.
Note that the lhs is also used as a label for each node.
non-terminal is represented by [.+], and terminal can be any strings
which does not starts with '[' and does not end with ']'.
Also, we cannot use '|||' as a terminal.
Each non-terminal in the rhs can embed an "index" separeted by comma,
which is represented by [.+,([0-9]+)].
The index is an indice to tail nodes (one-base).
For instance, [x,2] refers to the second position of a tail in a
hyperedge.

Each node consists of a list of edge.
Each edge consists of four fields:
The "tail" field is a list of node in a hypergraph.
No tails implies "source nodes."
The optional "feature" field is a list of key-value pair, consisting
of JSON string and FLOAT value.
The optional "attribuete" field is a list of key-value pair,
consisting of JSON string and attribute value.
The attribute value can take either 64bit integer, floaring point
value (double precision) or JSON string.
Internally, it is implemented as
`boost.variant <http://www.boost.org/doc/libs/release/libs/variant/>`_
which supports "enum" like storage in an efficient fashion.
The "rule" field is a map to a rule id in the rule list.

Example
-------

.. code:: json

  {"rules": ["[PRP] ||| I",
           "[NP] ||| [PRP]",
           "[MD] ||| 'd",
           "[VB] ||| like",
	   "[TO] ||| to",
	   "[VB] ||| have",
	   "[DT] ||| a",
	   "[NN] ||| glass",
	   "[NP] ||| [DT] [NN]",
	   "[IN] ||| of", 
	   "[NN] ||| water", 
	   "[NP] ||| [NN]", 
	   "[PP] ||| [IN] [NP]", 
	   "[NP] ||| [NP] [PP]", 
	   "[VP] ||| [VB] [NP]", 
	   "[VP] ||| [TO] [VP]", 
	   "[S] ||| [VP]", 
	   "[VP] ||| [VB] [S]",
	   "[VP] ||| [MD] [VP]", 
	   "[.] ||| .", 
	   "[S] ||| [NP] [VP] [.]", 
	   "[ROOT] ||| [S]"],
    "nodes": [[{"rule":1}],
           [{"tail":[0],"rule":2}],
	   [{"rule":3}],
	   [{"rule":4}], 
	   [{"rule":5}], 
	   [{"rule":6}], 
	   [{"rule":7}], 
	   [{"rule":8}],
	   [{"tail":[6,7],"rule":9}], 
	   [{"rule":10}],[{"rule":11}],
	   [{"tail":[10],"rule":12}], 
	   [{"tail":[9,11],"rule":13}],
	   [{"tail":[8,12],"rule":14}], 
	   [{"tail":[5,13],"rule":15}],
 	   [{"tail":[4,14],"rule":16}], 
	   [{"tail":[15],"rule":17}],
	   [{"tail":[3,16],"rule":18}], 
	   [{"tail":[2,17],"rule":19}],
	   [{"rule":20}], 
	   [{"tail":[1,18,19],"rule":21}],
	   [{"tail":[20],"rule":22}]],
    "goal": 21}

Tools
-----

`cicada_filter_penntreebank`

  A tool which transform Penn Treebank style constituency parse
  tree(s) into JSON hypergrpah format.

`cicada_filter_dependency`

  A tool which transforms dependency trees into a JSON hypergraph
  format. Currently, we support: MST, CoNLL, Malt, Cabocha and cicada
  native format.

`ciada_filter_parseforest`

  A tool which transforms parser forest output into a JSON
  hypergraph format.

`cicada_unite_forest`

  A tool to merge multiple hypergraphs into one. If the label of goal
  nodes differ, then, we will introduce an additional goal node,
  [goal].
