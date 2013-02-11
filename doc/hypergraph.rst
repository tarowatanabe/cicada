==========
hypergraph
==========

-------------------------------------------------------------------------
hypergraph foramt description and tools to convert from/to the hypergraph
-------------------------------------------------------------------------

:Author: Taro Watanabe <taro.watanabe@nict.go.jp>
:Date:   2013-2-11

FORMAT
------

The hypergraph, or forest in short, is represented by a `JSON data format <http://www.json.org>`_.
Strings must be escaped (see JSON specification). You may insert spaces at arbitrary positions.
One-line-per-single-hypergraph is prefered for easier preprocessing.
We assume topologically sorted hypergraph (or, node-id is ordered by post-traversal order), but
input/output can handle node ids in any orders.

HYPERGRAPH

  ::

    {
     "rules" : [ RULE, RULE, ... ], 
     "nodes" : [ NODE, NODE, ...],  
     "goal" : node-id-for-goal (optional)
    }

  The "rules" field is a list of rules associated with each edges. Each
  rule-id starts from one, and refered by its rule-id by each
  hyperedge. zero is reserved for no-rule, or errornous hypergraph.

  The "nodes" field is a list of nodes and each node is refered by its
  node id. The node-id starts from zero.

  The "goal" field points to the final goal node id of the
  hypergrpah. Unorder the topologically sorted order, the goal field
  usually points to the last node in the node list.
  No goal implies an invalid hypergraph.

RULE

  ::

    "LHS ||| RHS"

  Each rule is a JSON string, thus some characters should be escaped, i.e. " to \\", \\ to \\\\ etc.
  "LHS" is a left-hand-side of a rule. "RHS" is a list of symbols,
  either terminals or non-terminals. Note that the "LHS" is used as a label for
  each node.

NODE

  ::

    [EDGE, EDGE, ...]

  Each node consists of a list of edge.

EDGE

  :: 

    {
     "tail" : [node-id, node-id, ...], 
     "feature" : {"feature-name" : FLOAT, "feature-name2" : FLOAT, ... }, (optional)
     "attribute" : {"attribute-name" : ATTRIBUTE, "attribute-name2" : ATTRIBUTE, ...}, (optional)
     "rule" : rule-id (zero indicating no-rule == invalid hypergraph)
    }

  The "tail" field is a list of node in a hypergraph.
  The "feature" field is a list of key-value pair, consisting of
  JSON string and FLOAT value.
  The "attribuete" field is a list of key-value pair, consisting of
  JSON string and ATTRIBUTE value. The ATTRIBUTE can take either 64bit
  integer, floaring point value (double precision) or JSON string.
  The "rule" field is a map to a rule id in the rule list.

ATTRIBUTE

  :: 

    FLOAT(double precision) | INTEGER(64-bit) | JSON-string

  As described in the EDGE, ATTRIBUTE can be either 64bit integer,
  floating point value (double precision) or JSON string. Internally,
  it is implemented as `boost.variant <http://www.boost.org/doc/libs/release/libs/variant/>`_ which
  supports "enum" like storage in an efficient fashion.

TOOLS
-----

cicada_filter_penntreebank

  A tool which transform Penn Treebank style constituency parse
  tree(s) into JSON hypergrpah format.

cicada_filter_dependency

  A tool which transforms dependency trees into a JSON hypergraph
  format. Currently, we support: MST, CoNLL, Malt, Cabocha and cicada
  native format.

ciada_filter_charniak

  A tool which transforms Charniak's parser forest output into a JSON
  hypergraph format.
