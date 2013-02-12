Attributes in hypergraph

In addition to features, each edge in hypergraph can hold "attributes," meta data.
The format is:

key=value key=value

where key is non-space string, and value can take either integer, float or JASON string (or an instance of double quoted string).
See hypergraph format doc for possible values.

Actual implementation relies on boost.variant, a compact union storage.

The attributes can be associarted with each of the synchronous rules, like

[x] ||| source-rhs ||| target-lhs ||| features (||| attributes)

for synchronous-CFG and

source-xRS ||| target-xRS ||| features (||| attributes)

if using tree-transducer. Attributes are specified as key=value format as in features format.
Attributes are treated as an optional meta-values from the grammar.
For instance, you can encode lexicalized-reordering rules in attributes for phrasal translations.
