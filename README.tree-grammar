Tree-Grammar:

We can use tree-transduce grammar as in:

[X]([X,1] [X](a b c)) ||| [Y](A B [Y,1]) ||| feature1=0.5 feature2=-40

no '(' or ')' allowed in non-terminals or terminals, but you can still "escape" by '\(' '\)'
(You should escape '\' by '\\')

In cicada file-based grammar can be loaded in additin to pre-defined fallback-rules via

   --tree-grammar file-name:feature0="name-of-first-feature", feature1="name-of-second-feature", ...
   --tree-grammar fallback:non-terminal=[x]

For defailt, see --tree-grammar-list
