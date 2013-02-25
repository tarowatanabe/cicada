# Grammar for cicada

string-to-string rules (aka Hiero rules):

[lhs] ||| [x,1] abc bb ||| ddd [x,1] ||| feature=1.0 feature2=5.6 ||| attr1=1.6
[lhs] ||| [x,1] abc [x,2] bb ||| ddd [x,1] fff [x,2] ||| feature=-1 feature2=8 ||| attr2="bad"

cicada assume that rule's non-terminal indices are sorted wrt the source side, not target side.

tree-to-string/string-to-tree transduction rules:

[LHS]([NN,1](terminal [NP](xxx) [NP,2])) ||| [LHS]([VP](yyy [x,2] zzz) [x,1]) ||| feature=6.8 ||| attr2="attribute 2"

In cicada, various file-based gramamr can be loaded in addition to pre-defined glue-rules, by
   
	load a grammar from "file-name"
	--grammar file-name:max-span="maximum-span",feature0="name-of-first-feature", feature1="name-of-second-feature", ...

	glue-rule
	--grammar glue:straight=ture,inverted=true,non-terminal=[x],goal=[s],fallback=<fallback non-terminal symbol list>
	
	insertion-rule
	--grammar insertion:non-terminal=[x]

	deletion-rule
	--grammar deletion:non-terminal=[x]
		
	word pair grammar for word alignment
	--grammar pair:non-terminal=[x]
	
	ICU date formatter grammar
	--grammar format:non-terminal=[x],format="date:locale-source=en,locale-target=ja", remove-space={true or false}
	
	ICU number formatter grammar
	--grammar format:non-terminal=[x],format="number:locale-source=en,locale-target=ja", remove-space={true or false}

	unknown rule for used by monolingual parser
	--grammar unknown:signature=English,file=[grammar spec],character=[character model]

	pos rule for pos-annotated input
	--grammar pos

For detail see --grammar-list
