Phrase/Rule/Tree-fragment extraction and indexing

We support various phrase/rule/tree extraction. Basically, we will run in two steps:
   Collect counts by mapping input sentences/parse trees/alignment etc.
   Map/Reduce to collect target-side counts
   Map/Reduce to collect source/target counts and produce various counts

   They can operate similarly with underlying extraction specific map/reduce framework.

You can try "cicada-extract.py" which runs like phrase-extract in Moses, but starting from step 4 (and 6).
It is recommended to specify TMPDIR_SPEC (or TMPDIR) pointing to 
a large "local" storage for efficient extraction:

export TMPDIR_SPEC=<directory for temporary storage, better located at local disk space, not on NFS>

and/or specify --temporary-dir.

Then, you should run "cicada-index.py" to convert collected counts into a set of feature values, and then,
into a binary format for use with the cicada toolkit.

Internally, cicada-extract.py calls following binaries:

cicada_extract_phrase{,_mpi}
	Extract phrase
	Output scores:
	       lhs ||| rhs ||| count(lhs, rhs) \
	       	   count(prev, mono, lhs, rhs) \
		   count(prev, swap, lhs, rhs) \
		   count(next, mono, lhs, rhs) \
		   count(next, swap, lhs, rhs)

cicada_extract_rule{,_mpi}
	Extract synchronous-CFG + syntax augmentation (aka SAMT) when extracted with "span" data.
	You can generate span by "cicada_filter_penntreebank"
	Output scores:
	       root lhs ||| root rhs ||| count(lhs, rhs)

cicada_extract_ghkm{,_mpi}
	Extract tree-to-string or string-to-tree rules by GHKM
	Output scores:
	       lhs-xRS ||| rhs-xRS ||| count(lhs, rhs)

cicada_extract_tree{,_mpi}
	Extract tree-to-tree rules by GHKM
	Output scores:
	       lhs-xRS ||| rhs-xRS ||| count(lhs, rhs)

After counts collection, you can summarize them by cicada_extract_counts{,_mpi}
Which will output:

lhs ||| rhs ||| alignments ||| counts(lhs, rhs) ||| counts(lhs) ||| counts(rhs) ||| observed(lhs) observed(rhs) 

In addtion, we will dump root-joint.gz, root-source.gz and root-target.gz. root-joint.gz contains:

root(lhs)root(rhs) ||| counts(root(lhs)root(rsh)) ||| observed(lhs, rhs)

while root-source.gz looks like:

root(lhs) ||| counts(root(lhs)) ||| observed(lhs)

while root-target.gz looks like:

root(rhs) ||| counts(root(rhs)) ||| observed(rhs)

You can easily transform the counts into probabilities by maximum likelihood estimates, or use observed counts
to perform Dirichlet prior smoothing (default) by running "cicada-index.py".
The cicada-index.py transforms collected counts into a set of feature values, then, encodes the grammar into a binary format.
Internally, the indexer calls:

   cicada_filter_extract:
	extract only nbet of target variation for each source side,
   	measured by its joint frequency of lhs and rhs.

   cicada_filter_extract_phrase:
	Dump in moses or cicada format. Also, you can dump
	lexicalied reordering table.
   				 
   cicada_filter_extract_scfg:
	Dump in cicada format for synchronous-CFG. You can also 
	add features for lhs given root(lhs) and rhs given root(rhs)

   cicada_filter_extract_ghkm:
	Dump in cicada format for tree-to-string, string-to-tree, tree-to-tree.

and calls
    cicada_index_grammar
	for indexing phrase/scfg
    
    or

    cicada_index_tree_grammar
	for indexing tree-to-string, string-to-tree, tree-to-tree rules

    
Remarks on SCFG: (I don't know the default for Joshua... any clues?)
	If you want to simulate the Moses-style hierarchical rule, you can use "cicada-extract.py" with options:
	
	--max-span-source 15
	--max-span-target 15
	--min-hole-source 2
	--min-hole-target 1
	--max-length 5
	--max-fertility 0
	--exhaustive // exhaustively extract from all-possible phrases

	The Hiero-style extraction from Chiang (2007) will be:

	--max-span-source 10
	--max-span-target 0
	--min-hole-source 2
	--min-hole-target 1
	--max-length 5(?)
	--max-fertility 0
	--constrained // extract only from minimal phrases

	Our default:

	--max-span-source 15
	--max-span-target 15
	--min-hole-source 1
	--min-hole-target 1
	--max-length 7
	--max-fertility 4
	//
	// no --exhaustive, and no --constrained imply:
	//   extract from all-possible phrases, but uses minimal phrases as "holes"
	//

Remarks on GHKM:
	--max-nodes 15
	--max-height 4
	--max-compose 0 // no constraint for composition from minimal rules
	--max-scope 0   // no constriant for maximum scopes
	
	// other options...	
	--constrained  // even minimum rules should satisfy above max-constraints.
	               // RECOMMENDED when extracting rules from "forest"
	--exhaustive   // consider all possible attachment of unaligned words to any rules... NOT RECOMMENDED.
	               // default will attach unaligned words to rules "closer to root"

	--collapse-source
	--collapse-target
	               // treat source/target side as a "flat" rule. Recommended for string-to-x variant.
	--project      // project non-terminals symbols from tree-side to string-side.