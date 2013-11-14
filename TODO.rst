Do we allow state-less features applied during composition???
      - We need composition-based features for phrase-based MT
	to compute distortion/lexicalized reordering etc. not represented in
	hypergraph struct

      - State-less features are sometimes related to composition algorithm,
	such as shift-reduce etc.
	
      - We have "parse" variant of "compose". Do we allow stateless features during "parsing" ????

Systematic symbols:
    - So forth, we have abused many symbols, such as ``|:;+`` for
      binarization, permutation, parent, etc.
    - I think it's time to systemacically define symbols so that we
      can de-binarize, de-parentize etc.

Add head-finder?
    - Head finder requires child/parent node, therefore, requires the
      whole hypergraph...

Re-engineer kbest/viterbi framework...
  - I don't like current interface, separating weight and derivation.
  - Is there a workaround, like use iterator-based enumeration?

Support word alignment training!
  - We know feature-rich unsupervised training!
  - Use of syntax-tree combined with alignment
  - Cube-pruning to transform source-syntactic tree into alignment-derivation! (use of "alignment" as our output...)
  - Can we integrate supervised learning and unsupervised learning sharing the same code...?

      We can potentially share by modifying "intersection"

  - How to encode lexical feature...?

    use grammar! (source target word mapping with features!)

  - Intersection with the supplied alignment: supervised training!

    Can we do this, given that word alignment space is rather constrained...

Implement scalable training by L1 regularized feature selection:
  1. Collect features.. it is the hardest part...
     Decode-and-extract features
     or use simple featueres (word-pairs etc) from bitext?
  2. Compute importance
  3. Learn w/o introducing "new features": modify not to introduce new features

  Remarks: what kind of features?
     features easily collectable from bilingual data (with monolingual data), (pairs, ngrams etc.)
     or, 
     features requing complex hidden variables...?

POS annotated sentence input:
  - How to convert each POS into latent annotated POS?

    - Solution 1. Implemenet GrammarUnknownPOS which will give
      annotated POS when already in the model.
      Otherwise use POS as a "hint"?

Add more tests

A unified framework to collect statistics: Do we need stats from features?
    - We have a support to collect stats from operations, but not from
      individual features..

Cube summing?
     - I'm still not sure how to efficiently collect residuals... Do
       we include state-less features to serve as residual counts?

Add incremental version of MERT?

Correctly implement cutting-plane algorithms:
   - We need to merge feature vector into a bundle, then, iterative
     add new vectors...
   - Do we implement in maxlike families...?

Add MPI version of Model1/HMM alignment?

Model1 starting from dice-estimated parameters...??

Dictionary constrained Model1/HMM training?
   - As in GIZA++, we will constraint by the existence of dictionary
     items

Alignment constrained alignment
   - As observed in training, we will constraint the alignment space
     in training.
   - Do we also constrain for the viterbi/ITG/MaxMatch alignment?

Integrate cicada-alignment.py and cicada-extract.py???
   - Currently, NO, in order to encourange user to select an
     alternative paths.

Integrate MPI and non-MPI of cicada_learn{,kbest}{,mpi}
   - Use the same gradient,margin computers across mpi and non-mpi
     applications

Global lexicon learning by liblinear...?
   - We will add bias feature, then learn...

Collapse and binarize the source side of the ghkm/tree grammar
   - Perform "collapsing" and do binarization (either source or target)
   - Perform lexicalization?
   - Perform synchronous binarization? (we need to store all in
     memory...?) 
       
      And left-to-right and target-left-to-right binarization

   - How to represent intermedidate non-terminals?

      We need to keep track of both-side....???
     
Do we keep pair of source-target lhs???
  - I'm pretty not sure how symbols should be handled during the CKY
    algorithm.
  - Probably, we need to keep both side and maintain it during unary
    rules...?

Remove sub-tree sharing in compose-tree-cky and parse-tree-cky?
  - I'm not sure which is better..

Support joshua training?

Support cicada forest converter?

Implement succinct-rtrie:
  - This is a trie-like structure, but it's purpose is a reverse
    index: given id, query data of string.
  - strings are stored in a trie, with an index refering to the leaf
    position of trie of id.
  - string is uncovered by leaf-to-parent traversal (thus, indexing
    requires reversing, first)

remote grammar and tree-grammar:
  - set up server for grammar(s) and query from clients
  - performs composition/parse/generation at the server side with
    larger memory (meaning with larger grammar)
  
revise static create interface:
  - currently, we return reference, but it is safer to return
    pointer...?

Earley composition with skipping
   - Like phrase-based SMT, allow local skipping... very hard...

Earley composition with CFG!
   - Like phrasal composition, we perform CKY over the Earley
     generated forest of string.... very hard...

optimized variant of SGD (and xBLEU?)
  - xBLEU is impossible given that the combination is already
    weighted!

Parallel learning for PYP:
   - translit and segment

Add a shallow hiero rule by limiting the "depth" of rule instantiation...
    - level one is easier, but how to handle arbitrary depth? (or, at
      least, depth of 2?)

Revise syntactic alignment:
  - Binarize before processing
  - Use of pialign style bi-parsing algorithm
  - Allow arbitrary sub-tree alignment by a CFG-style sub-tree
    transformation
  - Paired with "phrases" in the target side.

An analysis tool based on error metrics?
 - Compute an error metric, i.e. BLEU, and visualize matched portions
   (such as ngrams, alignment etc.)
 - TOOD: API?

Mixture PYP-LM:
 - We employ multi-floor CRP for representing mixture of multiple LM.
 - Multiple: surface, prefix-4 and suffix-4! (+ class-LM or +POS-LM?)

Revise ngram-pyp so that we do not have to re-compute lower-order probabilities...

Implement ADMM for potentially better parallel training

Better sharing python code... HOW?

mpipe is buggy under mac osx...
  - This is probaly because of the interaction between openmpi and
    fork() with unmanagable file descriptors etc. The bug is clearly
    exhibited by the difference of the # of lines read and the actual
    read from stdin!

cicada_filter_kbest to support {file,directory}-to-{file,directory}

add error checking for codecs

msgpack for succinct storage

Rearrange learning code by splitting L2 and/or L1 projection
  - Especially for online-learning, they can be set up as additional
    "unified" code
  - Split learning rate scheduling from the code so that we can choose
    from exponential decay, adagrad etc.

Non-linear features
   - The decoder uses non-linear combination with hidden layers... How
     to implement on cicada?
   - Use a special dot-product function which project all the features
     to the hidden layers, then, perform combination

     - Training should perform hyper-edge-wise, not a simple sentence-wise training...
     - Thus, we will dump forest, and perform forest-wise training, or
       dump k-best trees represented as a set of hypergraphs.

Transform matrix into column-major to avoid confusion with fortran/BLAS etc. and use of "standard" matrix
  - Eigen, armagillo, vienna-cl?
  - I think eigen is easier since it can plug by copying headers...

Use of "float" not "double" for better integration with GPUs
  - And people will not care such precision...

Use of eigen for weight vector maintenance... (for potentially faster computation...)

Revise the neuron modules to support computation with more dimensions (not assuming single dimension...)
  - Add "embedding" structure to hold word-embedding features (a float-vector!)
  - Do we add "loss" module to compute losses?
  - Do we add "training" module?e

Implement LBFGS/CG by templates since this may conflict with float/double based implementation of liblbfgs

Unify the liblbfgs and cg_descent code...?

Lua integration:
  - Any use...????

Correctly implement alignment/distortion model estimation in lexicon_hmm/model4
  - Currenlty, it is very hacky, and gives non-optimal parameters...

Implement PCFG language model with arbitrary order
  - Format: a hyperedge is compactly represented following the
    tree-rule format(?) [lhs]([rhs],word,\,,word2,\))
  - Or, use `|||` as a separator?
  - Use of the tree-grammar way to encode tree structure
    * LHS is represented as a single node
    * A vector of children is represented as a single node

Diversified kbest/feature application + Sampled kbest/feature application

Port for other batch-queue systems, such as Torque.

Implement rejection sampling
 - How to approximately sample during non-local feature application?

Merge learning codes
 - There exists duplicated and experimental codes which I don't even
   remember.
 - Batch-based learner should follow online-learning style regularization.

Add head-node annotator for each edge:
 - This can be used to identify the pruning bin, so that we do not
   rely on the "prune-bin" attribute assigned during feature
   application.

NCE for autoencoding
- re-implement by perform classification and sampling using the lexicon model

Neural Network Alignment model
- re-implement by summation and Viterbi approximation
