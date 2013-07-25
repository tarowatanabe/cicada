non-terminals
=============

Non-terminals are represented by ``[.+(,[0-9]+)]``.
The ``(,[0-9]+)`` indicates an index to tail nodes, used in a
hypergraph or synchronous rules to indicate mapping of non-terminals
between the source and target sides. Besides the index, non-terminal
can have special flags:

-  If a non-terminal contains ``@`` it is a symbol-refixed non-terminal, and "deciam number" is
   assigned immediately aftr ``@``.
   
-  If a non-terminal contains ``^``, it is a binarized non-terminal so that we can erase
   it to reproduce original hyperedge.
      
   +  ``~`` is used to indicate left-binarized label sequence after ``^``.
   +  ``_`` is used to indicate right-binarized label sequence after ``^``.
   
   +  ``^L`` indicates dependency left-binarization.
   +  ``^R`` indicates dependency right-binarization.
   +  ``^H`` indicates dependency binarization with head.

   +  ``+`` singed indicate concatenation of labels, usually indicating binarized and merged labels.
   +  ``/`` indicates right substitution
   +  ``\`` indidates left substitution

There exists special curly bracketed "terminal" symbols, ``{.+}``, indicating dependency-head.
