Training Sample
===============

The sample data comes from the 30K sentences from Kyoto Free
Translation Task (http://www.phontron.com/kftt/), which was derived
from http://alaginrc.nict.go.jp/WikiCorpus/index_E.html.
The English side of the bilingual data is segmented with Moses
(http://statmt.org/moses/) tokenizer with a simple postprocessing to
match with Penntreebank style tokenizatin. The Japanese side was
segmented by Mecab (https://code.google.com/p/mecab/).
Then, training data with more than 30 words are removed, and took only
the first 30k sentences.

The scripts in the subdirectories are described in
`doc/training.rst` and `doc/traiing-stsg.rst`.

