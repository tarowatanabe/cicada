Training Sample
===============

The sample data comes from the 30K sentences from Kyoto Free
Translation Task (http://www.phontron.com/kftt/), which was derived
from http://alaginrc.nict.go.jp/WikiCorpus/index_E.html.
The sample data used in this service contains English contents which
is translated by the National Institute of Information and
Communications Technology (NICT) from Japanese sentences on
Wikipedia. Our use of this data is licensed by the Creative Comons
Attribution-Share-Alike License 3.0. Please refer to
http://creativecommons.org/licenses/by-sa/3.0/ or
http://alaginrc.nict.go.jp/WikiCorpus/ for details.

The English side of the bilingual data is segmented with Moses
(http://statmt.org/moses/) tokenizer with a simple postprocessing to
match with Penntreebank style tokenizatin. The Japanese side was
segmented by Mecab (https://code.google.com/p/mecab/).
Then, training data with more than 30 words are removed, and took only
the first 30k sentences.

The scripts in the subdirectories are described in
`doc/training.rst` and `doc/traiing-stsg.rst`.

