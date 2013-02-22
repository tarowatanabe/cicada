MT evaluations especially for MERT purpose

Implemented bootstrap resampling[1] and sign test[2] in "cicada_eval"

BLEU: IBM Bleu (TODO: implement NIST-bleu, average-bleu?) [3]
BLEUS: smoothed BLEU metric[8]
CDER: CDer (wer + shift in one side) [9]
WER: word error rate [4]
InvWER: Inversion word error rate [10]
PER: position independent word error rate [5]
TER: translation error rate [6]
RIBES: RIBES [7]
SB: skip bigram (by default, skip size is clipped to 4) [8]
WLCS: (weighted) longest common subsequence (by defaut, weight is one, meaning no weight, use alpha=2 for weight?) [8]
SK: string kernel (by default, decay factor = 0.8, spectrum p = 4)

cicada_eval [options]

configuration options:
  --tstset arg                 test set file(s)
  --tstset2 arg                test set file(s)   (required for sign test)
  --refset arg                 reference set file(s)
  --output arg (=-)            output file
  --scorer arg (=bleu:order=4) error metric
  --scorer-list                list of error metric
  --signtest                   sign test
  --bootstrap                  bootstrap resampling
  --samples arg                # of samples

reference/test set format is as follows:

segment-id |||  sentence (||| .... some information, such as features etc. ...)

Thus, you can directly feed k-best output from cicada(or cicada_mpi).

Also, you can compute oracle-BLEU from hypergraph(s) by cicada_oracle{,_mpi}

Tips: cicada_oracle{,_mpi} assumes "forests" as a input. If you want to compute
 oracle score for k-bests, you can use cicada_oracle_kbest{,_mpi}.

encode/decode API for scorer:

std::string score.encode();
score_ptr score::decode(std::string);

Evaluator score JSON format:
{"eval":"bleu", "reference":["base64" (length), "base64" (1gram), "base64", "base64", "base64"], "hypothesis":[...]}
{"eval":"bleus", "reference":["base64" (length), "base64" (1gram), "base64", "base64", "base64"], "hypothesis":[...]}
{"eval":"ter", "edits": ["base64", ... (insertion,deletion,substitution,shift,reference)]}
{"eval":"wer", "edits": ["base64", ... (insertion,deletion,substitution,reference)]}
{"eval":"inv-wer", "edits": ["base64", ... (insertion,deletion,substitution,inversion,reference)]}
{"eval":"cder", "edits": ["base64", ... (insertion,deletion,substitution,jump,reference)]}
{"eval":"per", "edits": ["base64", ... (insertion,deletion,substitution,reference)]}
{"eval":"ribes", "distance": "base64", "penalty": "base64"}
{"eval":"sk", "reference":["base64" (match), "base64" (norm)], "hypothesis":["base64", ... (match, norm)]}
{"eval":"sb", "reference":["base64" (match), "base64" (norm)], "hypothesis":["base64", ... (match, norm)]}
{"eval":"wlcs", "reference":["base64" (match), "base64" (norm)], "hypothesis":["base64", ... (match, norm)]}
{"eval":"combined", "score":[{"eval":"bleu", ...}, {"eval":"ter", ...}], "weight":["base64", "base64",.... list of weights]}
and doubles are encoded as base64 (string!). We can insert spaces, but encoder will generate a string w/o spaces.

References:

[1]
@inproceedings{koehn:2004:EMNLP,
  author    = {Koehn, Philipp},
  title     = {Statistical Significance Tests for Machine Translation Evaluation },
  booktitle = {Proceedings of EMNLP 2004},
  editor = {Dekang Lin and Dekai Wu},
  year      = 2004,
  month     = {July},
  address   = {Barcelona, Spain},
  publisher = {Association for Computational Linguistics},
  pages     = {388--395}
}

[2]
@inproceedings{1219906,
 author = {Collins, Michael and Koehn, Philipp and Ku\v{c}erov\'{a}, Ivona},
 title = {Clause restructuring for statistical machine translation},
 booktitle = {ACL '05: Proceedings of the 43rd Annual Meeting on Association for Computational Linguistics},
 year = {2005},
 pages = {531--540},
 location = {Ann Arbor, Michigan},
 doi = {http://dx.doi.org/10.3115/1219840.1219906},
 publisher = {Association for Computational Linguistics},
 address = {Morristown, NJ, USA},
 }

[3]
@InProceedings{papineni-EtAl:2002:ACL,
  author    = {Kishore Papineni  and  Salim Roukos  and  Todd Ward  and  Wei-Jing Zhu},
  title     = {Bleu: a Method for Automatic Evaluation of Machine Translation},
  booktitle = {Proceedings of 40th Annual Meeting of the Association for Computational Linguistics},
  month     = {July},
  year      = {2002},
  address   = {Philadelphia, Pennsylvania, USA},
  publisher = {Association for Computational Linguistics},
  pages     = {311--318},
  url       = {http://www.aclweb.org/anthology/P02-1040},
  doi       = {10.3115/1073083.1073135}
}

[4]
Sonja Nießn, Franz Josef Och, Gregor Leusch, Hermann Ney
"An Evaluation Tool for Machine Translation: Fast Evaluation for MT Research".
In Proc. 2nd International Conference on Language Resources and Evaluation, pp. 39-45, Athens, Greece, May-June 2000

[5]
@INPROCEEDINGS{Tillmann97accelerateddp,
    author = {C. Tillmann and S. Vogel and H. Ney and A. Zubiaga and H. Sawaf},
    title = {Accelerated Dp Based Search For Statistical Translation},
    booktitle = {In European Conf. on Speech Communication and Technology},
    year = {1997},
    pages = {2667--2670}
}

[6]
@INPROCEEDINGS{Snover06astudy,
    author = {Matthew Snover and Bonnie Dorr and Richard Schwartz and Linnea Micciulla and John Makhoul},
    title = {A study of translation edit rate with targeted human annotation},
    booktitle = {In Proceedings of Association for Machine Translation in the Americas},
    year = {2006},
    pages = {223--231}
}

[7]
@InProceedings{isozaki-EtAl:2010:EMNLP,
  author    = {Isozaki, Hideki  and  Hirao, Tsutomu  and  Duh, Kevin  and  Sudoh, Katsuhito  and  Tsukada, Hajime},
  title     = {Automatic Evaluation of Translation Quality for Distant Language Pairs},
  booktitle = {Proceedings of the 2010 Conference on Empirical Methods in Natural Language Processing},
  month     = {October},
  year      = {2010},
  address   = {Cambridge, MA},
  publisher = {Association for Computational Linguistics},
  pages     = {944--952},
  url       = {http://www.aclweb.org/anthology/D10-1092}
}

[8]
@inproceedings{lin-och:2004:ACL,
  author    = {Lin, Chin-Yew  and  Och, Franz Josef},
  title     = {Automatic Evaluation of Machine Translation Quality Using Longest Common Subsequence and Skip-Bigram Statistics},
  booktitle = {Proceedings of the 42nd Meeting of the Association for Computational Linguistics (ACL'04), Main Volume},
  year      = 2004,
  month     = {July},
  address   = {Barcelona, Spain},
  pages     = {605--612},
  url       = {http://www.aclweb.org/anthology/P04-1077},
  doi       = {10.3115/1218955.1219032}
}

[9]
@INPROCEEDINGS{Leusch06cder:efficient,
    author = {Gregor Leusch and Nicola Ueffing and Hermann Ney},
    title = {CDER: Efficient MT Evaluation Using Block Movements},
    booktitle = {In Proceedings of EACL},
    year = {2006},
    pages = {241--248}
}

[10]
@INPROCEEDINGS{Evaluation03anovel,
    author = {Machine Translation Evaluation and Gregor Leusch and Nicola Ueffing and Hermann Ney and Lehrstuhl Fﾃｼr Informatik},
    title = {A Novel String-to-String Distance Measure With Applications to},
    booktitle = {In Proceedings of MT Summit IX},
    year = {2003},
    pages = {240--247}
}