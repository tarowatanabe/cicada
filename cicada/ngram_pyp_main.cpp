#include "ngram_pyp.hpp"
#include "sentence.hpp"
#include "vocab.hpp"

#include "utils/mathop.hpp"

typedef cicada::Sentence sentence_type;
typedef cicada::Vocab    vocab_type;

int main(int argc, char** argv)
{
  if (argc < 2) {
    std::cout << argv[0] << " ngram-pyp-file" << std::endl;
    return 1;
  }
  
  cicada::NGramPYP lm(argv[1]);
  
  sentence_type sentence;
  sentence_type ngram(1, vocab_type::BOS);
  
  double logprob_total = 0.0;
  size_t num_word = 0;
  size_t num_oov = 0;
  size_t num_sentence = 0;
  
  while (std::cin >> sentence) {
    ngram.resize(1);
	  
    sentence_type::const_iterator siter_end = sentence.end();
    for (sentence_type::const_iterator siter = sentence.begin(); siter != siter_end; ++ siter) {
      ngram.push_back(*siter);
      
      const bool is_oov = ! lm.vocab().exists(*siter);
      const double prob = lm.prob(std::max(ngram.begin(), ngram.end() - lm.order()), ngram.end());
      
      if (! is_oov)
	logprob_total += std::log(prob);
      
      num_oov += is_oov;
    }
    
    ngram.push_back(vocab_type::EOS);
    
    const double prob = lm.prob(std::max(ngram.begin(), ngram.end() - lm.order()), ngram.end());
    
    logprob_total += std::log(prob);
    
    num_word += sentence.size();
    ++ num_sentence;
  }
  
  std::cout << "# of sentences: " << num_sentence
	    << " # of words: " << num_word
	    << " # of OOV: " << num_oov
	    << " order: " << lm.order()
	    << std::endl;
  
  std::cout << "logprob = " << logprob_total << std::endl;
  std::cout << "ppl     = " << utils::mathop::exp(- logprob_total / (num_word - num_oov + num_sentence)) << std::endl;
  std::cout << "ppl1    = " << utils::mathop::exp(- logprob_total / (num_word - num_oov)) << std::endl;
}
