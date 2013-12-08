#include "ngram_nn.hpp"
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
  
  cicada::NGramNN lm(argv[1]);
  
  sentence_type sentence;
  sentence_type ngram(1, vocab_type::BOS);
  
  double logprob_total = 0.0;
  double logprob = 0.0;
  size_t num_word = 0;
  size_t num_oov = 0;
  size_t num_sentence = 0;
  
  while (std::cin >> sentence) {
    ngram.resize(1);

    double logprob_sentence = 0.0;
	  
    sentence_type::const_iterator siter_end = sentence.end();
    for (sentence_type::const_iterator siter = sentence.begin(); siter != siter_end; ++ siter) {
      ngram.push_back(*siter);
      
      const bool is_oov = ! lm.vocab().exists(*siter);
      const double lp = lm(std::max(ngram.begin(), ngram.end() - lm.order()), ngram.end());
      
      std::cerr << "word logprob: " << lp << " oov: " << is_oov << std::endl;
      
      if (! is_oov)
	logprob_total += lp;
      logprob += lp;
      logprob_sentence += lp;
      
      num_oov += is_oov;
    }
    
    ngram.push_back(vocab_type::EOS);
    
    const double lp = lm(std::max(ngram.begin(), ngram.end() - lm.order()), ngram.end());
    
    logprob_total += lp;
    logprob += lp;
    logprob_sentence += lp;
    
    std::cerr << "logprob: " << logprob_sentence << std::endl;
    
    num_word += sentence.size();
    ++ num_sentence;
  }
  
  std::cout << "# of sentences: " << num_sentence
	    << " # of words: " << num_word
	    << " # of OOV: " << num_oov
	    << " order: " << lm.order()
	    << std::endl;
  
  std::cout << "logprob       = " << logprob_total << std::endl;
  std::cerr << "logprob(+oov) = " << logprob << std::endl;
  std::cout << "ppl           = " << utils::mathop::exp(- logprob_total / (num_word - num_oov + num_sentence)) << std::endl;
  std::cout << "ppl1          = " << utils::mathop::exp(- logprob_total / (num_word - num_oov)) << std::endl;
  std::cerr << "ppl(+oov)     = " << utils::mathop::exp(- logprob / (num_word + num_sentence)) << std::endl;
  std::cerr << "ppl1(+oov)    = " << utils::mathop::exp(- logprob / (num_word)) << std::endl;
}
