// -*- mode: c++ -*-

#ifndef __CICADA__TOKENIZER__HPP__
#define __CICADA__TOKENIZER__HPP__ 1

#include <string>

#include <cicada/sentence.hpp>

namespace cicada
{

  struct __tokenizer_true {};
  struct __tokenizer_false {};
  
  template <typename Sent>
  struct __tokenizer_is_sentence_type
  {
    static __tokenizer_false value;
  };
  
  template <>
  struct __tokenizer_is_sentence_type<cicada::Sentence>
  {
    static __tokenizer_true value;
  };
  
  class Tokenizer
  {
  public:
    typedef cicada::Sentence sentence_type;
    
    typedef sentence_type::value_type word_type;

  public:
    Tokenizer() {}
    virtual ~Tokenizer() {}
    
  private:
    Tokenizer(const Tokenizer& x) {}
    Tokenizer& operator=(const Tokenizer& x) { return *this; }

  public:
    template <typename Sent>
    void operator()(const Sent& source, Sent& tokenized) const
    {
      operator()(source, tokenized, __tokenizer_is_sentence_type<Sent>::value);
    }
    
    const std::string& algorithm() const { return __algorithm; }
    
  public:
    static Tokenizer&  create(const std::string& parameter);
    static const char* lists();
    
  protected:
    virtual void tokenize(const sentence_type& source, sentence_type& tokenized) const=0;
    
  private:
    template <typename Sent>
    void operator()(const Sent& source, Sent& tokenized, __tokenizer_false) const
    {
      sentence_type __tokenized;
      tokenize(sentence_type(source.begin(), source.end()), __tokenized);
      tokenized = Sent(__tokenized.begin(), __tokenized.end());
    }
    
    template <typename Sent>
    void operator()(const Sent& source, Sent& tokenized, __tokenizer_true) const
    {
      tokenized.clear();
      tokenize(source, tokenized);
    }

    
  private:
    std::string __algorithm;
  };
};

#endif
