// -*- mode: c++ -*-

#ifndef __CICADA__TOKENIZER__HPP__
#define __CICADA__TOKENIZER__HPP__ 1

#include <string>

#include <cicada/sentence.hpp>

namespace cicada
{
  class Tokenizer
  {
  public:
    typedef cicada::Sentence sentence_type;

  public:
    Tokenizer() {}
    virtual ~Tokenizer() {}
    
  private:
    Tokenizer(const Tokenizer& x) {}
    Tokenizer& operator=(const Tokenizer& x) { return *this; }
    
  private:
    struct __true {};
    struct __false {};
    
    template <typename Sent>
    struct is_sentence_type
    {
      static __false value;
    };
    
    template <>
    struct is_sentence_type<senence_type>
    {
      static __true value;
    };

  public:
    template <typename Sent>
    void operator()(const Sent& source, Sent& tokenized) const
    {
      operator()(source, tokenized, is_sentence_type<Sent>::value);
    }
    
    const std::string& algorithm() const { return __algorithm; }
    
  public:
    static Tokenizer&    create(const std::string& parameter);
    static const char* lists();
    
  protected:
    virtual void tokenize(const sentence_type& source, sentence_type& tokenized) const=0;
    
  private:
    template <typename Sent>
    void operator()(const Sent& source, Sent& tokenized, __false) const
    {
      sentence_type __tokenized;
      tokenize(sentence_type(source.begin(), source.end()), __tokenized);
      tokenized = Sent(__tokenized.begin(), __tokenized.end());
    }
    
    template <typename Sent>
    void operator()(const Sent& source, Sent& tokenized, __true) const
    {
      tokenized.clear();
      tokenize(source, tokenized);
    }

    
  private:
    std::string __algorithm;
  };
};

#endif
