//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "decode.hpp"
#include "encode.hpp"

#include <sstream>
#include <iterator>
#include <queue>

#include "inv_wer.hpp"

#include <boost/functional/hash.hpp>

#include "utils/hashmurmur.hpp"
#include "utils/bichart.hpp"
#include "utils/chart.hpp"
#include "utils/dense_hash_set.hpp"

namespace cicada
{
  namespace eval
  {

    std::string InvWER::description() const
    {
      std::ostringstream stream;
      stream << "inv-wer: " << score()
	     << ' ' << insertion << '|' << deletion << '|' << substitution << '|' << inversion
	     << " length: " << references;
      
      return stream.str();
    }
    
    std::string InvWER::encode() const
    {
      std::ostringstream stream;
      stream << '{' << "\"eval\":\"inv-wer\",";
      stream << "\"edits\":[";
      stream << escaper(insertion)
	     << ',' << escaper(deletion)
	     << ',' << escaper(substitution)
	     << ',' << escaper(inversion)
	     << ',' << escaper(references);
      stream << "]}";
      return stream.str();
    }

    typedef boost::fusion::tuple<double, double, double, double, double> inv_wer_parsed_type;
    
    template <typename Iterator>
    struct inv_wer_parser : boost::spirit::qi::grammar<Iterator, inv_wer_parsed_type(), boost::spirit::standard::space_type>
    {
      inv_wer_parser() : inv_wer_parser::base_type(inv_wer_parsed)
      {
	namespace qi = boost::spirit::qi;
	namespace standard = boost::spirit::standard;
	
	inv_wer_parsed %= (qi::lit('{')
			   >> qi::lit("\"eval\"") >> qi::lit(':') >> qi::lit("\"inv-wer\"") >> qi::lit(',')
			   >> qi::lit("\"edits\"") >> qi::lit(':')
			   >> qi::lit('[')
			   >> double_value >> qi::lit(',')
			   >> double_value >> qi::lit(',')
			   >> double_value >> qi::lit(',')
			   >> double_value >> qi::lit(',')
			   >> double_value
			   >> qi::lit(']')
			   >> qi::lit('}'));
      }
      
      typedef boost::spirit::standard::space_type space_type;
      
      double_base64_parser<Iterator> double_value;
      boost::spirit::qi::rule<Iterator, inv_wer_parsed_type(), space_type> inv_wer_parsed;
    };

    Score::score_ptr_type InvWER::decode(std::string::const_iterator& iter, std::string::const_iterator end)
    {
      typedef std::string::const_iterator iterator_type;
      
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      inv_wer_parser<iterator_type> parser;
      inv_wer_parsed_type           parsed;
      
      const bool result = qi::phrase_parse(iter, end, parser, standard::space, parsed);
      if (! result)
	return score_ptr_type();
      
      std::auto_ptr<InvWER> inv_wer(new InvWER());
      inv_wer->insertion    = boost::fusion::get<0>(parsed);
      inv_wer->deletion     = boost::fusion::get<1>(parsed);
      inv_wer->substitution = boost::fusion::get<2>(parsed);
      inv_wer->inversion    = boost::fusion::get<3>(parsed);
      inv_wer->references   = boost::fusion::get<4>(parsed);
      
      return score_ptr_type(inv_wer.release());
    }
    
    Score::score_ptr_type InvWER::decode(const utils::piece& encoded)
    {
      std::string::const_iterator iter(encoded.begin());
      std::string::const_iterator iter_end(encoded.end());
      
      return decode(iter, iter_end);
    }

    class InvWERScorerImpl
    {
    private:
      friend class InvWERScorer;
      
    public:
      typedef InvWERScorer inv_wer_scorer_type;
      typedef inv_wer_scorer_type::weights_type weights_type;

      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;
      
      typedef cicada::Sentence sentence_type;
      typedef cicada::Symbol   word_type;
      typedef cicada::Matcher  matcher_type;
      
      struct Score
      {
	double score;
	
	double insertion;
	double deletion;
	double substitution;
	double inversion;
	
	Score() : score(std::numeric_limits<double>::infinity()), insertion(0), deletion(0), substitution(0), inversion(0) {}
	
	Score(const double& __score,
	      const double& __insertion,
	      const double& __deletion,
	      const double& __substitution,
	      const double& __inversion)
	  : score(__score), insertion(__insertion), deletion(__deletion), substitution(__substitution), inversion(__inversion) {}

	Score& operator+=(const Score& x)
	{
	  score        += x.score;
	  insertion    += x.insertion;
	  deletion     += x.deletion;
	  substitution += x.substitution;
	  inversion    += x.inversion;
	  
	  return *this;
	}
	
	friend
	Score operator+(const Score& x, const Score& y)
	{
	  Score tmp(x);
	  tmp += y;
	  return tmp;
	}
      };
      
      typedef Score value_type;
      
      struct span_type
      {
	typedef size_t    size_type;
	typedef ptrdiff_t difference_type;
    
	size_type first;
	size_type last;
    
	span_type() : first(0), last(0) {}
	span_type(const size_type& __first, const size_type& __last)
	  : first(__first), last(__last) { }

	bool empty() const { return first == last; }
	size_type size() const { return last - first; }
    
	friend
	bool operator==(const span_type& x, const span_type& y)
	{
	  return x.first == y.first && x.last == y.last;
	}
    
	friend
	size_t hash_value(span_type const& x)
	{
	  typedef utils::hashmurmur<size_t> hasher_type;
      
	  return hasher_type()(x.first, x.last);
	}
      };

      struct span_pair_type
      {
	typedef span_type::size_type      size_type;
	typedef span_type::difference_type difference_type;

	span_type source;
	span_type target;
    
	span_pair_type() : source(), target() {}
	span_pair_type(const span_type& __source, const span_type& __target)
	  : source(__source), target(__target) {}
	span_pair_type(size_type s1, size_type s2, size_type t1, size_type t2)
	  : source(s1, s2), target(t1, t2) {}
    
	bool empty() const { return source.empty() && target.empty(); }
	size_type size() const { return source.size() + target.size(); }
    
	friend
	bool operator==(const span_pair_type& x, const span_pair_type& y)
	{
	  return x.source == y.source && x.target == y.target;
	}
    
	friend
	size_t hash_value(span_pair_type const& x)
	{
	  typedef utils::hashmurmur<size_t> hasher_type;
	  
	  return hasher_type()(x);
	}
      };

      typedef utils::bichart<value_type, std::allocator<value_type> > chart_type;
      
      typedef std::vector<span_pair_type, std::allocator<span_pair_type> > span_pair_set_type;
      typedef std::vector<span_pair_set_type, std::allocator<span_pair_set_type> > agenda_type;
      
      typedef utils::chart<double, std::allocator<double> > chart_mono_type;
      
      typedef std::pair<double, span_pair_type> score_span_pair_type;
      typedef std::vector<score_span_pair_type, std::allocator<score_span_pair_type> > heap_type;

      typedef std::vector<double, std::allocator<double> > alpha_type;
      typedef std::vector<double, std::allocator<double> > beta_type;

      // sort by greater so that we can pop from a small-cost item.
      struct heap_compare
      {
	bool operator()(const score_span_pair_type& x, const score_span_pair_type& y) const
	{
	  return x.first > y.first;
	}
      };
      
      InvWERScorerImpl() { }
      InvWERScorerImpl(const sentence_type& __ref) : ref(__ref) {}
      
      value_type operator()(const sentence_type& sentence, const weights_type& weights, const matcher_type* matcher)
      {
	const sentence_type& source = ref;
	const sentence_type& target = sentence;
	
	const difference_type T = source.size();
	const difference_type V = target.size();
	const difference_type L = T + V;
	
	if (T == 0) return value_type(V * weights.insertion, V * weights.insertion, 0, 0, 0);
	if (V == 0) return value_type(T * weights.deletion, 0, T * weights.deletion, 0, 0);
	
	const double infinity = std::numeric_limits<double>::infinity();
	
	// initialization...
	chart_type      chart(T + 1, V + 1, value_type());
	agenda_type     agenda(L + 1);
	chart_mono_type chart_source(T + 1, infinity);
	chart_mono_type chart_target(V + 1, infinity);
	
	alpha_type alpha_source(T + 1, infinity);
	alpha_type alpha_target(V + 1, infinity);
	beta_type  beta_source(T + 1, infinity);
	beta_type  beta_target(V + 1, infinity);

	// insertion...
	for (difference_type t = 0; t <= T; ++ t) {
	  double score = 0.0;
	  for (difference_type v = 1; v <= V; ++ v) {
	    score += weights.insertion;
	    
	    chart(t, t, 0, v).score     = score;
	    chart(t, t, 0, v).insertion = score;

	    chart_target(0, v) = std::min(chart_target(0, v), score);

	    agenda[v].push_back(span_pair_type(t, t, 0, v));
	  }
	}
	
	// deletion...
	for (difference_type v = 0; v <= V; ++ v) {
	  double score = 0.0;
	  for (difference_type t = 0; t <= T; ++ t) {
	    score += weights.deletion;
	    
	    chart(0, t, v, v).score    = score;
	    chart(0, t, v, v).deletion = score;
	    
	    chart_source(0, t) = std::min(chart_source(0, t), score);
	    
	    agenda[t].push_back(span_pair_type(0, t, v, v));
	  }
	}
	
	// one-to-one...
	for (difference_type t = 0; t != T; ++ t)
	  for (difference_type v = 0; v != V; ++ v) {
	    const double score = (source[t] == target[v] ? 0.0
				  : (matcher && matcher->operator()(source[t], target[v])
				     ? weights.match
				     : weights.substitution));
	    
	    chart(t, t + 1, v, v + 1).substitution = score;
	    chart(t, t + 1, v, v + 1).score        = score;
	    
	    chart_source(t, t + 1) = std::min(chart_source(t, t + 1), score);
	    chart_target(v, v + 1) = std::min(chart_target(v, v + 1), score);
	    
	    agenda[2].push_back(span_pair_type(t, t + 1, v, v + 1));
	  }
	
	// forward-backward to compute estiamtes...
	forward_backward(chart_source, alpha_source, beta_source);
	forward_backward(chart_target, alpha_target, beta_target);
	
	// forward...
	heap_type heap;
	double beam = 2.0;
	
	do {
	  for (difference_type l = 1; l != L; ++ l) {
	    span_pair_set_type& spans = agenda[l];
	    
	    heap.clear();
	    heap.reserve(spans.size());
	    
	    span_pair_set_type::const_iterator siter_end = spans.end();
	    for (span_pair_set_type::const_iterator siter = spans.begin(); siter != siter_end; ++ siter)  {
	      const double cost = (chart(siter->source.first, siter->source.last, siter->target.first, siter->target.last).score
				   + std::max(alpha_source[siter->source.first] + beta_source[siter->source.last],
					      alpha_target[siter->target.first] + beta_target[siter->target.last]));
	      
	      heap.push_back(score_span_pair_type(cost, *siter));
	      std::push_heap(heap.begin(), heap.end(), heap_compare());
	    }
	    
	    heap_type::iterator hiter_begin = heap.begin();
	    heap_type::iterator hiter       = heap.end();
	    heap_type::iterator hiter_end   = heap.end();

	    const double cost_threshold = hiter_begin->first + beam;
	    for (/**/; hiter_begin != hiter && hiter_begin->first < cost_threshold; -- hiter)
	      std::pop_heap(hiter_begin, hiter, heap_compare());

	    const double cost_str = 0.0;
	    const double cost_inv = weights.inversion;
	    
	    for (heap_type::iterator iter = hiter ; iter != hiter_end; ++ iter) {
	      const span_pair_type& span_pair = iter->second;
	      
	      const difference_type s = span_pair.source.first;
	      const difference_type t = span_pair.source.last;
	      const difference_type u = span_pair.target.first;
	      const difference_type v = span_pair.target.last;
	      
	      for (difference_type S = utils::bithack::max(s - l, difference_type(0)); S <= s; ++ S) {
		const difference_type L = l - (s - S);
		
		// straight
		for (difference_type U = utils::bithack::max(u - L, difference_type(0)); U <= u - (S == s); ++ U) {
		  // parent span: StUv
		  // span1: SsUu
		  // span2: stuv
		  
		  if (chart(S, s, U, u).score == infinity) continue;
		  
		  const span_pair_type  span1(S, s, U, u);
		  const span_pair_type& span2(span_pair);
		  const span_pair_type  span_head(S, t, U, v);

		  value_type& value = chart(S, t, U, v);
		  
		  if (value.score == infinity)
		    agenda[span_head.size()].push_back(span_head);
		  
		  const double cost = (cost_str
				       + chart(span1.source.first, span1.source.last, span1.target.first, span1.target.last).score
				       + chart(span2.source.first, span2.source.last, span2.target.first, span2.target.last).score);
		  
		  if (cost < value.score) {
		    value  = (chart(span1.source.first, span1.source.last, span1.target.first, span1.target.last)
			      + chart(span2.source.first, span2.source.last, span2.target.first, span2.target.last));
		    
		    value.score     += cost_str;
		    value.inversion += cost_str;
		  }
		}
		
		// inversion
		for (difference_type U = v + (S == s); U <= utils::bithack::min(v + L, V); ++ U) {
		  // parent span: StuU
		  // span1: SsvU
		  // span2: stuv

		  if (chart(S, s, v, U).score == infinity) continue;
		  
		  const span_pair_type  span1(S, s, v, U);
		  const span_pair_type& span2(span_pair);
		  const span_pair_type  span_head(S, t, u, U);
		  
		  value_type& value = chart(S, t, u, U);
		  
		  if (value.score == infinity)
		    agenda[span_head.size()].push_back(span_head);
		  
		  const double cost = (cost_inv
				       + chart(span1.source.first, span1.source.last, span1.target.first, span1.target.last).score
				       + chart(span2.source.first, span2.source.last, span2.target.first, span2.target.last).score);
		  
		  if (cost < value.score) {
		    value  = (chart(span1.source.first, span1.source.last, span1.target.first, span1.target.last)
			      + chart(span2.source.first, span2.source.last, span2.target.first, span2.target.last));
		    
		    value.score     += cost_inv;
		    value.inversion += cost_inv;
		  }
		}
	      }
	      
	      for (difference_type S = t; S <= utils::bithack::min(t + l, T); ++ S) {
		const difference_type L = l - (S - t);
		
		// inversion
		for (difference_type U = utils::bithack::max(u - L, difference_type(0)); U <= u - (S == t); ++ U) {
		  // parent span: sSUv
		  // span1: stuv
		  // span2: tSUu
		  
		  if (chart(t, S, U, u).score == infinity) continue;
	      
		  const span_pair_type& span1(span_pair);
		  const span_pair_type  span2(t, S, U, u);
		  const span_pair_type  span_head(s, S, U, v);
		  
		  value_type& value = chart(s, S, U, v);
		  
		  if (value.score == infinity)
		    agenda[span_head.size()].push_back(span_head);
		  
		  const double cost = (cost_inv
				       + chart(span1.source.first, span1.source.last, span1.target.first, span1.target.last).score
				       + chart(span2.source.first, span2.source.last, span2.target.first, span2.target.last).score);
		  
		  if (cost < value.score) {
		    value  = (chart(span1.source.first, span1.source.last, span1.target.first, span1.target.last)
			      + chart(span2.source.first, span2.source.last, span2.target.first, span2.target.last));
		    
		    value.score     += cost_inv;
		    value.inversion += cost_inv;
		  }
		}
		
		// straight
		for (difference_type U = v + (S == t); U <= utils::bithack::min(v + L, V); ++ U) {
		  // parent span: sSuU
		  // span1: stuv
		  // span2: tSvU
		  
		  if (chart(t, S, v, U).score == infinity) continue;
		  
		  const span_pair_type& span1(span_pair);
		  const span_pair_type  span2(t, S, v, U);
		  const span_pair_type  span_head(s, S, u, U);
		  
		  value_type& value = chart(s, S, u, U);
		  
		  if (value.score == infinity)
		    agenda[span_head.size()].push_back(span_head);
		  
		  const double cost = (cost_str
				       + chart(span1.source.first, span1.source.last, span1.target.first, span1.target.last).score
				       + chart(span2.source.first, span2.source.last, span2.target.first, span2.target.last).score);
		  
		  if (cost < value.score) {
		    value  = (chart(span1.source.first, span1.source.last, span1.target.first, span1.target.last)
			      + chart(span2.source.first, span2.source.last, span2.target.first, span2.target.last));
		    
		    value.score     += cost_str;
		    value.inversion += cost_str;
		  }
		}
	      }
	    }
	  }
	  	  
	  beam *= 2;
	} while (chart(0, T, 0, V).score == infinity);

	return chart(0, T, 0, V);
      }
      
      void forward_backward(const chart_mono_type& chart, alpha_type& alpha, beta_type& beta)
      {
	const size_type sentence_size = chart.size() - 1;

	// forward...
	alpha[0] = 0;
	for (size_type last = 1; last <= sentence_size; ++ last)
	  for (size_type first = 0; first != last; ++ first)
	    alpha[last] = std::min(alpha[last], alpha[first] + chart(first, last));
	
	// backward...
	beta[sentence_size] = 0;
	for (difference_type first = sentence_size - 1; first >= 0; -- first)
	  for (size_type last = first + 1; last <= sentence_size; ++ last)
	    beta[first] = std::min(beta[first], chart(first, last) + beta[last]);
      }
      
    private:
      sentence_type ref;
    };
   
    InvWERScorer::InvWERScorer(const InvWERScorer& x)
      : Scorer(static_cast<const Scorer&>(*this)),
	weights(x.weights),
	matcher(0)
    {
      for (impl_set_type::const_iterator iter = x.impl.begin(); iter != x.impl.end(); ++ iter)
	impl.push_back(new impl_type(*(*iter)));
      
      if (x.matcher)
	matcher = &matcher_type::create(x.matcher->algorithm());
    }
    
    InvWERScorer::~InvWERScorer()
    {
      clear();
    }
    
    InvWERScorer& InvWERScorer::operator=(const InvWERScorer& x)
    {
      clear();
      
      static_cast<Scorer&>(*this) = static_cast<const Scorer&>(x);
      
      for (impl_set_type::const_iterator iter = x.impl.begin(); iter != x.impl.end(); ++ iter)
	impl.push_back(new impl_type(*(*iter)));
      
      weights = x.weights;
      matcher = 0;
      if (x.matcher)
	matcher = &matcher_type::create(x.matcher->algorithm());
      
      return *this;
    }
    
    void InvWERScorer::clear()
    {
      for (impl_set_type::iterator iter = impl.begin(); iter != impl.end(); ++ iter)
	delete *iter;
      
      impl.clear();
    }
    
    void InvWERScorer::insert(const sentence_type& __sentence)
    {
      sentence_type sentence;
      tokenize(__sentence, sentence);

      impl.push_back(new impl_type(sentence));
    }
    
    InvWERScorer::score_ptr_type InvWERScorer::score(const sentence_type& __sentence) const
    {
      sentence_type sentence;
      tokenize(__sentence, sentence);
      
      double score_best = std::numeric_limits<double>::infinity();

      std::auto_ptr<InvWER> inv_wer(new InvWER());
      
      for (impl_set_type::const_iterator iter = impl.begin(); iter != impl.end(); ++ iter) {
	impl_type& evaluator = const_cast<impl_type&>(*(*iter));
	
	impl_type::value_type value = evaluator(sentence, weights, matcher);
	
	if (value.score < score_best) {
	  score_best = value.score;
	  
	  inv_wer->insertion    = value.insertion;
	  inv_wer->deletion     = value.deletion;
	  inv_wer->substitution = value.substitution;
	  inv_wer->inversion    = value.inversion;
	}
	
	inv_wer->references += evaluator.ref.size();
      }
      
      if (! impl.empty())
	inv_wer->references /= impl.size();

      return score_ptr_type(inv_wer.release());
    }
  };
};

