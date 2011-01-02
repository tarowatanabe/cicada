// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__SEMIRING__ENVELOPE__HPP__
#define __CICADA__SEMIRING__ENVELOPE__HPP__ 1

// @InProceedings{kumar-EtAl:2009:ACLIJCNLP,
//   author    = {Kumar, Shankar  and  Macherey, Wolfgang  and  Dyer, Chris  and  Och, Franz},
//   title     = {Efficient Minimum Error Rate Training and Minimum Bayes-Risk Decoding for Translation Hypergraphs and Lattices},
//   booktitle = {Proceedings of the Joint Conference of the 47th Annual Meeting of the ACL and the 4th International Joint Conference on Natural Language Processing of the AFNLP},
//   month     = {August},
//   year      = {2009},
//   address   = {Suntec, Singapore},
//   publisher = {Association for Computational Linguistics},
//   pages     = {163--171},
//   url       = {http://www.aclweb.org/anthology/P/P09/P09-1019}
//   }

#include <limits>

#include <cicada/semiring/traits.hpp>

#include <cicada/hypergraph.hpp>
#include <cicada/sentence.hpp>
#include <cicada/feature_vector.hpp>

#include <utils/simple_vector.hpp>

#include <boost/shared_ptr.hpp>

namespace cicada
{
  namespace semiring
  {
    
    class Envelope
    {
    public:
      typedef cicada::HyperGraph hypergraph_type;
      typedef cicada::Sentence sentence_type;
      

    public:
      struct Line
      {
	typedef Line line_type;
	typedef boost::shared_ptr<line_type>  line_ptr_type;
	
	Line()
	  : x(0.0), m(0.0), y(0.0), edge(0) {}
	Line(const double& __m, const double& __y, const hypergraph_type::edge_type& __edge)
	  : x(- std::numeric_limits<double>::infinity()), m(__m), y(__y), edge(&__edge) {}
	Line(const double& __x, const double& __m, const double& __y)
	  : x(__x), m(__m), y(__y), edge(0) {}
	Line(const double& __x, const double& __m, const double& __y, const line_ptr_type& __parent, const line_ptr_type& __antecedent)
	  : x(__x), m(__m), y(__y), edge(0), parent(__parent), antecedent(__antecedent) {}
	
	double x;
	double m;
	double y;
	
	// D part...
	const hypergraph_type::edge_type* edge;
	
	line_ptr_type parent;
	line_ptr_type antecedent;
	
	void yield(sentence_type& sentence) const;
      };
      
      typedef Line line_type;
      typedef boost::shared_ptr<line_type>  line_ptr_type;
      typedef std::vector<line_ptr_type, std::allocator<line_ptr_type> > line_ptr_set_type;

      typedef line_ptr_set_type::const_iterator const_iterator;

    public:      
      Envelope() : lines(), is_sorted(true) {}
      Envelope(const line_ptr_type& line) : lines(1, line), is_sorted(true) {}
      
    public:
      const Envelope& operator+=(const Envelope& x);
      const Envelope& operator*=(const Envelope& x);
      
    public:
      const_iterator begin() const { return lines.begin(); }
      const_iterator end()   const { return lines.end(); }
      
      void sort();
      
    private:
      line_ptr_set_type lines;
      bool is_sorted;
    };

    template <typename FeatureVector>
    struct EnvelopeFunction
    {
      typedef FeatureVector feature_set_type;
      
      typedef Envelope envelope_type;
      typedef envelope_type::line_type line_type;
      
      EnvelopeFunction(const feature_set_type& __origin,
		       const feature_set_type& __direction)
	: origin(__origin), direction(__direction) {}
      
      template <typename Edge>
      Envelope operator()(const Edge& edge) const
      {
	const double m = edge.features.dot(direction);
	const double y = edge.features.dot(origin);
	
	return Envelope(boost::shared_ptr<line_type>(new line_type(m, y, edge)));
      }
      
      const feature_set_type& origin;
      const feature_set_type& direction;
    };
    
    template <>
    struct traits<Envelope>
    {
      static inline Envelope zero() { return Envelope();  }
      static inline Envelope one()  { return Envelope(boost::shared_ptr<Envelope::line_type>(new Envelope::line_type())); }
    };

  };
};

#endif
