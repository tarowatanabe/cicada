

#include <iostream>
#include <vector>
#include <utility>
#include <stdexcept>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/fusion/adapted.hpp>
#include <boost/filesystem.hpp>

#include "utils/compress_stream.hpp"
#include "utils/lexical_cast.hpp"

#include "weight_vector.hpp"
#include "feature_vector.hpp"
#include "optimize_qp.hpp"

typedef cicada::WeightVector<double> weight_set_type;
typedef std::pair<int, double> feature_type;
typedef std::vector<feature_type, std::allocator<feature_type> > feature_set_type;
struct FeatureVector
{
  double label;
  feature_set_type features;
  
  FeatureVector(const double& __label, const feature_set_type& __features) : label(__label), features(__features) {} 
  FeatureVector() : label(), features() {}
};
typedef FeatureVector feature_vector_type;
typedef std::vector<feature_vector_type, std::allocator<feature_vector_type> > model_type;
typedef std::vector<double, std::allocator<double> > alpha_type;
typedef std::vector<double, std::allocator<double> > f_type;
typedef boost::filesystem::path path_type;

BOOST_FUSION_ADAPT_STRUCT(
			  feature_vector_type,
			  (double,           label)
			  (feature_set_type, features)
			  )

template <typename Iterator>
struct feature_vector_parser : boost::spirit::qi::grammar<Iterator, model_type(), boost::spirit::standard::blank_type>
{
  feature_vector_parser() : feature_vector_parser::base_type(model)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    features %= qi::double_ >> *(qi::int_ >> ':' >> qi::double_) >> (qi::eol | qi::eoi);
    model %= *features;
  }
  
  typedef boost::spirit::standard::blank_type blank_type;
  
  boost::spirit::qi::rule<Iterator, feature_vector_type(), blank_type> features;
  boost::spirit::qi::rule<Iterator, model_type(),          blank_type> model;
};

struct H
{
  H(const model_type& __model) : model(__model) {}
  
  double operator()(int i, int j) const
  {
    feature_set_type::const_iterator iter1     = model[i].features.begin();
    feature_set_type::const_iterator iter1_end = model[i].features.end();
    feature_set_type::const_iterator iter2     = model[j].features.begin();
    feature_set_type::const_iterator iter2_end = model[j].features.end();
    
    double dot = 0.0;
    
    const double factor1 = (model[i].label >= 0.0 ? 1.0 : - 1.0);
    const double factor2 = (model[j].label >= 0.0 ? 1.0 : - 1.0);
    const double factor = factor1 * factor2;
    
    while (iter1 != iter1_end && iter2 != iter2_end) {
      if (iter1->first < iter2->first)
	++ iter1;
      else if (iter2->first < iter1->first)
	++ iter2;
      else {
	dot += iter1->second * iter2->second * factor;
	++ iter1;
	++ iter2;
      }
    }
    
    return dot;
  }
  
  const model_type& model;
};

struct M
{

  M(const model_type& __model) : model(__model) {}

  template <typename W>
  void operator()(W& w, const alpha_type& alpha) const
  {
    for (size_t i = 0; i != model.size(); ++ i) {
      const double factor = (model[i].label >= 0.0 ? 1.0 : - 1.0);
      
      feature_set_type::const_iterator fiter_end = model[i].features.end();
      for (feature_set_type::const_iterator fiter = model[i].features.begin(); fiter != fiter_end; ++ fiter)
	w[fiter->first] += alpha[i] * fiter->second * factor;
    }
  }
  
  template <typename W>
  double operator()(const W& w, const size_t& i) const
  {
    const double factor = (model[i].label >= 0.0 ? 1.0 : - 1.0);
    
    double value = 0.0;
    feature_set_type::const_iterator fiter_end = model[i].features.end();
    for (feature_set_type::const_iterator fiter = model[i].features.begin(); fiter != fiter_end; ++ fiter)
      value += w[fiter->first] * fiter->second * factor;
      
    return value;
  }
  
  template <typename W>
  void operator()(W& w, const double& update, const size_t& i) const
  {
    const double factor = (model[i].label >= 0.0 ? 1.0 : - 1.0);

    feature_set_type::const_iterator fiter_end = model[i].features.end();
    for (feature_set_type::const_iterator fiter = model[i].features.begin(); fiter != fiter_end; ++ fiter)
      w[fiter->first] += update * fiter->second * factor;
  }
  
  const model_type& model;
};

int main(int argc, char** argv)
{
  // argv[1] for algorithm: 0 for dcd, 1 for smo, 2 for simplex
  // argv[2] for the feature set
  // argv[3] for output by alpha
  // argv[4] for output by weights
  
  if (argc != 5) {
    std::cerr << argv[0] << " [alg] [feature file] [outut alpha] [output weights]" << std::endl;
    return 1;
  }
  
  try {
    const int alg = utils::lexical_cast<int>(argv[1]);
    const path_type input_path(argv[2]);
    const path_type output_alpha_path(argv[3]);
    const path_type output_weights_path(argv[4]);

    typedef boost::spirit::istream_iterator iter_type;
    typedef feature_vector_parser<iter_type> parser_type;
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    model_type  model;
    parser_type parser;

    utils::compress_istream is(input_path, 1024 * 1024);
    is.unsetf(std::ios::skipws);
    
    iter_type iter(is);
    iter_type iter_end;
      
    if (! qi::phrase_parse(iter, iter_end, parser, standard::blank, model))
      throw std::runtime_error("failed model parsing");
    
    std::cerr << "# of features: " << model.size() << std::endl;
    
    // initialize alpha
    alpha_type alpha(model.size(), 0.0);
    f_type     f(model.size());
    for (size_t i = 0; i != model.size(); ++ i)
      f[i] = - std::fabs(model[i].label);
    
    H h(model);
    M m(model);
    
    double objective = 0.0;
    if (alg == 0) {
      cicada::optimize::QPSMO solver;
      objective = solver(alpha, f, h, m, 1.0 / (1e-4 *  model.size()), 1e-4);
    } else if (alg == 1) {
      cicada::optimize::QPDCD solver;
      objective = solver(alpha, f, h, m, 1.0 / (1e-4 *  model.size()), 1e-4);
    } else if (alg == 2) {
      cicada::optimize::QPSimplex solver;
      objective = solver(alpha, f, h, m, 1.0 / (1e-4 *  model.size()), 1e-4);
    } else
      throw std::runtime_error("algorithm can be 0 (dcd), 1(smo) or 2(simplex)");
    
    {
      utils::compress_ostream os(output_alpha_path);
      for (size_t i = 0; i != alpha.size(); ++ i)
	os << alpha[i] << '\n';
    }

    weight_set_type weights;
    
    m(weights, alpha);
    
    {
      utils::compress_ostream os(output_weights_path);
      for (size_t i = 0; i != weights.size(); ++ i)
	os << weights[i] << '\n';
    }
    
    size_t correct = 0;
    for (size_t i = 0; i != model.size(); ++ i) {
      const double predicted = m(weights, i);
      
      correct += (predicted * model[i].label) > 0.0;
    }

    std::cerr << "correct: " << correct << std::endl;
    
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}
