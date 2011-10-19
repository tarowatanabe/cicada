// MST related algorithms...

#ifndef __CICADA__OPTIMIZE_MST__HPP__
#define __CICADA__OPTIMIZE_MST__HPP__ 1

#include <algorithm>
#include <numeric>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>

#include <cicada/cicada/optimize/edmonds_optimum_branching.hpp>

#include <utils/mathop.hpp>

namespace cicada
{
  namespace optimize
  {
    struct MST
    {
      typedef double prob_type;
      typedef double    logprob_type;
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;

      typedef boost::numeric::ublas::matrix<prob_type> matrix_type;
      typedef boost::numeric::ublas::permutation_matrix<> permutation_matrix_type;

      typedef boost::property<boost::edge_weight_t, logprob_type>                               edge_property_type;
      typedef boost::adjacency_matrix<boost::directedS, boost::no_property, edge_property_type> graph_type;
      typedef boost::graph_traits<graph_type>::vertex_descriptor                                vertex_type;
      typedef boost::graph_traits<graph_type>::edge_descriptor                                  edge_type;
      
      typedef std::vector<vertex_type, std::allocator<vertex_type> > vertex_set_type;
      
      template <typename Dep>
      class insert_align_type
      {
	const graph_type* graph;
	Dep*              deps;
	logprob_type*     score;
      
      public:
	insert_align_type(const graph_type& _graph,
			  Dep& _deps,
			  logprob_type& _score)
	  : graph(&_graph),
	    deps(&_deps),
	    score(&_score) {}
      
	insert_align_type& operator=(const edge_type& edge)
	{	
	  deps->operator[](boost::target(edge, *graph) - 1) = boost::source(edge, *graph);
	
	  *score += boost::get(boost::get(boost::edge_weight_t(), *graph), edge);
	
	  return *this;
	}
      
	insert_align_type& operator*() { return *this; }
	insert_align_type& operator++() { return *this; }
	insert_align_type operator++(int) { return *this; }
      };
    
    
      template <typename Matrix, typename Vector>
      double viterbi_forest(Matrix& matrix, Vector& viterbi)
      {
	const size_type matrix_size = matrix.size1();
	const size_type vertex_size = matrix_size + 1;

	graph_type graph(vertex_size);
      
	vertices.clear();
	vertices.reserve(vertex_size);
      
	viterbi.clear();
	viterbi.resize(matrix_size);

      
	BOOST_FOREACH(vertex_type v, boost::vertices(graph))
	  vertices.push_back(v);

	const double logprob_zero = boost::numeric::bounds<typename Matrix::value_type>::lowest();
      
	for (int src = 0; src < matrix_size; ++ src)
	  for (int trg = 0; trg < matrix_size; ++ trg)  {
	    const int edge_source = (src != trg ? src + 1 : 0);
	    const int edge_target = trg + 1;

	    if (matrix(src, trg) == logprob_zero) continue;
	  
	    boost::add_edge(vertices[edge_source], vertices[edge_target], matrix(src, trg), graph);
	  }
      
	const boost::property_map<graph_type, boost::edge_weight_t>::type&  weights = boost::get(boost::edge_weight_t(), graph);
	const boost::property_map<graph_type, boost::vertex_index_t>::type& vertex_indices = boost::get(boost::vertex_index_t(), graph);

	double score = 0.0;
      
	edmonds_optimum_branching<true, true, true>(graph,
						    vertex_indices,
						    weights,
						    vertices.begin(),
						    vertices.begin() + 1,
						    insert_align_type<Vector>(graph, viterbi, score));
      
	return score;
      }
    
      template <typename Matrix, typename Vector>
      double viterbi_tree(Matrix& matrix, Vector& viterbi)
      {
	const size_type matrix_size = matrix.size1();
	const size_type vertex_size = matrix_size + 1;

	graph_type graph(vertex_size);
      
	vertices.clear();
	vertices.reserve(vertex_size);
      
	viterbi.clear();
	viterbi.resize(matrix_size);
      
	BOOST_FOREACH(vertex_type v, boost::vertices(graph))
	  vertices.push_back(v);

	const double logprob_zero = boost::numeric::bounds<typename Matrix::value_type>::lowest();

	logprob_type sum_weights = 1e5;
      
	for (int src = 0; src < matrix_size; ++ src)
	  for (int trg = 0; trg < matrix_size; ++ trg) 
	    if (src != trg) {
	      const int edge_source = src + 1;
	      const int edge_target = trg + 1;

	      if (matrix(src, trg) == logprob_zero) continue;
	    
	      boost::add_edge(vertices[edge_source], vertices[edge_target], matrix(src, trg), graph);
	    
	      sum_weights += matrix(src, trg);
	    }
      
	// diagonal... first compute sum...
	for (int i = 0; i < matrix_size; ++ i)
	  sum_weights += std::fabs(matrix(i, i));
      
	// then, subtract...
	for (int i = 0; i < matrix_size; ++ i)
	  boost::add_edge(vertices[0], vertices[i + 1], matrix(i, i) - sum_weights, graph);
      
	const boost::property_map<graph_type, boost::edge_weight_t>::type&  weights = boost::get(boost::edge_weight_t(), graph);
	const boost::property_map<graph_type, boost::vertex_index_t>::type& vertex_indices = boost::get(boost::vertex_index_t(), graph);

	double score = sum_weights;
      
	edmonds_optimum_branching<true, true, true>(graph,
						    vertex_indices,
						    weights,
						    vertices.begin(),
						    vertices.begin() + 1,
						    insert_align_type<Vector>(graph, viterbi, score));
      
	return score;
      }

    
      template <typename Matrix>
      double marginal_forest(Matrix& matrix, Matrix& marginal)
      {
	const size_type matrix_size = matrix.size1();
      
	double logsum_adjust = 0.0;
	for (int j = 0; j < matrix_size; ++ j) {
	  boost::numeric::ublas::matrix_column<matrix_type> column(matrix, j);
	
	  const double sum = std::accumulate(column.begin(), column.end(), 0.0);
	
	  std::transform(column.begin(), column.end(), column.begin(), std::bind2nd(std::multiplies<prob_type>(), 1.0 / sum));
	
	  logsum_adjust += utils::mathop::log(sum);
	}
      
	// create transposed laplacian... (dL)
      
	laplacian.resize(matrix_size, matrix_size);
	for (int j = 0; j < matrix_size; ++ j) {
	  double sum = 0.0;
	  for (int i = 0; i < matrix_size; ++ i) {
	    sum += matrix(i, j);
	    laplacian(j, i) = - matrix(i, j);
	  }
	  laplacian(j, j) = sum;
	}
      
	// compute inverse...
	const double logsum = inverse(laplacian, laplacian_inverted) + logsum_adjust;
      
	// compute again...
	marginal.resize(matrix_size, matrix_size);
	for (int j = 0; j < matrix_size; ++ j) {
	  const double laplacian_inverted_jj = laplacian_inverted(j, j);
	
	  double sum = 0.0;
	  for (int i = 0; i < matrix_size; ++ i) {
	    marginal(i, j) = std::max(matrix(i, j) * (laplacian_inverted_jj - laplacian_inverted(i, j) * (i != j)), boost::numeric::bounds<double>::smallest());
	    sum += marginal(i, j);
	  }
	
	  boost::numeric::ublas::matrix_column<Matrix> column(marginal, j);
	  std::transform(column.begin(), column.end(), column.begin(), std::bind2nd(std::multiplies<prob_type>(), 1.0 / sum));
	}
      
	return logsum;
      }
    
      template <typename Matrix>
      double marginal_tree(Matrix& matrix, Matrix& marginal)
      {
	const size_type matrix_size = matrix.size1();

	double logsum_adjust = 0.0;
	for (int j = 0; j < matrix_size; ++ j) {
	  boost::numeric::ublas::matrix_column<Matrix> column(matrix, j);
	
	  const double sum = std::accumulate(column.begin(), column.end(), 0.0);
	
	  std::transform(column.begin(), column.end(), column.begin(), std::bind2nd(std::multiplies<prob_type>(), 1.0 / sum));
	
	  logsum_adjust += utils::mathop::log(sum);
	}
      
	// create transposed laplacian... (dL)

	laplacian.resize(matrix_size, matrix_size);
	laplacian.clear();
      
	// recover root scores from diagonal
	for (int i = 0; i < matrix_size; ++ i)
	  laplacian(i, 0) = matrix(i, i);
      
	// a-bit confusing assignment
	for (int j = 0; j < matrix_size; ++ j) {
	  double sum = 0.0;
	
	  for (int i = 0; i < matrix_size; ++ i) 
	    if (i != j) { // ignore diagonal...
	      sum += matrix(i, j);
	    
	      // laplacian(j, 0) is simply copied from matrix...
	      if (i > 0)
		laplacian(j, i) = - matrix(i, j);
	    }
	
	  if (j > 0)
	    laplacian(j, j) = sum;
	}
      
	const double logsum = inverse(laplacian, laplacian_inverted) + logsum_adjust;
      
	marginal.resize(matrix_size, matrix_size);
	marginal.clear();

	// marginals for roots
	marginal(0, 0) = std::max(matrix(0, 0) * laplacian_inverted(0, 0), boost::numeric::bounds<double>::smallest());
      
	double sum = marginal(0, 0);
	for (int i = 1; i < matrix_size; ++ i) {
	  marginal(i, i) = std::max(  matrix(i, i) * laplacian_inverted(0, i), boost::numeric::bounds<double>::smallest());
	  marginal(0, i) = std::max(  matrix(0, i) * laplacian_inverted(i, i), boost::numeric::bounds<double>::smallest());
	  marginal(i, 0) = std::max(- matrix(i, 0) * laplacian_inverted(i, 0), boost::numeric::bounds<double>::smallest());
	
	  sum += marginal(i, 0);
	}
      
	// re-normalize...
	{
	  boost::numeric::ublas::matrix_column<Matrix> column(marginal, 0);
	  std::transform(column.begin(), column.end(), column.begin(), std::bind2nd(std::multiplies<prob_type>(), 1.0 / sum));
	}
	//for (int i = 0; i < matrix_size; ++ i)
	//  marginal(i, 0) /= sum;
      
	for (int j = 1; j < matrix_size; ++ j) {
	
	  const double laplacian_inverted_jj = laplacian_inverted(j, j);
	
	  double sum = marginal(j, j) + marginal(0, j);
	  for (int i = 1; i < matrix_size; ++ i) 
	    if (i != j) {
	      marginal(i, j) = std::max(matrix(i, j) * (laplacian_inverted_jj - laplacian_inverted(i, j)), boost::numeric::bounds<double>::smallest());
	      sum += marginal(i, j);
	    }
	
	  boost::numeric::ublas::matrix_column<Matrix> column(marginal, j);
	  std::transform(column.begin(), column.end(), column.begin(), std::bind2nd(std::multiplies<prob_type>(), 1.0 / sum));
	
	  //for (int i = 0; i < matrix_size; ++ i)
	  //  marginal(i, j) /= sum;
	}
      
	return logsum;
      }

      template <typename Matrix>
      double inverse(Matrix& matrix, Matrix& inverted)
      {
	permutation_matrix_type pm(matrix.size1());
      
	boost::numeric::ublas::lu_factorize(matrix, pm);
      
	double logprod = 0.0;
	for (int i = 0; i < matrix.size1(); ++ i)
	  logprod += utils::mathop::log(fabs(matrix(i, i))); // log of absolute value...
      
	inverted.resize(matrix.size1(), matrix.size2());
	inverted = boost::numeric::ublas::identity_matrix<prob_type>(matrix.size1());
      
	boost::numeric::ublas::lu_substitute(matrix, pm, inverted);
      
	return logprod;
      }
    
      matrix_type laplacian;
      matrix_type laplacian_inverted;

      vertex_set_type vertices;
    };
  };
};


#endif
