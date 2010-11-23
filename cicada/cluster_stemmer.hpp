// -*- mode: c++ -*-

#ifndef __CICADA__CLUSTER_STEMMER__HPP__
#define __CICADA__CLUSTER_STEMMER__HPP__ 1

#include <stdexcept>

#include <cicada/cluster.hpp>
#include <cicada/stemmer.hpp>

namespace cicada
{
  
  class ClusterStemmer
  {
  private:
    typedef Cluster cluster_type;
    typedef Stemmer stemmer_type;
    typedef Symbol  symbol_type;

  public:
    ClusterStemmer(const cluster_type* __cluster)
      : cluster(__cluster), stemmer(0) {}
    ClusterStemmer(const stemmer_type* __stemmer)
      : cluster(0), stemmer(__stemmer) {}
    ClusterStemmer(const ClusterStemmer& x)
      : cluster(0), stemmer(0) 
    {
      if (x.cluster)
	cluster = &cluster_type::create(x.cluster->path());
      else if (x.stemmer)
	stemmer = &stemmer_type::create(x.stemmer->algorithm());
    }

    ClusterStemmer& operator=(const ClusterStemmer& x)
    {
      if (x.cluster)
	cluster = &cluster_type::create(x.cluster->path());
      else if (x.stemmer)
	stemmer = &stemmer_type::create(x.stemmer->algorithm());
      return *this;
    }
    
  public:
    symbol_type operator()(const symbol_type& word) const { return operator[](word); }
    symbol_type operator[](const symbol_type& word) const
    {
      if (cluster)
	return cluster->operator()(word);
      else if (stemmer)
	return stemmer->operator()(word);
      else
	throw std::runtime_error("no cluster or stemmer");
    }

  private:
    const cluster_type* cluster;
    const stemmer_type* stemmer;
  };
};


#endif
