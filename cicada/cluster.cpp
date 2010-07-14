
#include <cicada/cluster.hpp>

#include <utils/tempfile.hpp>
#include <utils/compress_stream.hpp>
#include <utils/space_separator.hpp>

#include <boost/tokenizer.hpp>

namespace cicada
{
  
  void Cluster::open(const path_type& path)
  {
    typedef boost::tokenizer<utils::space_separator> tokenizer_type;
    typedef utils::repository repository_type;
    
    if (boost::filesystem::is_directory(path)) {
      repository_type rep(path, repository_type::read);
      
      vocab.open(rep.path("vocab"));
      clusters.open(rep.path("clusters"));
    } else {
      typedef std::vector<id_type, std::allocator<id_type> > cluster_map_type;

      cluster_map_type cluster_map;
      
      utils::compress_istream is(path, 1024 * 1024);
      
      std::string line;
      
      while (std::getline(is, line)) {
	tokenizer_type tokenizer(line);
	
	tokenizer_type::iterator iter = tokenizer.begin();
	if (iter == tokenizer.end()) continue;
	
	const std::string __cluster = *iter;
	++ iter;
	if (iter == tokenizer.end()) continue;
	
	const std::string __word = *iter;
	
	const word_type cluster(__cluster);
	const word_type word(__word);
	
	if (word.id() >= cluster_map.size())
	  cluster_map.resize(word.id() + 1, 0);
	
	cluster_map[word.id()] = cluster.id() + 1;
      }
      
      const path_type tmp_dir = utils::tempfile::tmp_dir();
      const path_type path_repository = utils::tempfile::directory_name(tmp_dir / "cicada.cluster.XXXXXX");
      
      utils::tempfile::insert(path_repository);
      
      repository_type rep(path_repository, repository_type::write);
      
      boost::iostreams::filtering_ostream os;
      os.push(utils::packed_sink<id_type, std::allocator<id_type> >(rep.path("clusters")), 1024 * 1024);
      
      cluster_map_type::const_iterator iter_end = cluster_map.end();
      for (cluster_map_type::const_iterator iter = cluster_map.begin(); iter != iter_end; ++ iter)
	os.write((char*) &(*iter), sizeof(id_type));
      
      os.pop();
      
      word_type::write(rep.path("vocab"));
      
      vocab.open(rep.path("vocab"));
      clusters.open(rep.path("clusters"));
    }
  }

  void Cluster::write(const path_type& path) const
  {
    typedef utils::repository repository_type;

    if (empty()) return;
    
    const path_type path_curr = clusters.path().parent_path();

    // do not write at the same file
    if (path_curr == path) return;
    
    repository_type rep(path, repository_type::write);
    
    vocab.write(rep.path("vocab"));
    clusters.write(rep.path("clusters"));
  }
};
