// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//
#ifndef __CICADA_NGRAM_DEPENDENCY_EXTRACT_IMPL__HPP__
#define __CICADA_NGRAM_DEPENDENCY_EXTRACT_IMPL__HPP__

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>

#include <vector>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/sentence.hpp>
#include <cicada/dependency.hpp>

#include <boost/filesystem.hpp>

#include <utils/trie.hpp>
#include <utils/dense_hash_map.hpp>
#include <utils/lockfree_list_queue.hpp>
#include <utils/compress_stream.hpp>

#include "cicada_output_impl.hpp"

//
// input format:
// sentence1 ||| dependency1
// sentence2 ||| dependency2
// ...
//
// output format:
//
// counts/
//   unigram/
//     vocab.gz
//     vocab_cs.gz
//
//   root/
//     counts-xxx.gz
//     ...
//     root.idx
//   left/
//     counts-xxx.gz
//     ...
//     left.idx
//
//   left-head/
//     counts-xxx.gz
//     ...
//     left-head.idx
//
//   right/
//     counts-xxx.gz
//     ...
//     right.idx
//
//   right-head/
//     counts-xxx.gz
//     ...
//     right-head.idx
//

struct DependencyCounts
{
  typedef size_t             size_type;
  typedef ptrdiff_t          difference_type;
  
  typedef uint64_t           count_type;
  typedef cicada::Symbol     word_type;
  typedef cicada::Vocab      vocab_type;
  typedef cicada::Sentence   sentence_type;
  typedef cicada::Dependency dependency_type;
  
  typedef boost::filesystem::path                                     path_type;
  typedef std::vector<path_type, std::allocator<path_type> >          path_set_type;

  typedef utils::dense_hash_map<word_type, count_type, boost::hash<word_type>, std::equal_to<word_type>,
				std::allocator<std::pair<const word_type, count_type> > >::type root_count_set_type;
  
  typedef utils::trie<word_type, count_type, boost::hash<word_type>, std::equal_to<word_type>,
		      std::allocator<std::pair<const word_type, count_type> > > trie_count_set_type;
  
  struct count_set_type
  {
    root_count_set_type unigram;
    root_count_set_type root;
    trie_count_set_type left;
    trie_count_set_type right;
    
    count_set_type()
    {
      unigram.set_empty_key(word_type());
      root.set_empty_key(word_type());
    }
    
    void clear()
    {
      unigram.clear();
      root.clear();
      left.clear();
      right.clear();
    }

    void shrink()
    {
      root_count_set_type(unigram).swap(unigram);
      root_count_set_type(root).swap(root);
      trie_count_set_type(left).swap(left);
      trie_count_set_type(right).swap(right);
    }
    
    typedef std::vector<size_type, std::allocator<size_type> > index_set_type;
    typedef std::vector<index_set_type, std::allocator<index_set_type> > dependency_map_type;

    dependency_map_type dependency_map;
    
    void accumulate(const sentence_type& sentence,
		    const dependency_type& dependency)
    {
      if (sentence.size() != dependency.size())
	throw std::runtime_error("invalid sentence with dependency");
      
      if (sentence.empty()) return;
      
      dependency_map.clear();
      dependency_map.resize(sentence.size() + 1);
      
      for (size_type i = 0; i != dependency.size(); ++ i)
	dependency_map[dependency[i]].push_back(i + 1);
      
      if (dependency_map.front().empty())
	throw std::runtime_error("invalid dependency structure without root");

      // unigram...
      sentence_set_type::const_iterator siter_end = sentene.end();
      for (sentence_set_type::const_iterator siter = sentene.begin(); siter != siter_end; ++ siter)
	++ unigram[*siter];
      
      // root...
      index_set_type::const_iterator iiter_end = dependency_map[0].end();
      for (index_set_type::const_iterator iiter = dependency_map[0].begin(); iiter != iiter_end; ++ iiter)
	++ root[sentence[*iiter - 1]];
      
      // others
      for (size_type id = 1; id != dependency_map.size(); ++ id) 
	if (! dependency_map[id].empty()) {
	  index_set_type::const_iterator iiter_begin = dependency_map[id].begin();
	  index_set_type::const_iterator iiter_end   = dependency_map[id].end();
	  index_set_type::const_iterator iiter_lex   = std::lower_bound(iiter_begin, iiter_end, id);
	  
	  const word_type& head = sentence[id - 1];
	  
	  // left
	  if (iiter_begin != iiter_lex) {
	    // is it correct?
	    index_set_type::const_reverse_iterator liter_lex(iiter_lex);
	    index_set_type::const_reverse_iterator liter_end(iiter_begin);
	    
	    trie_count_set_type::id_type id_head = left_head.insert(left_head.root(), head);
	    
	    id_head = left_head.insert(id_head, sentence[*(liter_lex) - 1]);
	    ++ left_head[id_head];
	    
	    if (liter_lex + 1 != liter_begin) {
	      id_head = left_head.insert(id_head, sentence[*(liter_lex + 1) - 1]);
	      ++ left_head[id_head];
	    }
	    
	    // non-head
	    for (index_set_type::const_reverse_iterator liter = liter_lex + 2; liter <= liter_end; ++ liter) {
	      trie_count_set_type::id_type id = left.insert(left.root(), sentence[*(liter - 2) - 1]);
	      id = left.insert(id, sentence[*(liter - 1) - 1]);
	      id = left.insert(id, sentence[*(liter) - 1]);
	      ++ left[id];
	    }
	  }
	  
	  // right
	  if (iiter_lex != iiter_end) {
	    index_set_type::const_iterator riter_lex(iiter_lex);
	    index_set_type::const_iterator riter_end(iiter_end);
	    
	    // is it correct?
	    trie_count_set_tyep::id_type id_head = right_head.insert(right_head.root(), head);
	    
	    id_head = right_head.insert(id_head, sentene[*(riter_lex) - 1]);
	    ++ right_head[id_head];
	    
	    if (riter_lex + 1 != riter_end) {
	      id_head = right_head.insert(id_head, sentene[*(riter_lex + 1) - 1]);
	      ++ right_head[id_head];
	    }
	    
	    // non-head
	    for (index_set_type::const_iterator riter = riter_lex + 2; riter <= riter_end; ++ riter) {
	      trie_count_set_tyep::id_type id = right.insert(right.root(), sentence[*(riter - 2) - 1]);
	      id = right.insert(id, sentence[*(riter - 1) - 1]);
	      id = right.insert(id, sentence[*(riter) - 1]);
	      
	      ++ right[id];
	    }
	  }
	}
    }
  };

  struct count_path_set_type
  {
    path_set_type    unigram;
    path_set_type    root;
    path_set_type    left;
    path_set_type    right;
    path_set_type    left_head;
    path_set_type    right_head;

    count_path_set_type& operator+=(const count_path_set_type& x)
    {
      unigram.insert(unigram.end(), x.unigram.begin(), x.unigram.end());
      root.insert(root.end(), x.root.begin(), x.root.end());
      left.insert(left.end(), x.left.begin(), x.left.end());
      right.insert(right.end(), x.right.begin(), x.right.end());
      left_head.insert(left_head.end(), x.left_head.begin(), x.left_head.end());
      right_head.insert(right_head.end(), x.right_head.begin(), x.right_head.end());
    }
  };

  static
  bool parse(const utils::piece& line, sentence_type& sentence, dependency_type& dependency)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;

    sentence.clear();
    dependency.clear();
    
    std::string::const_iterator liter     = line.begin();
    std::string::const_iterator liter_end = line.end();
    
    if (! sentence.assign(liter, liter_end))
      return false;
    
    if (! qi::phrase_parse(liter, liter_end, "|||", standard::space))
      return false;
    
    if (! dependency.assign(liter, liter_end))
      return false;
    
    return liter == liter_end;
  }

  static
  path_type counts_file(const path_type& dir)
  {
    const path_type file_tmp = utils::tempfie::file_name(dir / "counts-XXXXXX");
    
    utils::tempfile::insert(file_tmp);
    
    return file_tmp.string() + ".gz";
  }
  
  template <typename Tp>
  struct less_pfirst_string
  {
    bool operator()(const Tp* x, const Tp* y) const
    {
      const std::string& wordx = static_cast<const std::string&>(x->first);
      const std::string& wordy = static_cast<const std::string&>(y->first);
      
      return wordx < wordy;
    }
  };

  template <typename Tp>
  struct greater_psecond
  {
    bool operator()(const Tp* x, const Tp* y) const
    {
      return x->second > y->second;
    }
  };

  template <typename Counts>
  static
  void dump(const path_type& path,
	    const Counts& counts)
  {
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;

    // sort, then dump...
    typedef typename Counts::value_type value_type;
    typedef std::vector<const value_type*, std::allocator<const value_type*> > sorted_type;

    sorted_type sorted(counts.size());
    
    {
      sorted_type::iterator siter = sorted.begin();
      typename Counts::const_iterator citer_end = counts.end();
      for (typename Counts::const_iterator citer = counts.begin(); citer != citer_end; ++ citer, ++ siter)
	*siter = &(*citer);
    }
    
    std::sort(sorted.begin(), sorted.end(), less_pfirst_string<value_type>());

    karma::uint_generator<count_type> count_generator;
    
    utils::compress_ostream os(path, 1024 * 1024);
    
    typename sorted_type::const_iterator siter_end = sorted.end();
    for (typename sorted_type::const_iterator siter = sorted.begin(); siter != siter_end; ++ siter)
      if (! karma::generate(std::ostream_iterator<char>(os),
			    standard::string << '\t' << count_generator << '\n',
			    (*siter)->first,
			    (*siter)->second))
	throw std::runtime_error("generation failed");
  }

  template <typename Prefix, typename Counts, typename Iterator, typename StreamIterator>
  void dump(const Prefix& prefix,
	    const Counts& counts,
	    Iterator first,
	    Iterator last,
	    StreamIterator stream_iter)
  {
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    typedef typename std::iterator_traits<Iterator>::value_type value_type;
    typedef std::vector<const value_type*, std::allocator<const value_type*> > sorted_type;
    
    sorted_type sorted;
    for (Iterator iter = first; iter != last; ++ iter)
      sorted.push_back(&(*iter));
    
    Prefix prefix_new(prefix.size() + 1);
    std::copy(prefix.begin(), prefix.end(), prefix_new.begin());

    karma::uint_generator<count_type> count_generator;
    
    typename sorted_type::const_iterator siter_end = sorted.end();
    for (typename sorted_type::const_iterator siter = sorted.begin(); siter != siter_end; ++ siter) {
      
      if (*stream_iter) {
	const count_type count = counts[(*siter)->second];
	
	if (count > 0)
	  if (! karma::generate(std::ostream_iterator<char>(*(*stream_iter)),
				*(standard::string << ' ') << standard::string << '\t' << count_generator << '\n',
				prefix,
				(*siter)->first,
				count))
	    throw std::runtime_error("generation failed");
      }
      
      if (! counts.empty((*siter)->second)) {
	prefix_new.back() = (*siter)->first;
	
	dump(prefix_new, counts, counts.begin((*siter)->second), counts.end((*siter)->second), stream_iter + 1);
      }
    }
  }
  
  template <typename Counts>
  static
  void dump(const path_type& path,
	    const path_type& path_head,
	    const Counts& counts)
  {
    typedef typename Counts::const_iterator const_iterator;
    typedef typename std::iterator_traits<const_iterator>::value_type value_type;
    typedef std::vector<const value_type*, std::allocator<const value_type*> > sorted_type;

    utils::compress_ostream os1(path, 1024 * 1024);
    utils::compress_ostream os2(path_head, 1024 * 1024);
    
    std::vector<std::ostream*, std::allocator<std::ostream*> > ostreams(3, 0);
    ostreams[1] = &os1;
    ostreams[2] = &os2;

    sentence_type prefix;
    
    dump(prefix, counts, counts.begin(), counts.end(), ostreams.begin());
  }
  
  template <typename Data>
  static
  void dump(const count_set_type& counts,
	    const path_type& path,
	    Data& data)
  {
    // data.path is a directory...
    
    const path_type dir_unigram       = path / "unigram";
    const path_type dir_root          = path / "root";
    const path_type dir_left          = path / "left";
    const path_type dir_right         = path / "right";
    const path_type dir_left_head  = path / "left-head";
    const path_type dir_right_head = path / "right-head";
    
    if (! boost::filesystem::exists(dir_unigram))
      throw std::runtime_error("no unigram directory?");
    if (! boost::filesystem::exists(dir_root))
      throw std::runtime_error("no root directory?");
    if (! boost::filesystem::exists(dir_left))
      throw std::runtime_error("no left directory?");
    if (! boost::filesystem::exists(dir_right))
      throw std::runtime_error("no right directory?");
    if (! boost::filesystem::exists(dir_left_head))
      throw std::runtime_error("no left-head directory?");
    if (! boost::filesystem::exists(dir_right_head))
      throw std::runtime_error("no right-head directory?");
    
    if (! counts.unigram.empty()) {
      data.unigram.push_back(counts_file(dir_unigram));
      
      dump(data.unigram.back(), counts.unigram);
    }
    
    if (! counts.root.empty()) {
      data.root.push_back(counts_file(dir_root));
      
      dump(data.root.back(), counts.root);
    }
    
    if (! counts.left.empty()) {
      data.left.push_back(counts_file(dir_left));
      data.left_head.push_back(counts_file(dir_left_head));
      
      dump(data.left.back(), data.left_head.back(), counts.left);
    }

    if (! counts.right.empty()) {
      data.right.push_back(counts_file(dir_right));
      data.right_head.push_back(counts_file(dir_right_head));
      
      dump(data.right.back(), data.right_head.back(), counts.right);
    }    
  }

  static
  void preprocess(const path_type& path)
  {
    prepare_directory(path);
    
    boost::filesystem::create_directory(path / "unigram");
    boost::filesystem::create_directory(path / "root");
    boost::filesystem::create_directory(path / "left");
    boost::filesystem::create_directory(path / "right");
    boost::filesystem::create_directory(path / "left-head");
    boost::filesystem::create_directory(path / "right-head");
  
    ::sync();
  }

  static
  void create_files(const path_type& path, const path_set_type& paths)
  {
    utils::compress_ostream os(path);
    path_set_type::const_iterator piter_end = paths.end();
    for (path_set_type::const_iterator piter = paths.begin(); piter != piter_end; ++ piter)
      os << piter->filename() << '\n';
  }

  static
  void postprocess(const path_type& path, const count_path_set_type& count_paths)
  {
    namespace qi = boost::spirit::qi;
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;

    //
    // merge unigram and vocab.gz and vocab_cs.gz
    //
    
    typedef utils::unordered_map<word_type,
				 count_type,
				 boost::hash<word_type>,
				 std::equal_to<word_type>,
				 std::allocator<std::pair<const word_type, count_type> > >::type word_set_type;

    typedef std::vector<const word_set_type::value_type*, std::allocator<const word_set_type::value_type*> > sorted_type;
    
    word_set_type words;
    
    qi::uint_parser<count_type> count_parser;
    
    typename path_set_type::const_iterator piter_end = count_paths.unigram.end();
    for (typename path_set_type::const_iterator piter = count_paths.unigram.begin(); piter != piter_end; ++ piter) {
      if (! boost::filesystem::exists(*piter))
	throw std::runtime_error(std::string("no unigramcounts? ") + piter->string());

      utils::tempfile::insert(*piter);
      
      typedef boost::spirit::istream_iterator iter_type;
      
      utils::compress_istream is(*piter, 1024 * 1024);
      is.unsetf(std::ios::skipws);

      iter_type iter(is);
      iter_type iter_end;

      std::string word;
      count_type  count;
      
      while (iter != iter_end) {
	word.clear();
	if (! qi::phrase_parse(iter, iter_end,
			       qi::lexeme[+(standard::char - standard::blank)] >> count_parser >> (qi::eol || qi::eoi),
			       standard::blank,
			       word,
			       count))
	  throw std::runtime_error("parsing failed");
	
	words[word] += count;
      }
    }

    sorted_type sorted(words.size());
    sorted_type::iterator siter = sorted.begin();
    word_set_type::const_iterator witer_end = words.end();
    for (word_set_type::const_iterator witer = words.begin(); witer != witer_end; ++ witer, ++ siter)
      *siter = &(*witer);

    karma::uint_generator<count_type> count_generator;
    
    // sort by word....
    {
      std::sort(sorted.begin(), sorted.end(), less_pfirst_string<word_set_type::value_type>());
      
      utils::compress_ostream os(path / "unigram" / "vocab.gz", 1024 * 1024);
      
      sorted_type::const_iterator siter_end = sorted.end();
      for (sorted_type::const_iterator siter = sorted.begin(); siter != siter_end; ++ siter)
	if (! karma::generate(std::ostream_iterator<char>(os),
			      standard::string << '\t' << count_generator << '\n',
			      (*siter)->first,
			      (*siter)->second))
	  throw std::runtime_error("generation failed");
    }
    
    // sort by count...
    {
      std::sort(sorted.begin(), sorted.end(), greater_psecond<word_set_type::value_type>());
      
      utils::compress_ostream os(path / "unigram" / "vocab_cs.gz", 1024 * 1024);
      
      sorted_type::const_iterator siter_end = sorted.end();
      for (sorted_type::const_iterator siter = sorted.begin(); siter != siter_end; ++ siter)
	if (! karma::generate(std::ostream_iterator<char>(os),
			      standard::string << '\t' << count_generator << '\n',
			      (*siter)->first,
			      (*siter)->second))
	  throw std::runtime_error("generation failed");
    }
    
    //
    // crate index files...
    //
    create_files(path / "root" / "root.idx",       counts_paths.root);
    create_files(path / "left" / "left.idx",       counts_paths.left);
    create_files(path / "right" / "right.idx",      counts_paths.right);
    create_files(path / "left-head" / "left-head.idx",  counts_paths.left_head);
    create_files(path / "right-head" / "right-lead.idx", counts_paths.right_head);
  }

  struct TaskFile
  {
    typedef utils::lockfree_list_queue<path_type, std::allocator<path_type> > queue_type;
    
    queue_type&      queue;
    path_type        path;
    double           max_malloc;
    
    count_path_set_type paths;
    
    TaskFile(queue_type&      _queue,
	     const path_type& _path,
	     const double     _max_malloc)
      : queue(_queue),
	path(_path),
	max_malloc(_max_malloc) {}
    
    void operator()()
    {
      count_set_type counts;
      path_type file;
      std::string line;
      
      sentence_type   sentence;
      dependency_type dependency;
      
      // every 1024 * 4 iterations, we will check for memory boundary
      const size_type iteration_mask = (1 << 12) - 1;
      size_type iteration = 0;
      while (1) {
	queue.pop_swap(file);
	if (file.empty()) break;
	
	if (file != "-" && ! boost::filesystem::exists(file))
	  throw std::runtime_error(std::string("no file? ") + file.string());
	
	utils::compress_istream is(file, 1024 * 1024);
	
	for (/**/; std::getline(is, line); ++ iteration) {
	  if (! DependencyCounts::parse(line, sentence, dependency))
	    throw std::runtime_error("invalid depdendency format");

	  counts.accumulate(sentence, dependency);
	  
	  if ((iteration & iteration_mask) == iteration_mask
	      && utils::malloc_stats::used() > size_t(max_malloc * 1024 * 1024 * 1024)) {
	    DependencyCounts::dump(counts, path, paths);
	    counts.clear();
	    counts.shrink();
	  }
	}
      }
      
      if (! counts.empty()) {
	DependencyCounts::dump(counts, path, paths);
	counts.clear();
      }
    }
  };
  
  struct TaskLine
  {
    typedef std::string line_type;
    typedef std::vector<line_type, std::allocator<line_type> > line_set_type;
    
    typedef utils::lockfree_list_queue<line_set_type, std::allocator<line_set_type> >  queue_type;
    
    queue_type&      queue;
    path_type        path;
    double           max_malloc;
    
    count_path_set_type paths;
    
    TaskLine(queue_type&      _queue,
	     const path_type& _path,
	     const double     _max_malloc)
      : queue(_queue),
	path(_path),
	max_malloc(_max_malloc) {}
    
    void operator()()
    {
      count_set_type counts;
      line_set_type  lines;
      
      sentence_type   sentence;
      dependency_type dependency;

      // every 1024 * 4 iterations, we will check for memory boundary
      const size_type iteration_mask = (1 << 12) - 1;
      size_type iteration = 0;
      while (1) {
	queue.pop_swap(lines);
	if (lines.empty()) break;
	
	line_set_type::const_iterator liter_end = lines.end();
	for (line_set_type::const_iterator liter = lines.begin(); liter != liter_end; ++ liter, ++ iteration) {
	  if (! DependencyCounts::parse(*liter, sentence, dependency))
	    throw std::runtime_error("invalid depdendency format");
	  
	  counts.accumulate(sentence, dependency);
	  
	  if ((iteration & iteration_mask) == iteration_mask
	      && utils::malloc_stats::used() > size_t(max_malloc * 1024 * 1024 * 1024)) {
	    DependencyCounts::dump(counts, path, paths);
	    counts.clear();
	    counts.shrink();
	  }
	}
      }
      
      if (! counts.empty()) {
	DependencyCounts::dump(counts, path, paths);
	counts.clear();
      }
    }
  };
};

#endif
