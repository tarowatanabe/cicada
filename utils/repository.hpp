// -*- mode: c++ -*-

#ifndef __UTILS__REPOSITORY__HPP__
#define __UTILS__REPOSITORY__HPP__ 1



// repository handling.
// each directory must hold
//   prop.list
// which describes key->value mapping (string->string mapping...)
// with additional directories/files (user defined...)
//


//
// repository will keep three things:
//
// prop.list file
//    <key> = <value>
//   pairs, but no spaces allowed for key/value
//
// other files
// other directories (possiblly sub-repositories...)
//
// it is user's responsibility how to create files/directories...
//

//
// close() will simply dump the prop.list again...
//
// open() with write_mode will create directory, and clear the content
// open() with read_mode will open, and read the contents 
//


//
// files (or directory), iterator is supported
// find_file(), begin_file(), end_file()
// 
// directory() return the directory of the repository
// 
// create_directory(name) // yes, we need this...???
//   NOTE: we can create sub-repository by
//
//    repository.open(a-repository)
//    repository_sub.open(repository.directory() / sub_repository)
//
// all the names are relative to the repository...
// 
//
// for prop.list, map-like operation is supported
// find(), begin(), end(), insert(), erase(), operator[]
// the content is dumpped after close() "when modified" and/or when no prop.list exsists
//

#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <stdexcept>

#include <boost/shared_ptr.hpp>
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>

#include <utils/filesystem.hpp>
#include <utils/space_separator.hpp>

namespace utils
{
  
  //
  // support for read/write...
  // currently, read is supported...
  // how to set up write operations???
  //
  
  class repository
  {
  private:
    typedef std::map<std::string, std::string> props_type;
    
  public:
    typedef boost::filesystem::path path_type;
    typedef boost::filesystem::directory_iterator directory_iterator;
    
  public:
    typedef props_type::const_iterator  const_iterator;
    typedef props_type::iterator        iterator;
    typedef props_type::const_reference const_reference;
    typedef props_type::reference       reference;
    
  public:
    
    // this is exclusive... no read/write simultaneously
    typedef enum { 
      read,
      write,
    } mode_type;
    
  private:
    path_type repository_dir;
    props_type props;
    mode_type repository_mode;
    bool modified;
    
  public:
    repository() : repository_mode(read), modified(false) {}
    repository(const path_type& dir, const mode_type mode=read): repository_mode(read), modified(false)  { open(dir, mode); }
    repository(const repository& x) : repository_mode(read), modified(false)
    {
      if (x.repository_mode == write) 
	throw std::runtime_error("cannot initialize with write mode");
      if (! x.repository_dir.empty() && boost::filesystem::exists(x.repository_dir) && boost::filesystem::is_directory(x.repository_dir))
	open_read(x.repository_dir);
    }
    ~repository() throw() { close(); }
    
    repository& operator=(const repository& x)
    {
      close();
      if (x.repository_mode == write) 
	throw std::runtime_error("cannot copy with write mode");
      if (! x.repository_dir.empty() && boost::filesystem::exists(x.repository_dir) && boost::filesystem::is_directory(x.repository_dir))
	open_read(x.repository_dir);
      return *this;
    }
    
  public:
    
    // repository operations
    void open(const path_type& dir, const mode_type mode=read)
    {
      // close first
      close();
      
      // then, open depending on its mode
      if (mode == read)
	open_read(dir);
      else
	open_write(dir);
    }
    
    
    // clear the contents
    void clear() {
      // remove props,
      if (repository_mode == read)
	throw std::runtime_error("clearing read mode repository");
      if (repository_dir.empty())
	throw std::runtime_error("no repository");
      
      // clear properties
      props.clear();
      
      // if not directory, remove
      if (boost::filesystem::exists(repository_dir) && ! boost::filesystem::is_directory(repository_dir))
	boost::filesystem::remove_all(repository_dir);
      
      // create directory
      boost::filesystem::create_directories(repository_dir);
      
      // remove all the files under the repository_dir
      boost::filesystem::directory_iterator iter_end;
      for (boost::filesystem::directory_iterator iter(repository_dir); iter != iter_end; ++ iter)
	boost::filesystem::remove_all(*iter);
    }
    
    // close the repository... this will dump the contents of props at props.list if write-mode
    void close()
    {
      
      if ((repository_mode == write || modified) && ! repository_dir.empty() && boost::filesystem::exists(repository_dir) && boost::filesystem::is_directory(repository_dir)) {
	
	// here, dump contents of props at "prop.list"
	std::ofstream os((repository_dir / "prop.list").string().c_str());
	props_type::const_iterator piter_end = props.end();
	for (props_type::const_iterator piter = props.begin(); piter != piter_end; ++ piter)
	  os << piter->first << " = " << piter->second << "\n";
      }
      
      // clear props and repository_dir
      props.clear();
      repository_dir = "";
      modified = false;
    }
    
    // repository access
    std::string& operator[](const std::string& key) 
    { 
      if (repository_mode == read)
	modified = true;

      return props[key];
    }
    
    const_iterator find(const std::string& key) const { return props.find(key); }
    iterator find(const std::string& key) { return props.find(key); }
    
    const_iterator begin() const { return props.begin(); }
    iterator begin() { return props.begin(); }
    const_iterator end() const { return props.end(); }
    iterator end() { return props.end(); }
    
    // for directory access...
    const path_type& base() const { return repository_dir; }

    path_type path(const std::string& file) const
    {
      if (file == "prop.list") throw std::runtime_error("not allowed file/directory");
      return repository_dir / file;
    }
    
    directory_iterator begin_file() const { return directory_iterator(repository_dir); }
    directory_iterator end_file() const { return directory_iterator(); }
    
    bool is_reader() const { return repository_mode == read; }
    bool is_writer() const { return repository_mode == write; }
    
  private:
    
    void open_read(const path_type& path)
    {
      typedef boost::tokenizer<utils::space_separator> tokenizer_type;

      modified = false;
      
      if (path.empty())
	throw std::runtime_error("empty path");
      if (! boost::filesystem::exists(path) || ! boost::filesystem::is_directory(path))
	throw std::runtime_error("invalid directory");
      
      repository_dir = path;
      repository_mode = read;
      
      const path_type prop_path = repository_dir / "prop.list";
      
      if (! boost::filesystem::exists(prop_path) || boost::filesystem::is_directory(prop_path))
	throw std::runtime_error("invalid property list");
      
      std::ifstream is(prop_path.string().c_str());
      std::string line;
      
      while (std::getline(is, line)) {
	
	tokenizer_type tokenizer(line);
	tokenizer_type::iterator iter = tokenizer.begin();
	if (iter == tokenizer.end()) continue;
	
	const std::string key = *iter;
	++ iter;
	
	if (iter == tokenizer.end() || *iter != "=") continue;
	++ iter;
	
	if (iter == tokenizer.end()) continue;
	const std::string value = *iter;
	
	//++ iter;
	//if (iter != tokenizer.end()) continue;
	
	props[key] = value;
      }
    }
    
    void open_write(const path_type& path)
    {
      if (path.empty())
	throw std::runtime_error("no directory");
      
      repository_dir = path;
      repository_mode = write;
      modified = false;
      
      clear();
    }    
  };
};

#endif
