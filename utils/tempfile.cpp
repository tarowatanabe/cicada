//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <unistd.h>
#include <signal.h>

#include <boost/thread.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string/replace.hpp>

#include <utils/filesystem.hpp>

#include "tempfile.hpp"

namespace utils
{
  namespace tempfile_impl
  {
    typedef boost::mutex              mutex_type;
    typedef boost::mutex::scoped_lock lock_type;
    
    static mutex_type __mutex;
    
    // file management...
    struct SignalBlocker
    {
      SignalBlocker()
      {
	sigemptyset(&mask);
	
	sigaddset(&mask, SIGHUP);
	sigaddset(&mask, SIGINT);
	sigaddset(&mask, SIGQUIT);
	sigaddset(&mask, SIGILL);
	sigaddset(&mask, SIGABRT);
	sigaddset(&mask, SIGKILL);
	sigaddset(&mask, SIGSEGV);
	sigaddset(&mask, SIGTERM);
	sigaddset(&mask, SIGBUS);
	
	sigprocmask(SIG_BLOCK, &mask, &mask_saved);
      }
      
      ~SignalBlocker()
      {
	sigprocmask(SIG_SETMASK, &mask_saved, 0);
      }

      sigset_t mask;
      sigset_t mask_saved;
    };

    struct PathSet
    {
      typedef boost::filesystem::path                            path_type;
      typedef std::vector<path_type, std::allocator<path_type> > path_set_type;
      
      PathSet() {}
      ~PathSet() { clear(); }
      
      void insert(const path_type& path)
      {
	if (path.empty()) return;
	
	SignalBlocker __blocker;
	
	lock_type lock(__mutex);
	
	path_set_type::iterator piter = std::find(paths.begin(), paths.end(), path);
	if (piter == paths.end())
	  paths.push_back(path);
      }
      
      void erase(const path_type& path)
      {
	if (path.empty()) return;
	
	SignalBlocker __blocker;
	
	lock_type lock(__mutex);

	path_set_type::iterator piter = std::find(paths.begin(), paths.end(), path);
	if (piter != paths.end())
	  paths.erase(piter);
      }

      void clear()
      {
	SignalBlocker __blocker;
	
	lock_type lock(__mutex);
	
	if (paths.empty()) return;
	
	for (path_set_type::const_iterator piter = paths.begin(); piter != paths.end(); ++ piter) {
	  try {
	    if (boost::filesystem::exists(*piter))
	      utils::filesystem::remove_all(*piter);
	    errno = 0;
	  }
	  catch (...) { }
	} 
	paths.clear();
      }
      
      path_set_type paths;
    };
    
    // signal management...
    struct SignalSet
    {
      typedef struct sigaction sigaction_type;
      //typedef std::vector<sigaction_type, std::allocator<sigaction_type> > sigaction_set_type;
      
      SignalSet() {}

      sigaction_type& operator[](size_t pos)
      {
	return signals[pos];
      }
      
      const sigaction_type& operator[](size_t pos) const
      {
	return signals[pos];
      }
      
      //void reserve(size_t __size) { signals.reserve(__size); }
      //void resize(size_t __size) { signals.resize(__size); }
      
      sigaction_type signals[NSIG];
    };

    static PathSet   __paths;
    static SignalSet __signals;
    
    static boost::once_flag signal_installer_once = BOOST_ONCE_INIT;

    // callback functions...
    static void callback(int sig) throw()
    {
      // actual clear...
      __paths.clear();
      
      // is this safe without mutex???
      ::sigaction(sig, &__signals[sig], 0);
      ::kill(getpid(), sig);
    }
    
    struct SignalInstaller
    {
      static void initializer()
      {
	typedef struct sigaction sigaction_type;
	
	sigaction_type sa;
	sa.sa_handler = callback;
	sa.sa_flags = SA_RESTART;
	sigemptyset(&sa.sa_mask);
      
	::sigaction(SIGHUP,  &sa, &__signals[SIGHUP]);
	::sigaction(SIGINT,  &sa, &__signals[SIGINT]);
	::sigaction(SIGQUIT, &sa, &__signals[SIGQUIT]);
	::sigaction(SIGILL,  &sa, &__signals[SIGILL]);
	::sigaction(SIGABRT, &sa, &__signals[SIGABRT]);
	//::sigaction(SIGKILL, &sa, &__signals[SIGKILL]);
	::sigaction(SIGSEGV, &sa, &__signals[SIGSEGV]);
	::sigaction(SIGTERM, &sa, &__signals[SIGTERM]);
	::sigaction(SIGBUS,  &sa, &__signals[SIGBUS]);

      }

      SignalInstaller() 
      {
	boost::call_once(signal_installer_once, initializer);
      }
    };

    static SignalInstaller __signal_installer;
    
  };
  
  void tempfile::insert(const path_type& path)
  {
    tempfile_impl::__paths.insert(path);
  }
  
  void tempfile::erase(const path_type& path)
  {
    tempfile_impl::__paths.erase(path);
  }
  
  inline 
  std::string get_hostname()
  {
    char name[1024];
    int result = ::gethostname(name, sizeof(name));
    if (result < 0)
      throw std::runtime_error("gethostname()");
    
    return std::string(name);
  }
  
  inline
  std::string get_hostname_short()
  {
    std::string hostname = get_hostname();

    std::string::size_type pos = hostname.find('.');
    if (pos != std::string::npos)
      return hostname.substr(0, pos);
    else
      return hostname;
  }
  
  
  tempfile::path_type tempfile::tmp_dir()
  {
    const char* tmpdir_spec_env = getenv("TMPDIR_SPEC");
    if (tmpdir_spec_env) {
      std::string tmpdir_spec(tmpdir_spec_env);
      
      const std::string tmpdir_spec_short = boost::algorithm::replace_all_copy(tmpdir_spec, "%host", get_hostname_short());
      const std::string tmpdir_spec_long = boost::algorithm::replace_all_copy(tmpdir_spec, "%host", get_hostname());

      
      if (boost::filesystem::exists(tmpdir_spec_short) && boost::filesystem::is_directory(tmpdir_spec_short))
	return tmpdir_spec_short;
      
      if (boost::filesystem::exists(tmpdir_spec_long) && boost::filesystem::is_directory(tmpdir_spec_long))
	return tmpdir_spec_long;
    }
    
    // wired to /tmp... is this safe?
    std::string tmpdir("/tmp");
    
    const char* tmpdir_env = getenv("TMPDIR");
    if (! tmpdir_env)
      return path_type(tmpdir);
      
    const path_type tmpdir_env_path(tmpdir_env);
    if (boost::filesystem::exists(tmpdir_env_path) && boost::filesystem::is_directory(tmpdir_env_path))
      return tmpdir_env_path;
    else
      return path_type(tmpdir);
    
    
  }
  
};
