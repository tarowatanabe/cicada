// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__SUBPROCESS__HPP__
#define __UTILS__SUBPROCESS__HPP__ 1

#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <fcntl.h>
#include <signal.h>

#include <cerrno>
#include <cstring>

#include <string>
#include <stdexcept>

#include <boost/filesystem.hpp>

namespace utils
{
  
  class subprocess
  {
    
  public:
    explicit subprocess(const std::string& sh_command)
      : __pid(-1), __pread(-1), __pwrite(-1) { open(sh_command); }
    explicit subprocess(const boost::filesystem::path& command)
      : __pid(-1), __pread(-1), __pwrite(-1) { open(command); }
    ~subprocess() { close(); }
  private:
    subprocess(const subprocess& x) {}
    subprocess& operator=(const subprocess& x) { return *this; }
    
  public:
    int desc_read() const { return __pread; }
    int desc_write() const { return __pwrite; }
    int& desc_read() { return __pread; }
    int& desc_write() { return __pwrite; }
    

    void terminate()
    {
      if (__pid >= 0)
	::kill(__pid, SIGTERM);
    }

    void open(const std::string& sh_command)
    {
      int pin[2] = {-1, -1};
      int pout[2] = {-1, -1};
    
      if (::pipe(pin) < 0)
	throw std::runtime_error(std::string("pipe(): ") + strerror(errno));
      if (::pipe(pout) < 0) {
	::close(pin[0]);
	::close(pin[1]);
	throw std::runtime_error(std::string("pipe(): ") + strerror(errno));
      }
      
      const pid_t pid = ::fork();
      if (pid < 0) {
	::close(pin[0]);
	::close(pin[1]);
	::close(pout[0]);
	::close(pout[1]);
	throw std::runtime_error(std::string("fork(): ") + strerror(errno));
      }
      
      if (pid == 0) {
	// child process...
	// redirect input...
	::close(pin[1]);
	::dup2(pin[0], STDIN_FILENO);
	::close(pin[0]);
	
	// redirect output...
	::close(pout[0]);
	::dup2(pout[1], STDOUT_FILENO);
	::close(pout[1]);
	
	::execlp("sh", "sh", "-c", sh_command.c_str(), (char*) 0);
	
	::_exit(errno);  // not exit(errno)!
      } else {
	// parent process...
	::close(pin[0]);
	::close(pout[1]);
	
	__pid = pid;
	__pread = pout[0];
	__pwrite = pin[1];
      }
    }

    void open(const boost::filesystem::path& command)
    {     
      int pin[2] = {-1, -1};
      int pout[2] = {-1, -1};
    
      if (::pipe(pin) < 0)
	throw std::runtime_error(std::string("pipe(): ") + strerror(errno));
      if (::pipe(pout) < 0) {
	::close(pin[0]);
	::close(pin[1]);
	throw std::runtime_error(std::string("pipe(): ") + strerror(errno));
      }
      
      const pid_t pid = ::fork();
      if (pid < 0) {
	::close(pin[0]);
	::close(pin[1]);
	::close(pout[0]);
	::close(pout[1]);
	throw std::runtime_error(std::string("fork(): ") + strerror(errno));
      }
      
      if (pid == 0) {
	// child process...
	// redirect input...
	::close(pin[1]);
	::dup2(pin[0], STDIN_FILENO);
	::close(pin[0]);
	
	// redirect output...
	::close(pout[0]);
	::dup2(pout[1], STDOUT_FILENO);
	::close(pout[1]);
	
	::execlp(command.file_string().c_str(), command.file_string().c_str(), (char*) 0);
	
	::_exit(errno);  // not exit(errno)!
      } else {
	// parent process...
	::close(pin[0]);
	::close(pout[1]);
	
	__pid = pid;
	__pread = pout[0];
	__pwrite = pin[1];
      }
    }
    
    void close() {
      if (__pwrite >= 0) {
	::close(__pwrite);
	__pwrite = -1;
      }
      
      if (__pread >= 0) {
	::close(__pread);
	__pread  = -1;
      }
      
      if (__pid >= 0) {
	int status = 0;
	::waitpid(__pid, &status, 0);
	__pid = -1;
      }
    }
    
    void __close_on_exec(int fd)
    {
      int flags = ::fcntl(fd, F_GETFD, 0);
      if (flags == -1)
	throw std::runtime_error(std::string("fcntl(): ") + strerror(errno));
      flags |= FD_CLOEXEC;
      if (::fcntl(fd, F_SETFD, flags) == -1)
	throw std::runtime_error(std::string("fcntl(): ") + strerror(errno));
    }
    
  public:
    pid_t __pid;
    int   __pread;
    int   __pwrite;
  };
  
};

#endif
