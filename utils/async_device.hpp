// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__ASYNC_DEVICE__HPP__
#define __UTILS__ASYNC_DEVICE__HPP__ 1

#include <cstring>
#include <cerrno>

#include <sys/types.h>
#include <unistd.h>

#include <stdexcept>

#include <boost/iostreams/categories.hpp>
#include <boost/shared_ptr.hpp>

namespace utils
{
  
  class async_sink
  {
  public:
    typedef char char_type;
    struct category : public boost::iostreams::sink_tag,
		      public boost::iostreams::closable_tag {};
    
    async_sink(int fd, bool close_on_exit=false) : pimpl(new impl()) { open(fd, close_on_exit); }
    async_sink() : pimpl(new impl()) { }
    
    void open(int fd, bool close_on_exit=false) { pimpl->open(fd, close_on_exit); }
    void close() { pimpl->close(); }
    std::streamsize write(const char_type* s, std::streamsize n) { return pimpl->write(s, n); }
    
  private:
    struct impl
    {
      impl() : fd(-1), close_on_exit(false) {}
      ~impl() { close(); }

      std::streamsize write(const char_type* s, std::streamsize n)
      {
	std::streamsize offset = 0;
	
	while (offset < n) {
	  const std::streamsize written = ::write(fd, s + offset, n - offset);
	  
	  if (written < 0)
	    throw std::runtime_error(std::string("write(): ") + strerror(errno));
	  offset += written;
	}
	
	return offset;
      }

      void open(int _fd, bool _close_on_exit=false)
      {
	close();
	
	fd = _fd;
	close_on_exit = _close_on_exit;
      }

      void close()
      {
	if (fd >= 0)
	  ::fsync(fd);
	
	if (fd >= 0 && close_on_exit)
	  ::close(fd);
	fd = -1;
	close_on_exit = false;
      }
      
      int  fd;
      bool close_on_exit;
    };

    boost::shared_ptr<impl> pimpl;
  };
  

  class async_source
  {
  public:
    typedef char char_type;
    struct category : public boost::iostreams::source_tag,
		      public boost::iostreams::closable_tag {};
    
    async_source(int fd, bool close_on_exit=false) : pimpl(new impl()) { open(fd, close_on_exit); }
    async_source() : pimpl(new impl()) { }
    
    void open(int fd, bool close_on_exit=false) { pimpl->open(fd, close_on_exit); }
    void close() { pimpl->close(); }
    std::streamsize read(char_type* s, std::streamsize n) { return pimpl->read(s, n); }
    
  private:
    struct impl
    {
      impl() : fd(-1), close_on_exit(false) {}
      ~impl() { close(); }
      
      std::streamsize read(char_type* s, std::streamsize n)
      {
	if (fd < 0) return -1;

	int offset = 0;
	
	while (offset < n) {
	  const int read_size = ::read(fd, s + offset, n - offset);
	  
	  if (read_size == -1) {
	    if (errno == EAGAIN) {
	      
	      continue;
	    } else
	      throw std::runtime_error(std::string("read(): ") + strerror(errno));
	  } else if (read_size == 0) {
	    close();
	    
	    return (offset == 0 ? -1 : offset);
	  } 
	  offset += read_size;
	}

	return offset;
      }

      void open(int _fd, bool _close_on_exit=false)
      {
	close();
	
	fd = _fd;
	close_on_exit = _close_on_exit;
      }

      void close()
      {
	if (fd >= 0)
	  ::fsync(fd);
	
	if (fd >= 0 && close_on_exit)
	  ::close(fd);
	fd = -1;
	close_on_exit = false;
      }
      
      int  fd;
      bool close_on_exit;
    };

    boost::shared_ptr<impl> pimpl;
  };
};


#endif
