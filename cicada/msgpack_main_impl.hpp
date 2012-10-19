//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __MSGPACK_MAIN_IMPL__HPP__
#define __MSGPACK_MAIN_IMPL__HPP__ 

#include <iostream>

#include "utils/config.hpp"

#ifdef HAVE_MSGPACK_HPP
#include <msgpack.hpp>
#endif

template <typename Obj>
inline
void msgpack_test(const Obj& obj)
{
#ifdef HAVE_MSGPACK
  {
    // packing...
    msgpack::sbuffer sbuf;
    msgpack::pack(sbuf, obj);

    std::cerr << "packed size: " << sbuf.size() << std::endl;
    
    // deserialize it.
    msgpack::unpacked msg;
    msgpack::unpack(&msg, sbuf.data(), sbuf.size());
    
    // get an object...
    Obj back;
    msg.get().convert(&back);
    
    if (back != obj)
      std::cerr << "different msgpacked?" << std::endl;
  }
  
  {
    // streaming...
     msgpack::sbuffer buffer;
 
     msgpack::packer<msgpack::sbuffer> pk(&buffer);
     pk.pack(obj);

     std::cerr << "streamed packed size: " << buffer.size() << std::endl;
     
     // deserializes these objects using msgpack::unpacker.
     msgpack::unpacker pac;
     
     // feeds the buffer.
     pac.reserve_buffer(buffer.size());
     memcpy(pac.buffer(), buffer.data(), buffer.size());
     pac.buffer_consumed(buffer.size());
     
     // now starts streaming deserialization.
     Obj back;
     msgpack::unpacked result;
     while(pac.next(&result)) {
       result.get().convert(&back);

       if (back != obj)
	 std::cerr << "different msgpacked?" << std::endl;
     }
  }
#endif
}

#endif
