#!/usr/bin/env python
#
#  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
#

import time
import sys
import os, os.path
import string

from optparse import OptionParser, make_option

opt_parser = OptionParser(
    option_list=[
    # output directory/filename prefix
    make_option("--byte-size",  action="store", type="int"),
    ])

(options, args) = opt_parser.parse_args()

bit_size_max = options.byte_size * 8

print "#ifndef __UTILS__BITPACK%d_IMPL__HPP__" %(bit_size_max)
print "#define __UTILS__BITPACK%d_IMPL__HPP__ 1" %(bit_size_max)
print "namespace utils {"
print "namespace bitpack {"

for bit_size in range(1, bit_size_max):
    
    unit_size = options.byte_size * 8
    bits_total = options.byte_size * 8 * bit_size
    loop_total = options.byte_size * 8
    mask = (1L << bit_size) - 1

    mask_mod = ""
    if mask >= 0xffffffffL:
        mask_mod = "ull"
    
    print "template <typename Tp>"
    print "struct __struct_bitpack_impl<Tp,%d,%d>" %(options.byte_size, bit_size)
    print "{"
    print "  static const size_t unit_size = %d;"  %(unit_size)
    print "  static const size_t bits_total = %d;" %(bits_total)
    print "  static const size_t loop_total = %d;" %(loop_total)
    print "  static const size_t bit_size = %d;" %(bit_size)
    print "  static const Tp     mask = 0x%02x%s;" %(mask, mask_mod)
    print "  static inline"
    print "  void pack(const Tp* source, Tp* destination)"
    print "  {"

    for i in range(loop_total):
        shift_size_total = bits_total - bit_size * (i + 1)
        code_pos = bit_size - 1 - (shift_size_total / unit_size)
        shift_amount = shift_size_total & (unit_size - 1)

        if (unit_size - shift_amount) < bit_size:
            print "    destination[%d] |= (source[%d] & 0x%02x%s) >> %d;" %(code_pos - 1, i, mask, mask_mod, unit_size - shift_amount)

        if (unit_size - shift_amount) <= bit_size:
            print "    destination[%d]  = (source[%d] & 0x%02x%s) << %d;" %(code_pos, i, mask, mask_mod, shift_amount)
        else:
            print "    destination[%d] |= (source[%d] & 0x%02x%s) << %d;" %(code_pos, i, mask, mask_mod, shift_amount)
    
    print "  }"
    print "  static inline"
    print "  void unpack(const Tp* source, Tp* destination)"
    print "  {"
    
    for i in range(loop_total):
        shift_size_total = bits_total - bit_size * (i + 1)
        code_pos = bit_size - 1 - (shift_size_total / unit_size)
        shift_amount = shift_size_total & (unit_size - 1)

        if (unit_size - shift_amount) < bit_size:
            print "    destination[%d]  = (source[%d] << %d) & 0x%02x%s;" %(i, code_pos - 1, unit_size - shift_amount, mask, mask_mod)

        if (unit_size - shift_amount) < bit_size:
            print "    destination[%d] |= (source[%d] >> %d) & 0x%02x%s;" %(i, code_pos, shift_amount, mask, mask_mod)
        else:
            print "    destination[%d]  = (source[%d] >> %d) & 0x%02x%s;" %(i, code_pos, shift_amount, mask, mask_mod)
    
    print "  }"
    print "};"

print "};"
print "};"
print "#endif"
