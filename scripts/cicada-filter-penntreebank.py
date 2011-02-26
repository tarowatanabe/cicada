#!/usr/bin/env python

import sys
import re

patterns = [
    re.compile(r"<DOC>"),
    re.compile(r"</DOC>"),
    re.compile(r"<BODY>"),
    re.compile(r"</BODY>"),
    re.compile(r"<HEADER>"),
    re.compile(r"</HEADER>"),
    re.compile(r"<HEADLINE>"),
    re.compile(r"</HEADLINE>"),
    re.compile(r"<TEXT>"),
    re.compile(r"</TEXT>"),
    re.compile(r"<TURN>"),
    re.compile(r"</TURN>"),
    re.compile(r"<P.*?>"),
    re.compile(r"</P>"),
    re.compile(r"<S.*?>"),
    re.compile(r"</S>"),
    re.compile(r"<seg.*?>"),
    re.compile(r"</seg>"),
    re.compile(r"<segment.*?>"),
    re.compile(r"</segment>"),
    re.compile(r"<DATE.*?>.*?</DATE>"),
    re.compile(r"<DATE_TIME.*?>.*?</DATE_TIME>"),
    re.compile(r"<END_TIME.*?>.*?</END_TIME>"),
    re.compile(r"<ENDTIME.*?>.*?</ENDTIME>"),
    re.compile(r"<DOCID.*?>.*?</DOCID>"),
    re.compile(r"<DOCNO.*?>.*?</DOCNO>"),
    re.compile(r"<DOCTYPE.*?>.*?</DOCTYPE>"),
    re.compile(r"<DATETIME.*?>.*?</DATETIME>")
    ]

for line in sys.stdin:
    for pattern in patterns:
        line = pattern.sub('', line)
    if not line.strip(): continue
    
    sys.stdout.write(line)
