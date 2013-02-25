=======
lattice
=======

# A description for JSON lattice format

:Author: Taro Watanabe <taro.watanabe@nict.go.jp>
:Date:   2013-2-11

## JLF Format

The native lattice format is JLF (JSON Lattice Format) which is
represented by the [JSON data format](http://www.json.org>).
Strings must be escaped like a C string, (see JSON specification), and
integers and floating points are discriminated by whether a fractional
point "." is included, or not.
	
	lattice    ::= [arcs *(, arcs)]
	arcs       ::= [arc *(, arc)]
	arc        ::= [label, features(, attributes), distance]
	label      ::= JSON-STRING
	features   ::= {} | { JSON-STRING : JSON-FLOAT *(, JSON-STRING : JSON-FLOAT) }
	attributes ::= { JSON-STRING : value *(, JSON-STRING : value) }
	distance   ::= JSON-INT
	value      ::= JSON-STRING | JSON-FLOAT | JSON-INT
	

## PLF Format

cicada can read PLF (Python Lattice Format) used in [Moses](http://statmt.org/moses/), but do not
support writing.

	lattice  ::= (+(arcs,))
	arcs     ::= (+(arc,))
	arc      ::= (label, cost, distance)
	label    ::= PYTHON-STRING
	cost     ::= PYTHON-FLOAT
	distance ::= PYTHON-INT

The cost in the PLF is interpreted as a "lattice-cost" feature.

## Examples:

JLF:

	[[["ein'\"en", {"lattice-cost": 1.0}, 1]], [["wettbewerbsbedingten", {"lattice-cost": 0.5}, 2], ["wettbewerbs", {"lattice-cost": 0.25}, 1], ["wettbewerb", {"lattice-cost": 0.25}, 1]], [["bedingten", {"lattice-cost": 1.0}, 1]], [["preissturz", {"lattice-cost": 0.5}, 2], ["preis", {"lattice-cost": 0.5}, 1]], [["sturz", {"lattice-cost": 1.0}, 1]]]

PLF:

	((('ein\\"en',1.0,1),),(('wettbewerbsbedingten',0.5,2),('wettbewerbs',0.25,1), ('wettbewerb',0.25, 1),),(('bedingten',1.0,1),),(('preissturz',0.5,2), ('preis',0.5,1),),(('sturz',1.0,1),),)


For detail on PLF, see Moses web pages. In JLF, we can easily add extra features encoded in the lattice,
by inserting new one in {}, such as

	{"lattice-cost": 0.5, "accoustic": -5000.9}

Remark, the costs are interpreted as logarithmic value so that we can compute score by "weight \cdot feature-function".

We allow <epsilon> node which allows us to directly handle confusion nework.
Inside cicada, it is recommended to perform remove-epsilon, before intersecting with grammar 

## Tools

cicada\_unite\_lattice

> Merge multiple lattices (or sentences) in one. Here, we simply compute the union sharing start/goal states.
> Optionally dump output in graphviz format (--graphviz flag).

cicada\_unite\_sentence

> Merge multiple sentences into one confusion-network via TER alignment.
> Support incrementat merging by confusion-network-TER. (--merge flag)
> TER-conputation by lower-cased word (--lower flag)
> Dump in graphviz format (--graphviz flag)

