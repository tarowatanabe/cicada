Option format
=============

Each feature function or grammar file (and many other components) can be specified via
option format
::

  name : option1 = value1, option2=value2, option3=value3, option4=value4

You may want to escape by double quotes when you want to specify names
or options with ``:``, ``=``, ``,``:
::

  "name" : "option1"="value1", ....

In this case, the double quoted strings follow "C" escaping,
i.e. ``\n`` for new-line.
