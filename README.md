HOT
===

*Notice: Houdini versions >= 12.5 now include functionality equivalent/superior to the HOT and hence this code won't be supported
from this point on. I'll leave it up on github as the source code may be useful for perusal.* It's been a blast!

Houdini Ocean Toolkit

This version of the HOT is strictly for versions of Houdini >= 12.0.
For legacy versions go to the old repository on the Google Code site.

There is currently no support for building on Windows, I would welcome
the contribution of a compile.bat script that does the equivalent
of compile.sh.

Compile and install with

> ./compile.sh

Note that
=======

> ./compile.sh -f

will rebuild and install without a time consuming recompile of the 3rd party libs (of course you must have compiled them once).

(Read the source of the compile.sh script for more build options.)

Test by loading the examples, start with SOP_Simple.hip.
