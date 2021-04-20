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

Sightings
========= 

The HOT was used in Houdini, but also ported to many (most?) of the other 3D systems. E.g.

* Rising Sun Pictures used the HOT on Terminator Salvation, see Gregory Yepes's web site.
* Sam Loxton's Ocean Storm.
* Some blender guys have ported the code to C.
* It looks like this XSI guy has also used the HOT as a base.
* Framestore used the HOT in the award winning Smirnoff Sea commercial - cool!
* Ocean whirlpool WIP on the forum.
* The HOT has been ported to Lightwave, see this newsletter for an interview with the author.
* Kenneth A Huff used the HOT for his art work, see his blog entry.
* Houdini ocean toolkit videos on vimeo.
* The HOT for MAX, some example animations on this youtube channel.
* The HOT for Maya, see this forum and this download.
* More HOT for Maya.
* A HOT for Lightwave.
* HOT for Maya tutorial video.
* Naiad plus HOT video by Igor Zanic.
* 15 Tips for impressive fluid effects from 3D World Magazine.
* Water Foam Generation Tests by Christian Schnellhammer, using his own port to Maya.
