The Houdini Ocean Toolkit (HOT)
===============================
 

_Notice: Houdini versions >= 12.5 now include functionality (Houdini Ocean **Tools**) equivalent/superior to the HOT and hence this code won't be supported
from this point on. I'll leave it up on github as the source code may be useful for perusal._

_The first port of the code outside of Houdini was to 3DStudio Max. The guy asked me if it would be ok, and what he should call it. I jokingly suggested HOT for Max, and it became Hot4Max, then Hot4LW etc. Eventually I was seeing references all over the place to the Houdini Ocean Toolkit for <some 3d software>, which was a bit embarrasing, but the SideFX people thankfully never seemed to have a problem with it. I hoped it would make people more curious to try Houdini which back in 2011 was nowhere near as popular as it is now. All up it was a blast!_

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

The HOT was used quite a bit in Houdini, but also ported to many (most?) of the other 3D systems. 

Here are a few references. 

PS I'd love to hear of any other places it was used. (email: Drew.Whitehouse@gmail.com)

* Happy Feet2 https://vimeo.com/98038762 (0:36)
* Storm Studio's Kon-Tiki VFX https://vimeo.com/48314160, https://vimeo.com/25545997
* Rising Sun Pictures used the HOT on Terminator Salvation, see Gregory Yepes's web site. https://gregoryyepes.com/#projects
* Sam Loxton's Ocean Storm.
* Some blender guys have ported the code to C https://mattebb.com/blog/weblog/ocean-simulation/, https://docs.blender.org/manual/en/latest/modeling/modifiers/physics/ocean.html
* It looks like this XSI guy has also used the HOT as a base.
* Framestore used the HOT in the award winning Smirnoff Sea commercial - cool! https://www.cgw.com/Press-Center/Web-Exclusives/2007/Framestore-CFC-Works-on-Smirnoff-Shot.aspx, https://www.broadcastnow.co.uk/british-companies-snare-vfx-awards/658050.article, https://www.digitalartsonline.co.uk/news/creative-lifestyle/framestore-cleans-up-sea-for-smirnoff/
* Ocean whirlpool WIP on the forum.
* The HOT has been ported to Lightwave, see this newsletter for an interview with the author.https://www.facebook.com/groups/lightwiki/permalink/1773119119460746/
* Kenneth A Huff used the HOT for his art work, see his blog entry.
* Houdini ocean toolkit videos on vimeo.
* The HOT for MAX, some example animations on this youtube channel.
* The HOT for Maya, see this forum and this download.
* More HOT for Maya http://nico-rehberg.de/shader.html
* HOT for Maya tutorial video.
* Naiad plus HOT video by Igor Zanic.
* 15 Tips for impressive fluid effects from 3D World Magazine.
* Water Foam Generation Tests by Christian Schnellhammer, using his own port to Maya.
* Cinema 4D https://www.studiodaily.com/2015/02/framestores-graphics-send-evironmental-message-mission-blue/, https://www.cgrecord.net/2013/01/download-houdini-ocean-toolkit-v03-for.html

(link to the original code archive https://code.google.com/archive/p/houdini-ocean-toolkit/)

