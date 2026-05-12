The intent is to make as much as possible in one file -- put the solar
system in one file with all the spice necessary to make it work.

* VoyagerScene.inc - Set up time and camera variables, then include
  this file as the last line in your specific scene
* Voyager2Uranus.pov - time and camera variables for the Voyager 2
  Uranus flyby

## Imported
* SpiceQuat.inc - contains QuatTrans, which uses a quaternion from pxform
  to transform an object
* KwanMath.inc - Lots of stuff, but PrintQuat for now.
* StarsRight.inc - draw the stars in the sky. Stars are in the ICRF/J2000 
  frame, in a right handed frame with 0hRA at +X, 6hRA at +Y and +90dec at +Z. 
  * Catalog.inc - Yale Bright Star Catalog in the form of a
    string for each line in the catalog, sorted by brightness 