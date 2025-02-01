Images from the identified sequences as shown in the Voyager Uranus
Travel Guide, and highlighted in Blinn's animation.

These are all downloaded from OPUS with the following boxes checked:
* Voyager ISS-Specific Products
   * ( ) Raw image
   * ( ) Cleaned Image
   * ( ) Calibrated Image
   * (/) Geometrically Corrected Image
   * (/) Reseau Table
   * (/) Geometric Tiepoint Table
   * ( ) Extra preview (raw)
   * ( ) Extra preview (cleaned)
   * ( ) Extra preview (calibrated)
   * (/) Extra preview (geometrically corrected)
   * ( ) Documentation [only selected for VOBEST, would otherwise be 20+MB of the same thing every time]
* Metadata products [All selected for VOBEST and deselected for others. All the same for any image on volume VGISS_7206, which all of the near-encounter stuff is]
   * ( ) RMS Node Augmented Index
   * ( ) Raw Image Index (contains a column describing the intent of each image)
   * ( ) Supplemental Index
   * ( ) Target Body Inventory
   * ( ) Planet Geometry Index [For poor unfortunate souls with no access to Spice]
   * ( ) Moon Geometry Index [Same]
   * ( ) Ring Geometry Index [Same. It says ring, means ring plane geometry. No actual data on rings which are physically present]
* Browse Products
   * ( ) Browse Image (thumbnail)
   * ( ) Browse Image (small)
   * ( ) Browse Image (medium)
   * ( ) Browse Image (full)

Each individual folder includes:

* `manifest.csv` - files downloaded in this shopping cart. OPUS uses
  the concept of a "shopping cart" to select and download only selected
  images.
* `data.csv` - File containing one row for each image in the cart.
  This indluces the intended target, observation time, and exposure duration.
* `url.csv` - URL of each file in the package that isn't generated specificall for this download.
* <basename>.LBL - Each 
* Cnnnnnnn_<stuff>.<ext> - Each observation is identified by a 7-digit number
  which monotonically increases with observation time.
  * Cnnnnnnn_full.jpg - Browse image of raw data. Full resolution but uncalibrated/uncorrected
  * Cnnnnnnn_GEOMA.DAT - Resseau markings in binary format as used by Voyager science pipeline.
    This shows the observed and expected locations of each Resseau mark.
    * Cnnnnnnn_GEOMA.TAB - Same data in human-readable ASCII 
  * Cnnnnnnn_GEOMED.IMG - Calibrated and geometrically corrected image. Pixel
    values are proportional to best estimate of true I/F albedo and reprojected
    into a space with known geometry.
  * Cnnnnn_RESLOC.DAT 

# VOBEST
Best resolution and color images of Oberon. 
* C2683623,1986-01-24T08:48:44.080
* C2683625,1986-01-24T08:50:18.160
* C2683627,1986-01-24T08:51:54.160
* C2683629,1986-01-24T08:53:30.160

# VTCOLOR
Color images of Titania, not featured with a zoom-in in Blinn's movie
* C2683649,1986-01-24T09:09:31.080
* C2683651,1986-01-24T09:11:05.160
* C2683653,1986-01-24T09:12:41.160
* C2683655,1986-01-24T09:14:17.160

# VUBEST
Best resolution of Umbriel
* C2684004,1986-01-24T11:45:31.080
* C2684006,1986-01-24T11:47:05.160


