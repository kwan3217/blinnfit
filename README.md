# blinnfit
Find the camera position and pointing by 
fitting to the stars and planets

# VoyagerUranusHD
This is a rip by ILoveSaturn of the best copy of the Voyager Uranus flyby animation that I have ever seen. 
It contains 6906 frames (projected at 30FPS) numbered from 0001 to 6906. Some areas are harder to fit than
others because there are so few stars. These are generally zoom-ins but I have also been 
having problems with some zoom-outs.

* 1-128 Zoom out, including orbit rings and moon names
* 487-984 Occultation of SigSag by rings
* 1314-1596 VOBEST (Imaging Oberon)
* 2264-2545 VUBEST (Imaging Umbriel)
* 3020-3288 VTBEST (Imaging Titania)
* 3768-4084 VABEST (Imaging Ariel)
* 4204-4664 VMBEST (Imaging Miranda)
* 5264-6906 - Only every 10th image (last digit 0) is done to save time

# Example: Fitting the Voyager 1 Jupiter encounter
There is a view on YouTube collected by ILikeSaturn. I think he and I could be friends. The video is at:
[Voyager 1 flyby Animation - Jupiter and Moons (1979)](https://www.youtube.com/watch?v=KfCij7iTz3U):

* Download it with `yt-dlp`. Get this program from github https://github.com/yt-dlp/yt-dlp/releases:
```
cd data/frames
yt-dlp https://www.youtube.com/watch?v=KfCij7iTz3U
```
* Unpack it into frames
```
mkdir Voyager1Jupiter
ffmpeg -i Voyager\ 1\ flyby\ Animation\ -\ Jupiter\ and\ Moons\ \(1979\)\ \[KfCij7iTz3U\].mkv -y Voyager1Jupiter/frame%04d.png
```
* Run starfit:
```
python src/starfit/__init__.py -f 1000 Voyager1Jupiter
```
The first few will fail because it tries to interpolate or extrapolate, and can't because there isn't enough data yet
to do so.

It doesn't look like the stars match up well with reality for V1J, so we will time it
with crossings of moons (mainly Amalthea) for timing, and positions of moons for 
