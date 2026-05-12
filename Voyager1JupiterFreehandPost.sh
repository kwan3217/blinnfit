# Postprocess
composite data/frames/Voyager1Jupiter/frame$(printf "%04d" $1).png data/frames/pov/Voyager1JupiterFreehand/frame_$(printf "%04d" $1).png -blend 100x50 data/frames/composite/Voyager1Jupiter/composite_$(printf "%04d" $1).png
