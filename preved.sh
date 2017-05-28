mogrify -format jpg -quality 90 $1*.png
ffmpeg -i $1%d.jpg -qscale 4 $2.avi
