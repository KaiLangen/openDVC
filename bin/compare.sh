#!/bin/bash

sourcevid=$1
recvid=$2

ffmpeg -y -f rawvideo -vcodec rawvideo -pix_fmt yuv420p -s 176x144 -i $sourcevid \
          -f rawvideo -vcodec rawvideo -pix_fmt yuv420p -s 176x144 -i $recvid -lavfi \
"setpts=PTS-STARTPTS,psnr" -f null - 2>&1 \
| grep -i parsed
#"setpts=PTS-STARTPTS,split=2[a][b],[a][1:v]ssim;[b][1:v]psnr" -f null - 2>&1 \
