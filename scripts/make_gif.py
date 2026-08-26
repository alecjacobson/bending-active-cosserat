#!/usr/bin/env python3
"""Assemble a GIF from a sequence of PNG frames.
Usage: make_gif.py <out.gif> <fps> <frame1.png> <frame2.png> ...
"""
import sys
from PIL import Image

out = sys.argv[1]
fps = float(sys.argv[2])
paths = sys.argv[3:]
if not paths:
    print("no frames")
    sys.exit(1)

frames = [Image.open(p).convert("RGB") for p in paths]
dur = int(1000 / fps)
frames[0].save(out, save_all=True, append_images=frames[1:], duration=dur, loop=0,
               optimize=True)
print(f"wrote {out} ({len(frames)} frames @ {fps} fps)")
