#!/bin/bash
for filename in NoisyFingerprints/*.bmp; do
    ./feature_extraction $filename
done