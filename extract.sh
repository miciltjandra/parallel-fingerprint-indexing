#!/bin/bash
for ((i=0; i<5; i++)); do
    for filename in CleanFingerprints/*.bmp; do
        ./feature_extraction $filename
    done

    for filename in NoisyFingerprints/*.bmp; do
        ./feature_extraction $filename
    done
done