#!/bin/bash

URL_ROOT="http://archive.gersteinlab.org/proj/encodenets"
DISTAL_URL="${URL_ROOT}/Distal.txt"
PROXIMAL_URL="${URL_ROOT}/Proximal_filtered.txt"

DISTAL_FILE="Distal.txt"
PROXIMAL_FILE="Proximal_filtered.txt"

OUTPUT_FILE="tf-targets.txt"

curl "$DISTAL_URL"   | perl -pe 's/\r\n/\n/g; s/\s+(.)/\t$1/gm;' > "$DISTAL_FILE"
curl "$PROXIMAL_URL" | perl -pe 's/\r\n/\n/g; s/\s+(.)/\t$1/gm'  > "$PROXIMAL_FILE"

cat "$DISTAL_FILE" "$PROXIMAL_FILE" > "$OUTPUT_FILE"
