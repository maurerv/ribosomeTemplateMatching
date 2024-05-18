#!/bin/bash

rsync -am --progress \
  --include='coordinates_annotated.tsv' \
  --include='*/' \
  --exclude='*' \
  "embl:/g/kosinski/vmaurer/ribosomePaper/templates_*" ../
