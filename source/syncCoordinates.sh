#!/bin/bash

rsync -am --progress \
  --include='coordinates.tsv' \
  --include='*/' \
  --exclude='*' \
  embl:/g/kosinski/vmaurer/ribosomePaper/templates ../
