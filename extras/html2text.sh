#!/bin/sh
w3m -dump -cols 78 $1 \
  | sed -e '1,4d' \
  | sed -e '/SourceForge Logo/,$d' \
  | recode utf8..ascii