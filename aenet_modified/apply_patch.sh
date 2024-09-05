#!/bin/bash
sleep 1
cp -r src src_modified
patch -u -p1 -d src_modified < src.patch
