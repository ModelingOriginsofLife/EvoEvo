#!/bin/bash

bash cmake_clean.sh
cmake -DCMAKE_BUILD_TYPE=Debug -DWITH_GRAPHICS_CONTEXT=1 ..
make
