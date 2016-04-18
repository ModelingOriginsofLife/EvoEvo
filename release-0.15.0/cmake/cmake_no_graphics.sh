#!/bin/bash

bash cmake_clean.sh
cmake -DCMAKE_BUILD_TYPE=Release -DWITH_GRAPHICS_CONTEXT=0 ..
make
