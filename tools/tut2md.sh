#! /bin/bash

tail -n +6 "$1" | awk '{if ( $1 == "#") {print substr($0, 3)} else {print $0}}'
