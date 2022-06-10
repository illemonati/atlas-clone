
#!/bin/bash

echo -e "checking for GCC version..."
version=$(gcc --version | awk 'NR==1 {print $NF}')
echo "GCC version: $version"
versionMajor=$(echo $version | cut -d. -f1)
if [ $versionMajor -gt 7 ]; then echo "GCC version is sufficient for ATLAS"; else echo "ERR: please use GCC version 8 or higher."; exit 1 ;fi


mkdir -p build
cd build

cmake=cmake


# define he help
usage='
\n \******************************COMPILING ATLAS*****************************************\
\n
\n\t This script compiles the atlas executable.
\n\t Usage: bash compile_atlas.sh [options]
\n\t 
\n\t 
\n\t Optional parameters:
\n\t  -h \t open this page
\n\t 
\n\t  -c \t specify the full path to a cmake executable. 
\n\t     \t Required version: 2.14 or higher
\n\t     \t Default: cmake
\n\t 
\n\t  -a \t specify the full path to a locally installed Armadillo library.
\n\t 
\n\t 
\n \**************************************************************************************\
'

while getopts ":hc:a:" OPTION
do
    case $OPTION in
    h)
      echo -e $usage
	    exit 1
	    ;;
    c)
      cmake=${OPTARG}
      if [[ -f $cmake ]]; then
        echo -e "Using custom cmake version: $cmake \n"
      else
        echo -e "ERROR: File $cmake specified with -c does not exist. Please specify the full path to a cmake executable \n"
        echo -e $usage
        exit 1
      fi
		;;
    a)
      arm=${OPTARG}
      if [[ -d "$arm" ]]; then
        echo -e "Using custom armadillo version from: $arm \n"
      else
        echo -e "ERROR: Directory $arm specified with -a does not exist. Please specify the full path to the directory of your armadillo library \n" 
        echo -e $usage
        exit 1
      fi
    ;;
    \?)
      echo -e "ERROR: no valid parameters defined! \n\n"
      echo -e $usage
      exit 1
      ;;
    :)
      echo "option requires an argument -- $OPTARG" 
      echo -e $usage
      exit 1
      ;;
    esac
done

echo -e "checking for cmake version..."
version=$($cmake --version | awk 'NR==1 {print $NF}')
echo "cmake version: $version"
version1=$(echo $version | cut -d. -f1)
version2=$(echo $version | cut -d. -f2)
if ( [ $version1 == 2 ] && [ $version2 -gt 13 ]) || [ $version1 -gt 2 ] ; then echo "Cmake version is sufficient for ATLAS"; else echo "ERR: please use cmake version 2.14 or higher."; exit 1 ;fi



if [ -z $arm ]; then
  $cmake ../
  make
else
  $cmake -DCMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH};$arm"
  make
fi

echo "done compiling. \n"

cd ..