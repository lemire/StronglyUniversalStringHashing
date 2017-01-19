cd "${0%/*}"
./powerpolicy.sh  performance
./disablehyperthreading.sh 
./turboboost.sh off
