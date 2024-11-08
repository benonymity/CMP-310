now="$(date +'%m-%d-%y')"
# dir="$(pwd)"
# cd $now && g++ -std=c++11 $now.cpp
g++ -std=c++17 strassen.cpp
./a.out
rm a.out
# cd $dir