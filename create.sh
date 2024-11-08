# Only run this at the start of each class!
now="$(date +'%m-%d-%y')"
mkdir $now
cp -n templates/hello.cpp $now/$now.cpp
cp -n templates/notes.md $now/