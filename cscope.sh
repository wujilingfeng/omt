find ./ /home/libo/Documents/c++/Viewer/include /home/libo/Documents/c++/libcell/include -name '*.h' -o -name '*.c' -o -name '*.cpp' > ./cscope.files&&cscope -Rbq -P `pwd`
