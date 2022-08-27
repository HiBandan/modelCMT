# modelCMT

## install gcc/g++

  ### windows
  
    download the latest version of Mingw: https://nuwen.net/mingw.html
  
    run mingw-xx.exe -> this will create a folder MinGw
      
    copy MinGw to C:\
       
    add environment variable: C:\MinGW or C:\MinGW\bin

  ### linux
  
    sudo apt install build-essential

## install eigen:
download eigen: https://drive.google.com/file/d/1VFLlKJI9EaSP9VUaeV8vmPQYkqdIIJPg/view?usp=sharing
  
  extract "eigen-3.4.0" folder and rename it to "Eigen"
  
  copy "Eigen" folder to C:\Program Files (windows) OR /usr/local/include (linux)

## install boost:
download boost: https://drive.google.com/file/d/1apu5_am2kJj7HvXJNPhn30k_ryi3gJDf/view?usp=sharing

  extract "boost_1_77_0" folder and rename it to "Boost"
  
  copy "Boost" to C:\Program Files (windows) OR /usr/local/include (linux)

## run 

  ### code::Block
  
  1. settings -> compiler -> search directories (compiler) -> add 

      -> C:\Program Files\Eigen
  
      -> C:\Program Files\Boost

  2. project -> set program’s arguments
