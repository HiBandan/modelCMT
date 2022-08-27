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
  
  ### windows
  
    copy "Eigen" folder to C:\Program Files 
  
  ### linux 
  
    copy "Eigen" folder to /usr/local/include

## install boost:
download boost: https://drive.google.com/file/d/1apu5_am2kJj7HvXJNPhn30k_ryi3gJDf/view?usp=sharing

  extract "boost_1_77_0" folder and rename it to "Boost"
  
  ### windows
  
    copy "Boost" to C:\Program Files 
  
  ### linux 
  
    copy "Boost" to /usr/local/include
    
## run 

  1. open -> compile.CB (in code-block)
  
  2. settings -> compiler -> search directories (compiler) -> add 

      -> C:\Program Files\Eigen
  
      -> C:\Program Files\Boost

  3. project -> set program’s arguments
  
  4. run 
