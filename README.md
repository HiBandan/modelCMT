# modelCMT

## windows

  ### install gcc/g++ 
  
  download the latest version of Mingw: https://nuwen.net/mingw.html
  
      run mingw-xx.exe -> this will create a folder MinGw
      
      copy MinGw to C:\
       
      add environment variable: C:\MinGW or C:\MinGW\bin

  ### install eigen:
  
  download eigen: https://drive.google.com/file/d/1VFLlKJI9EaSP9VUaeV8vmPQYkqdIIJPg/view?usp=sharing
  
      extract "eigen-3.4.0" folder and rename it to "Eigen"
  
      copy "Eigen" folder to C:\Program Files 
  
  ### install boost:
      download boost: https://drive.google.com/file/d/1apu5_am2kJj7HvXJNPhn30k_ryi3gJDf/view?usp=sharing

      extract "boost_1_77_0" folder and rename it to "Boost"
  
      copy "Boost" to C:\Program Files 
  
  ### run 

    open -> compile.CB (in code-block)
  
    settings -> compiler -> search directories (compiler) -> add 
  
        -> C:\Program Files\Eigen 
  
        -> C:\Program Files\Boost 

    project -> set program’s arguments: parameters_ARRAY.txt
  
    run -> 
