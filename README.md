# modelCMT

## eigen:
download the latest version of eigen:https://drive.google.com/file/d/1VFLlKJI9EaSP9VUaeV8vmPQYkqdIIJPg/view?usp=sharing
  
  extract "eigen-3.4.0" folder and rename it to "Eigen"
  
  copy "Eigen" folder to C:\Program Files 

## boost:
download the latest version of boost: https://drive.google.com/file/d/1apu5_am2kJj7HvXJNPhn30k_ryi3gJDf/view?usp=sharing

  extract "boost_1_77_0" folder and rename it to "Boost"
  
  copy "Boost" to C:\Program Files 


## code::Block
1. install codeBlock
2. create a new project with console application
3. link additional libraries (from #include directory): 
   (a) settings -> compiler -> search directories (compiler) -> add -> C:\Program Files\Eigen
   (b) settings -> compiler -> search directories (compiler) -> add -> C:\Program Files\Boost
4. Error: [skipping 2 instantiation contexts, use -ftemplate-backtrace-limit=0 to disable]
    >> go to project head and right click to access >> build options
    >> choose project head (skip debug/release) >> compiler settings 
    >> other compiler option >> add: -ftemplate-backtrace-limit=0
5. supply command line arguments: project >> set program’s arguments
6. compiling:
    (a) General >> ……. [-std=gnu++17]
     (b) Warning >> …… [-Wall]
