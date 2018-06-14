
SET(MY_LIBRARY_SRCS  
MSc_Math_Finance_2018.cpp
MSc_projects_mVector.cpp
MSc_projects_gnuplot.cpp
MSc_projects_pde_solvers.cpp
)

ADD_LIBRARY(MSc_Math_Finance_2018LIB SHARED
${MY_LIBRARY_SRCS} 
)

TARGET_LINK_LIBRARIES( MSc_Math_Finance_2018LIB boost_system pthread)

ADD_LIBRARY(MSc_Math_Finance_2018LIB_static STATIC
${MY_LIBRARY_SRCS} 
)

TARGET_LINK_LIBRARIES( MSc_Math_Finance_2018LIB_static boost_system pthread)

