ADD_EXECUTABLE( MSc_Math_Finance_2018.out main )
 TARGET_LINK_LIBRARIES( MSc_Math_Finance_2018.out MSc_Math_Finance_2018LIB )

ADD_EXECUTABLE( MSc_projects_4.out project-4-optimal-liquidation.cpp )
 TARGET_LINK_LIBRARIES( MSc_projects_4.out MSc_Math_Finance_2018LIB boost_system pthread)

ADD_EXECUTABLE( MSc_projects_5.out project-5-hedge-funds.cpp )
 TARGET_LINK_LIBRARIES( MSc_projects_5.out MSc_Math_Finance_2018LIB )
