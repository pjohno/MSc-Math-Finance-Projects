ADD_EXECUTABLE( MSc_Math_Finance_2018.out main )
 TARGET_LINK_LIBRARIES( MSc_Math_Finance_2018.out MSc_Math_Finance_2018LIB )

ADD_EXECUTABLE( MSc_projects_4.out project-4-optimal-liquidation-all-in-one.cpp )
 TARGET_LINK_LIBRARIES( MSc_projects_4.out MSc_Math_Finance_2018LIB)

ADD_EXECUTABLE( MSc_projects_5.out project-5-hedge-funds.cpp )
 TARGET_LINK_LIBRARIES( MSc_projects_5.out MSc_Math_Finance_2018LIB )

ADD_EXECUTABLE( MSc_projects_6.out project-6-revenue-management.cpp )
 TARGET_LINK_LIBRARIES( MSc_projects_6.out MSc_Math_Finance_2018LIB )

ADD_EXECUTABLE( MSc_projects_7.out project-7-exotic-options.cpp )
 TARGET_LINK_LIBRARIES( MSc_projects_7.out MSc_Math_Finance_2018LIB )
