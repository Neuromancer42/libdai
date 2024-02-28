add_executable(createfg createfg.cpp)
target_link_libraries(createfg dai Boost::program_options)

add_executable(fg2dot fg2dot.cpp)
target_link_libraries(fg2dot dai)

add_executable(fginfo fginfo.cpp)
target_link_libraries(fginfo dai)

add_executable(uai2fg uai2fg.cpp)
target_link_libraries(uai2fg dai)


